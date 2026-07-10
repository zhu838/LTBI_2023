"""latent_TB data preparation

Purpose
-------
Convert raw GBD export CSV/ZIP files under ./Data/RawGBD into cleaned, de-duplicated
CSV files under ./Data/CleanGBD that downstream scripts can consume directly.

This script splits data **by usage** (not by measure/metric/sex):

Outputs
-------
1) Global+Regional+SDI trend inputs (combined in a single table per metric):
    - global_regional_number.csv
    - global_regional_rate.csv

2) National-level files split by year:
    - national_by_year/national_YYYY.csv

Notes
-----
- Keeps `location_id` but renames it to `location` (matching the old pipeline).
- Drops other *_id columns and drops `cause_name` / `population_group_name`.
- Only keeps measures in TARGET_MEASURES if present (defaults: Prevalence).
"""

from __future__ import annotations

import argparse
import re
import sys
import zipfile
from pathlib import Path
from typing import Iterable

import pandas as pd


TARGET_MEASURES_DEFAULT = ["Prevalence"]


def _project_root() -> Path:
    # script is ./latent_TB/script/1_data_prepare.py
    return Path(__file__).resolve().parents[1]


def _clean_token(text: str) -> str:
    # mimic old naming convention: lowercase + underscores + strip parentheses content
    text = re.sub(r"\(.*$", "", str(text)).strip()
    text = text.lower()
    text = re.sub(r"[\s\-]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text


def _load_country_location_ids(root: Path) -> set[int]:
    """Load national-level (country/territory) location_ids from iso_code.csv.

    Returns an empty set if the file doesn't exist or can't be parsed.
    """

    iso_path = root / "Data" / "iso_code.csv"
    if not iso_path.exists():
        print(f"WARNING: {iso_path} not found; falling back to name heuristics for national split")
        return set()

    iso = pd.read_csv(iso_path, low_memory=False)
    if "location_id" not in iso.columns:
        print(f"WARNING: {iso_path} missing location_id column; falling back to name heuristics")
        return set()

    ids = pd.to_numeric(iso["location_id"], errors="coerce").dropna().astype(int)
    return set(ids.tolist())


def _iter_raw_tables(input_dir: Path) -> Iterable[pd.DataFrame]:
    # Accept both raw CSV and ZIP (containing same-named CSV).
    # If both are present for the same export, prefer the CSV and skip the ZIP to avoid double-loading.
    csv_basenames = {p.name for p in input_dir.iterdir() if p.is_file() and p.suffix.lower() == ".csv"}

    for p in sorted(input_dir.iterdir()):
        if not p.is_file():
            continue
        if p.suffix.lower() == ".csv":
            yield pd.read_csv(p, low_memory=False)
        elif p.suffix.lower() == ".zip":
            preferred = p.with_suffix(".csv").name
            if preferred in csv_basenames:
                continue
            with zipfile.ZipFile(p, "r") as zf:
                csv_candidates = [n for n in zf.namelist() if n.lower().endswith(".csv")]
                if not csv_candidates:
                    continue
                # Prefer a CSV whose basename matches the zip name; else take first.
                name = preferred if preferred in csv_candidates else csv_candidates[0]
                with zf.open(name) as f:
                    yield pd.read_csv(f, low_memory=False)


def _keep_and_rename_columns(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    if "location_id" in df.columns:
        df = df.rename(columns={"location_id": "location"})

    # Drop *_id columns (keep `location` which is not *_id anymore)
    drop_id_cols = [c for c in df.columns if c.endswith("_id")]
    if drop_id_cols:
        df = df.drop(columns=drop_id_cols)

    # Drop columns not used by downstream scripts
    for col in ("cause_name", "population_group_name"):
        if col in df.columns:
            df = df.drop(columns=[col])

    return df


def _deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    # Remove exact duplicate rows (common when users download multiple overlapping exports).
    before = len(df)
    df = df.drop_duplicates()
    after = len(df)
    print(f"Deduplicate: {before} -> {after} rows")
    return df


def _select_measures(df: pd.DataFrame, target_measures: list[str]) -> tuple[pd.DataFrame, list[str], list[str]]:
    if "measure_name" not in df.columns:
        raise ValueError("Missing required column: measure_name")

    available = sorted(df["measure_name"].dropna().unique().tolist())
    keep = [m for m in target_measures if m in available]
    missing = [m for m in target_measures if m not in available]
    if not keep:
        raise ValueError(f"None of target_measures exist. Available measures: {available}")
    if missing:
        print(f"Missing measures (skipped): {missing}")

    return df[df["measure_name"].isin(keep)].copy(), keep, missing


def _write_global_regional(df: pd.DataFrame, out_dir: Path, national_location_ids: set[int]) -> None:
    # Convenience file for downstream global/regional/SDI scripts:
    # keep only non-national locations (Global, WHO regions, SDI groups, etc.).
    if "location" not in df.columns:
        raise ValueError("Missing required column: location")
    if "metric_name" not in df.columns:
        raise ValueError("Missing required column: metric_name")

    location_ids = pd.to_numeric(df["location"], errors="coerce").astype("Int64")
    if national_location_ids:
        df_gr = df.loc[~location_ids.isin(list(national_location_ids))].copy()
    else:
        # Fallback: heuristic based on naming (less accurate but works without iso_code.csv)
        if "location_name" not in df.columns:
            raise ValueError("Missing required column: location_name")
        loc = df["location_name"].astype(str)
        is_aggregate = (loc == "Global") | loc.str.contains("SDI", case=False, regex=False) | loc.str.contains("Region", case=False, regex=False)
        df_gr = df.loc[is_aggregate].copy()

    df_gr[df_gr["metric_name"] == "Number"].to_csv(out_dir / "global_regional_number.csv", index=False)
    df_gr[df_gr["metric_name"] == "Rate"].to_csv(out_dir / "global_regional_rate.csv", index=False)


def _write_national_by_year(df: pd.DataFrame, out_dir: Path, national_location_ids: set[int]) -> None:
    if "year" not in df.columns:
        raise ValueError("Missing required column: year")
    if "location" not in df.columns:
        raise ValueError("Missing required column: location")
    if "location_name" not in df.columns:
        raise ValueError("Missing required column: location_name")

    location_ids = pd.to_numeric(df["location"], errors="coerce").astype("Int64")

    if national_location_ids:
        df_national = df.loc[location_ids.isin(list(national_location_ids))].copy()

        # If iso_code.csv misses some national ids, include them using heuristics.
        # This avoids accidental loss of locations like "Hong Kong" if not present in iso_code.
        locname = df["location_name"].astype(str)
        looks_aggregate = (locname == "Global") | locname.str.contains("SDI", case=False, regex=False) | locname.str.contains("Region", case=False, regex=False)
        missing_ids = set(location_ids.dropna().astype(int).unique().tolist()) - set(national_location_ids)
        if missing_ids:
            cand = df.loc[location_ids.isin(list(missing_ids)) & ~looks_aggregate]
            if not cand.empty:
                print(f"NOTE: Adding {cand['location_name'].nunique()} locations not found in iso_code.csv")
                df_national = pd.concat([df_national, cand], ignore_index=True)
    else:
        locname = df["location_name"].astype(str)
        looks_aggregate = (locname == "Global") | locname.str.contains("SDI", case=False, regex=False) | locname.str.contains("Region", case=False, regex=False)
        df_national = df.loc[~looks_aggregate].copy()

    out_year_dir = out_dir / "national_by_year"
    out_year_dir.mkdir(parents=True, exist_ok=True)

    # Remove stale outputs to avoid mixing different year ranges.
    for old in out_year_dir.glob("national_*.csv"):
        old.unlink()

    years = sorted(pd.unique(df_national["year"]))
    if years:
        print(
            f"Writing national-by-year files: {len(years)} years ({int(min(years))}-{int(max(years))}) -> {out_year_dir}"
        )
    else:
        print(f"Writing national-by-year files: 0 years -> {out_year_dir}")

    for y in years:
        part = df_national[df_national["year"] == y]
        # Avoid Excel scientific notation surprises: ensure year is int-like in filename
        y_int = int(y)
        part.to_csv(out_year_dir / f"national_{y_int}.csv", index=False)


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default=None, help="Raw data folder (default: ./Data/RawGBD under latent_TB)")
    parser.add_argument("--output-dir", default=None, help="Output folder (default: ./Data/CleanGBD under latent_TB)")
    parser.add_argument(
        "--measures",
        nargs="*",
        default=TARGET_MEASURES_DEFAULT,
        help="Measures to keep (default: Prevalence)",
    )
    args = parser.parse_args(argv)

    root = _project_root()
    if args.input_dir:
        input_dir = Path(args.input_dir)
    else:
        input_dir = root / "Data" / "RawGBD"
    output_dir = Path(args.output_dir) if args.output_dir else (root / "Data" / "CleanGBD")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_dir.exists():
        raise FileNotFoundError(f"Input dir not found: {input_dir}")

    tables = list(_iter_raw_tables(input_dir))
    if not tables:
        raise FileNotFoundError(f"No CSV/ZIP found in: {input_dir}")

    df = pd.concat(tables, ignore_index=True)
    print(f"Loaded raw rows: {len(df)}")

    df = _keep_and_rename_columns(df)
    df, kept, _missing = _select_measures(df, list(args.measures))
    df = _deduplicate(df)

    national_location_ids = _load_country_location_ids(root)

    # Ensure column order is friendly/consistent
    preferred_order = [
        "measure_name",
        "location",
        "location_name",
        "sex_name",
        "age_name",
        "metric_name",
        "year",
        "val",
        "lower",
        "upper",
    ]
    cols = [c for c in preferred_order if c in df.columns] + [c for c in df.columns if c not in preferred_order]
    df = df[cols]

    print(f"Writing purpose-based outputs to: {output_dir}")
    _write_global_regional(df, output_dir, national_location_ids)
    _write_national_by_year(df, output_dir, national_location_ids)

    print(f"Done. Measures kept: {kept}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main(sys.argv[1:]))
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        raise
