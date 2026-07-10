## fig axis
scientific_10 <- function(x) {
  ifelse(x == 0, 0, parse(text = gsub("[+]", "", gsub("e", "%*%10^", scales::scientific_format()(x)))))
}

## get AAPC from jp_model
get_aapc <- function(jp_model) {
  data <- jp_model$aapc |>
    dplyr::mutate(
      dplyr::across(c(aapc, aapc_c_i_low, aapc_c_i_high), ~formatC(., format = "f", digits = 2)),
      p_value = as.numeric(p_value),
      p_value_label = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      legend = paste0(
        start_obs, "~", end_obs, "\n",
        aapc, "(", aapc_c_i_low, "~", aapc_c_i_high, ")", p_value_label
      ),
      Year = paste(start_obs, end_obs, sep = "~"),
      Value = paste0(aapc, " (", aapc_c_i_low, "~", aapc_c_i_high, ")")
    )

  return(data)
}

## get segment-specific APC from a joinpoint model
get_apc <- function(jp_model) {
  jp_model$apc |>
    dplyr::mutate(
      dplyr::across(c(apc, apc_95_lcl, apc_95_ucl), ~formatC(., format = "f", digits = 2)),
      p_value = as.numeric(p_value),
      p_value_label = dplyr::case_when(
        p_value < 0.001 ~ "<0.001",
        TRUE ~ formatC(p_value, format = "f", digits = 3)
      ),
      Interval = paste(segment_start, segment_end, sep = " to "),
      `Estimate (95% CI)` = paste0(apc, " (", apc_95_lcl, " to ", apc_95_ucl, ")")
    ) |>
    dplyr::select(Interval, `Estimate (95% CI)`, p_value_label)
}

## visual apc
plot_apc <- function(jp_model, data, use_scientific_10 = TRUE, y_divisor = 1) {
  if (!is.numeric(y_divisor) || length(y_divisor) != 1 || is.na(y_divisor) || y_divisor <= 0) {
    stop("y_divisor must be a single positive number")
  }

  if (y_divisor != 1) {
    data <- data |>
      dplyr::mutate(dplyr::across(c(val, lower, upper), ~ .x / y_divisor))
  }

  df_jp_apc <- jp_model$apc |>
    dplyr::mutate(
      dplyr::across(c(apc, apc_95_lcl, apc_95_ucl), ~round(., 2)),
      p_value = as.numeric(p_value),
      p_value_label = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01 ~ "**",
        p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      legend = paste0(
        segment_start, "~", segment_end, "\n",
        apc, "(", apc_95_lcl, "~", apc_95_ucl, ")", p_value_label
      )
    )

  # get breaks of y axis
  breaks <- pretty(c(data$val, data$lower, data$upper))

  # Use the bottom 15% of the y-range to display APC rectangles.
  y_min <- min(breaks, na.rm = TRUE)
  y_max <- max(breaks, na.rm = TRUE)
  y_apc_top <- y_min + 0.15 * (y_max - y_min)

  # APC fill colors: pick from a continuous palette by APC value,
  # but keep a discrete legend (one entry per segment) in each panel.
  palette_base <- as.character(paletteer::paletteer_d("MetBrewer::Paquin", direction = -1))
  palette_cont <- grDevices::colorRampPalette(palette_base)(256)
  apc_vals <- df_jp_apc$apc
  apc_limits <- c(-3, 3)
  apc_vals_clamped <- pmin(apc_limits[[2]], pmax(apc_limits[[1]], apc_vals))
  idx <- floor(scales::rescale(apc_vals_clamped, to = c(1, 256), from = apc_limits))
  idx <- pmin(256, pmax(1, idx))
  colors <- stats::setNames(palette_cont[idx], df_jp_apc$legend)

  fig <- ggplot2::ggplot(data) +
    ggplot2::geom_vline(
      data = df_jp_apc,
      mapping = ggplot2::aes(xintercept = segment_end),
      alpha = 0.5,
      color = "grey50"
    ) +
    ggplot2::geom_rect(
      data = df_jp_apc,
      ggplot2::aes(
        xmin = segment_start,
        xmax = segment_end,
        fill = legend
      ),
      ymin = y_min,
      ymax = y_apc_top,
      alpha = 0.5
    ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, x = year, y = val),
      fill = "#00798CFF",
      alpha = 0.35
    ) +
    ggplot2::geom_line(ggplot2::aes(x = year, y = val), color = "#00798CFF", size = 0.7) +
    ggplot2::geom_point(ggplot2::aes(x = year, y = val), color = "#00798CFF", size = 1.2) +
    ggplot2::scale_x_continuous(
      limits = range(data$year),
      breaks = scales::pretty_breaks(n = 7),
      expand = ggplot2::expansion(add = c(0, 1))
    ) +
    ggplot2::scale_fill_manual(values = colors, breaks = df_jp_apc$legend) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title.position = "plot",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.justification.bottom = "right",
      legend.title.position = "top",
      legend.key.spacing.y = grid::unit(0.35, "cm")
    )

  if (use_scientific_10) {
    fig <- fig +
      ggplot2::scale_y_continuous(
        limits = range(breaks),
        breaks = breaks,
        labels = scientific_10,
        expand = ggplot2::expansion(mult = c(0, 0))
      )
  } else {
    fig <- fig +
      ggplot2::scale_y_continuous(
        limits = range(breaks),
        breaks = breaks,
        expand = ggplot2::expansion(mult = c(0, 0))
      )
  }

  return(fig)
}

# visualize val
plot_val <- function(data, measure, filter_col, filter_val, ylab) {
  data_filtered <- data |>
    dplyr::filter(measure_name == measure)

  breaks <- pretty(c(0, range(data_filtered$upper)), n = 5)

  data_filtered <- data_filtered |>
    dplyr::filter(!!rlang::sym(filter_col) == filter_val)

  p <- ggplot2::ggplot(data_filtered, ggplot2::aes(x = year, y = val)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.3) +
    ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 7), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = breaks, limits = range(breaks), expand = c(0, 0)) +
    ggplot2::labs(title = filter_val, x = NULL, y = ylab) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title.position = "plot",
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  return(p)
}

write_markdown_table <- function(data, file_path) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }

  df_chr <- data |>
    dplyr::mutate(dplyr::across(dplyr::everything(), ~ {
      x <- as.character(.x)
      x[is.na(x)] <- ""
      gsub("\\|", "\\\\|", x)
    }))

  header <- paste0("|", paste(names(df_chr), collapse = "|"), "|")
  separator <- paste0("|", paste(rep(":---", ncol(df_chr)), collapse = "|"), "|")

  rows <- apply(df_chr, 1, function(row) {
    paste0("|", paste(row, collapse = "|"), "|")
  })

  writeLines(c(header, separator, rows), con = file_path)
}
