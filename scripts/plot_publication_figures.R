#!/usr/bin/env Rscript

# Publication-quality rendering for the CMS WCT paper.
#
# IMPORTANT:
#   This script is presentation-only. It does NOT refit omega, phase, amplitude,
#   backgrounds, or null models. It reads the frozen Python outputs already
#   written to CSV/JSON and renders them with ggplot2.
#
# Default inputs:
#   results/dimuon_certified_bootstrap500          (H1 discovery)
#   results/dimuon_replication_file2_frozen        (H2 frozen-omega replication)
#   results/dimuon_run2016g_frozen                 (G1 cross-period replication)
#   results/run2016g_file2_phase_locked             (G2 frozen omega/phase/sign holdout)
#
# Default behavior:
#   Rewrites the PNG filenames already referenced by the LaTeX manuscript and
#   writes matching PDF versions. It also creates summary figures in
#   results/publication_figures/.
#
# Usage:
#   Rscript scripts/plot_publication_figures.R
#   Rscript scripts/plot_publication_figures.R --root /path/to/wct-cms
#
# Required R packages:
#   ggplot2, jsonlite, scales
#
# Install once if needed:
#   install.packages(c("ggplot2", "jsonlite", "scales"))

required_packages <- c("ggplot2", "jsonlite", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing R packages: ", paste(missing_packages, collapse = ", "),
      "\nInstall with:\n  install.packages(c(",
      paste(sprintf('"%s"', missing_packages), collapse = ", "),
      "))"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(jsonlite)
  library(scales)
})

parse_root <- function(args) {
  root <- "."
  if (length(args) == 0) return(root)

  i <- 1
  while (i <= length(args)) {
    if (args[[i]] == "--root") {
      if (i == length(args)) stop("--root requires a path", call. = FALSE)
      root <- args[[i + 1]]
      i <- i + 2
      next
    }
    if (grepl("^--root=", args[[i]])) {
      root <- sub("^--root=", "", args[[i]])
      i <- i + 1
      next
    }
    stop(paste("Unknown argument:", args[[i]]), call. = FALSE)
  }
  root
}

repo_root <- normalizePath(
  parse_root(commandArgs(trailingOnly = TRUE)),
  winslash = "/",
  mustWork = TRUE
)
setwd(repo_root)

cat("CMS publication figure renderer\n")
cat("Repository root:", repo_root, "\n\n")

# -----------------------------------------------------------------------------
# Visual system
# -----------------------------------------------------------------------------

COL_DATA   <- "#202020"
COL_MODEL  <- "#2F6B9A"
COL_FROZEN <- "#B54A3A"
COL_BEST   <- "#6A6A6A"
COL_NULL   <- "#8AA6B8"
COL_BOOT   <- "#B7A57A"
COL_GRID   <- "#D9D9D9"

publication_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, color = "#444444", hjust = 0),
      plot.caption = element_text(size = base_size - 2, color = "#555555", hjust = 0),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(color = "#222222"),
      axis.line = element_line(linewidth = 0.45),
      axis.ticks = element_line(linewidth = 0.4),
      panel.grid.major.y = element_line(color = COL_GRID, linewidth = 0.3),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.key.width = grid::unit(1.3, "cm"),
      plot.margin = margin(8, 10, 8, 8)
    )
}

save_publication_plot <- function(plot, png_path, width = 8.8, height = 5.4) {
  dir.create(dirname(png_path), recursive = TRUE, showWarnings = FALSE)
  pdf_path <- sub("\\.png$", ".pdf", png_path, ignore.case = TRUE)
  if (identical(pdf_path, png_path)) pdf_path <- paste0(png_path, ".pdf")

  ggsave(
    filename = png_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 320,
    bg = "white"
  )
  ggsave(
    filename = pdf_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
  cat("  wrote:", png_path, "\n")
  cat("  wrote:", pdf_path, "\n")
}

read_json_required <- function(path) {
  if (!file.exists(path)) stop(paste("Missing required file:", path), call. = FALSE)
  fromJSON(path, simplifyVector = FALSE)
}

read_csv_required <- function(path) {
  if (!file.exists(path)) stop(paste("Missing required file:", path), call. = FALSE)
  read.csv(path, check.names = FALSE)
}

num <- function(x) as.numeric(unlist(x, use.names = FALSE)[1])

has_value <- function(x) {
  !is.null(x) && length(x) > 0 && !all(is.na(unlist(x, use.names = FALSE)))
}

fmt <- function(x, digits = 4) formatC(x, digits = digits, format = "f")

# -----------------------------------------------------------------------------
# Base-pipeline plots: H1 / H2 / G1
# -----------------------------------------------------------------------------

plot_omega_scan <- function(run_dir, title, subtitle = NULL) {
  scan_path <- file.path(run_dir, "omega_scan.csv")
  summary_path <- file.path(run_dir, "summary.json")
  if (!file.exists(scan_path) || !file.exists(summary_path)) {
    warning(paste("Skipping omega scan; missing outputs in", run_dir), call. = FALSE)
    return(invisible(NULL))
  }

  scan <- read_csv_required(scan_path)
  s <- read_json_required(summary_path)

  best_omega <- num(s$best_scan$omega)
  best_delta <- num(s$best_scan$delta_chi2)

  vlines <- data.frame(
    omega = best_omega,
    kind = "Exploratory maximum",
    stringsAsFactors = FALSE
  )

  if (has_value(s$frozen_scan)) {
    frozen_omega <- num(s$frozen_scan$omega)
    vlines <- rbind(
      vlines,
      data.frame(omega = frozen_omega, kind = "Frozen frequency")
    )
  }

  p <- ggplot(scan, aes(x = omega, y = delta_chi2)) +
    geom_line(linewidth = 0.8, color = COL_DATA) +
    geom_vline(
      data = vlines,
      aes(xintercept = omega, linetype = kind, color = kind),
      linewidth = 0.75,
      show.legend = TRUE
    ) +
    scale_color_manual(
      values = c(
        "Exploratory maximum" = COL_BEST,
        "Frozen frequency" = COL_FROZEN
      )
    ) +
    scale_linetype_manual(
      values = c(
        "Exploratory maximum" = "dashed",
        "Frozen frequency" = "solid"
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = expression(omega ~ "in" ~ ln(m[mumu] / m[0])),
      y = expression(Delta * chi^2),
      caption = paste0(
        "Exploratory maximum: omega = ", fmt(best_omega, 6),
        ", Delta chi^2 = ", fmt(best_delta, 2),
        if (has_value(s$frozen_scan))
          paste0(". Frozen test shown separately at omega = ", fmt(num(s$frozen_scan$omega), 6), ".")
        else "."
      )
    ) +
    publication_theme()

  save_publication_plot(p, file.path(run_dir, "omega_scan.png"))
}

plot_spectrum_background <- function(run_dir, title) {
  spectrum_path <- file.path(run_dir, "spectrum.csv")
  if (!file.exists(spectrum_path)) {
    warning(paste("Skipping spectrum plot; missing", spectrum_path), call. = FALSE)
    return(invisible(NULL))
  }

  d <- read_csv_required(spectrum_path)
  d <- d[is.finite(d$mass_center) & is.finite(d$count) & is.finite(d$background), ]

  p <- ggplot(d, aes(x = mass_center)) +
    geom_step(aes(y = count, color = "CMS dimuon data"), linewidth = 0.5, direction = "mid") +
    geom_line(aes(y = background, color = "Smooth continuum"), linewidth = 0.9) +
    scale_x_log10(
      breaks = c(2, 3, 5, 10, 20, 50, 100),
      labels = c("2", "3", "5", "10", "20", "50", "100")
    ) +
    scale_y_log10(labels = comma) +
    scale_color_manual(values = c("CMS dimuon data" = COL_DATA, "Smooth continuum" = COL_MODEL)) +
    labs(
      title = title,
      x = expression(m[mumu] ~ "[GeV]"),
      y = "Counts per bin",
      caption = "Rendering only: counts and fitted continuum are read from the frozen Python output spectrum.csv."
    ) +
    publication_theme()

  save_publication_plot(p, file.path(run_dir, "spectrum_background.png"), width = 8.8, height = 5.6)
}

plot_residuals <- function(run_dir, title) {
  spectrum_path <- file.path(run_dir, "spectrum.csv")
  summary_path <- file.path(run_dir, "summary.json")
  if (!file.exists(spectrum_path) || !file.exists(summary_path)) {
    warning(paste("Skipping residual plot; missing outputs in", run_dir), call. = FALSE)
    return(invisible(NULL))
  }

  d <- read_csv_required(spectrum_path)
  s <- read_json_required(summary_path)
  d <- d[
    d$analysis_mask == 1 &
      is.finite(d$mass_center) &
      is.finite(d$residual),
  ]

  p <- ggplot(d, aes(x = mass_center, y = residual)) +
    geom_hline(yintercept = 0, linewidth = 0.45, color = "#777777") +
    geom_point(size = 1.45, alpha = 0.78, color = COL_DATA) +
    scale_x_log10(
      breaks = c(2, 3, 5, 10, 20, 50, 100),
      labels = c("2", "3", "5", "10", "20", "50", "100")
    ) +
    labs(
      title = title,
      subtitle = paste0(
        "Residual RMS = ", fmt(num(s$residual_rms), 3),
        "; max |r| = ", fmt(num(s$residual_max_abs), 3)
      ),
      x = expression(m[mumu] ~ "[GeV]"),
      y = "Pearson residual",
      caption = "Only bins in the recorded analysis mask are shown."
    ) +
    publication_theme()

  save_publication_plot(p, file.path(run_dir, "residuals.png"), width = 8.8, height = 5.3)
}

plot_global_null <- function(run_dir, title) {
  null_path <- file.path(run_dir, "permutation_global_max.csv")
  summary_path <- file.path(run_dir, "summary.json")
  if (!file.exists(null_path) || !file.exists(summary_path)) {
    warning(paste("Skipping global null; missing outputs in", run_dir), call. = FALSE)
    return(invisible(NULL))
  }

  d <- read_csv_required(null_path)
  s <- read_json_required(summary_path)
  observed <- num(s$best_scan$delta_chi2)
  p_global <- if (has_value(s$global_permutation_p)) num(s$global_permutation_p) else NA_real_

  p <- ggplot(d, aes(x = max_delta_chi2)) +
    geom_histogram(bins = 40, fill = COL_NULL, color = "white", linewidth = 0.25) +
    geom_vline(xintercept = observed, color = COL_FROZEN, linewidth = 0.9) +
    labs(
      title = title,
      subtitle = if (is.finite(p_global)) paste0("Permutation global p = ", format(p_global, scientific = TRUE, digits = 3)) else NULL,
      x = expression("Permutation maximum " * Delta * chi^2),
      y = "Permutation count",
      caption = "The vertical line is the observed exploratory scan maximum."
    ) +
    publication_theme()

  save_publication_plot(p, file.path(run_dir, "global_null_distribution.png"), width = 8.6, height = 5.2)
}

# -----------------------------------------------------------------------------
# G2 prospective phase-locked holdout
# -----------------------------------------------------------------------------

plot_g2_phase_locked <- function(run_dir) {
  summary_path <- file.path(run_dir, "summary.json")
  spectrum_path <- file.path(run_dir, "spectrum.csv")
  if (!file.exists(summary_path) || !file.exists(spectrum_path)) {
    warning(paste("Skipping G2 phase-locked plot; missing outputs in", run_dir), call. = FALSE)
    return(invisible(NULL))
  }

  s <- read_json_required(summary_path)
  d <- read_csv_required(spectrum_path)

  omega <- num(s$omega)
  phase <- num(s$phase)
  amplitude <- num(s$observed$signed_amplitude)
  c0 <- num(s$observed$c0)
  delta <- num(s$observed$delta_chi2)
  local_p <- num(s$observed$local_p_one_sided_chi2_1dof)
  perm_p <- if (has_value(s$permutation_p)) num(s$permutation_p) else NA_real_
  boot_p <- if (has_value(s$parametric_bootstrap_p)) num(s$parametric_bootstrap_p) else NA_real_
  m0 <- if (!is.null(s$config) && has_value(s$config$m0)) num(s$config$m0) else 1.0

  d <- d[
    d$analysis_mask == 1 &
      is.finite(d$mass_center) &
      is.finite(d$residual) &
      d$mass_center > 0,
  ]
  if (nrow(d) < 3) stop("Too few valid G2 analysis bins", call. = FALSE)

  d$x <- log(d$mass_center / m0)
  d$centered_residual <- d$residual - mean(d$residual)

  curve <- data.frame(
    x = seq(min(d$x), max(d$x), length.out = 2000)
  )
  curve$locked_fit <- c0 + amplitude * cos(omega * curve$x - phase)

  annotation <- paste0(
    "Frozen omega = ", fmt(omega, 10), "\n",
    "Frozen phi = ", fmt(phase, 10), " rad\n",
    "A = ", fmt(amplitude, 4), "   Delta chi^2 = ", fmt(delta, 2), "\n",
    "local analytic p = ", format(local_p, scientific = TRUE, digits = 3), "\n",
    if (is.finite(perm_p)) paste0("permutation floor = ", format(perm_p, scientific = TRUE, digits = 3), "\n") else "",
    if (is.finite(boot_p)) paste0("refit-bootstrap floor = ", format(boot_p, scientific = TRUE, digits = 3)) else ""
  )

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "#777777", linewidth = 0.45) +
    geom_point(
      data = d,
      aes(x = x, y = centered_residual, color = "G2 centered Pearson residuals"),
      size = 1.55,
      alpha = 0.70
    ) +
    geom_line(
      data = curve,
      aes(x = x, y = locked_fit, color = "Frozen phase-locked waveform"),
      linewidth = 1.05
    ) +
    scale_color_manual(
      values = c(
        "G2 centered Pearson residuals" = COL_DATA,
        "Frozen phase-locked waveform" = COL_FROZEN
      )
    ) +
    annotate(
      "label",
      x = min(d$x) + 0.02 * diff(range(d$x)),
      y = max(d$centered_residual, na.rm = TRUE) * 0.96,
      label = annotation,
      hjust = 0,
      vjust = 1,
      size = 3.15,
      label.size = 0.25,
      fill = "white"
    ) +
    labs(
      title = "CMS Run2016G G2 prospective phase-locked holdout",
      subtitle = "Frequency, phase, sign, file rule, and analysis configuration were frozen before target inspection",
      x = expression(x == ln(m[mumu] / "1 GeV")),
      y = "Centered Pearson residual",
      caption = paste0(
        "Visualization only; no parameters are refit in R. ",
        "The local analytic p-value is conditional on the implemented fixed-waveform null; ",
        "the permutation/bootstrap values are finite Monte Carlo resolution floors."
      )
    ) +
    publication_theme()

  save_publication_plot(
    p,
    file.path(run_dir, "phase_locked_fit.png"),
    width = 9.2,
    height = 5.9
  )
}

plot_g2_nulls <- function(run_dir) {
  summary_path <- file.path(run_dir, "summary.json")
  perm_path <- file.path(run_dir, "permutation_locked.csv")
  boot_path <- file.path(run_dir, "bootstrap_locked.csv")
  if (!file.exists(summary_path) || !file.exists(perm_path) || !file.exists(boot_path)) {
    warning(paste("Skipping G2 null figure; missing outputs in", run_dir), call. = FALSE)
    return(invisible(NULL))
  }

  s <- read_json_required(summary_path)
  observed <- num(s$observed$delta_chi2)
  perm <- read_csv_required(perm_path)
  boot <- read_csv_required(boot_path)

  dp <- data.frame(delta_chi2 = perm[[1]], null = "Residual permutation")
  db <- data.frame(delta_chi2 = boot[[1]], null = "Poisson + background refit")
  d <- rbind(dp, db)

  p <- ggplot(d, aes(x = delta_chi2, fill = null)) +
    geom_histogram(bins = 38, color = "white", linewidth = 0.2, alpha = 0.82) +
    geom_vline(xintercept = observed, color = COL_FROZEN, linewidth = 0.95) +
    facet_wrap(~null, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c(
      "Residual permutation" = COL_NULL,
      "Poisson + background refit" = COL_BOOT
    )) +
    labs(
      title = "G2 frozen-waveform empirical null ensembles",
      subtitle = paste0("Observed locked Delta chi^2 = ", fmt(observed, 2)),
      x = expression(Delta * chi^2 ~ "for frozen waveform"),
      y = "Pseudoexperiment count",
      caption = "Zero exceedances in the recorded 1000 permutation and 500 refit-bootstrap ensembles; these are resolution floors, not resolved tail probabilities."
    ) +
    publication_theme() +
    theme(legend.position = "none")

  save_publication_plot(
    p,
    file.path(run_dir, "phase_locked_nulls.png"),
    width = 8.8,
    height = 7.0
  )
}

# -----------------------------------------------------------------------------
# Cross-stage summary figures (read recorded Python results only)
# -----------------------------------------------------------------------------

collect_stage_summary <- function(h1_dir, h2_dir, g1_dir, g2_dir) {
  h1 <- read_json_required(file.path(h1_dir, "summary.json"))
  h2 <- read_json_required(file.path(h2_dir, "summary.json"))
  g1 <- read_json_required(file.path(g1_dir, "summary.json"))
  g2 <- read_json_required(file.path(g2_dir, "summary.json"))

  out <- data.frame(
    sample = c("H1 discovery", "H2 frozen omega", "G1 frozen omega", "G2 frozen omega+phi+sign"),
    stage = c("Discovery", "Replication", "Replication", "Prospective holdout"),
    amplitude = c(
      num(h1$best_scan$amplitude),
      num(h2$frozen_scan$amplitude),
      num(g1$frozen_scan$amplitude),
      num(g2$observed$signed_amplitude)
    ),
    phase = c(
      num(h1$best_scan$phase),
      num(h2$frozen_scan$phase),
      num(g1$frozen_scan$phase),
      num(g2$phase)
    ),
    delta_chi2 = c(
      num(h1$best_scan$delta_chi2),
      num(h2$frozen_scan$delta_chi2),
      num(g1$frozen_scan$delta_chi2),
      num(g2$observed$delta_chi2)
    ),
    stringsAsFactors = FALSE
  )

  out$sample <- factor(out$sample, levels = out$sample)
  out
}

plot_replication_summary <- function(d, out_dir) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(d, file.path(out_dir, "replication_stage_summary.csv"), row.names = FALSE)

  p_amp <- ggplot(d, aes(x = sample, y = amplitude, shape = stage)) +
    geom_hline(yintercept = 0, color = "#888888", linewidth = 0.35) +
    geom_point(size = 3.4, color = COL_DATA) +
    geom_text(
      aes(label = sprintf("%.4f", amplitude)),
      nudge_y = 0.035,
      size = 3.3,
      color = "#333333"
    ) +
    scale_shape_manual(values = c("Discovery" = 16, "Replication" = 17, "Prospective holdout" = 18)) +
    coord_cartesian(ylim = c(0.68, 1.05)) +
    labs(
      title = "CMS fixed-component amplitude across the replication sequence",
      subtitle = "H1 is exploratory; H2/G1 use frozen omega; G2 uses frozen omega, phase, and positive sign",
      x = NULL,
      y = "Residual amplitude",
      caption = "No uncertainty bars are invented here; values are the recorded point estimates from the Python summaries."
    ) +
    publication_theme() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))

  save_publication_plot(
    p_amp,
    file.path(out_dir, "replication_amplitudes.png"),
    width = 8.8,
    height = 5.5
  )

  p_phase <- ggplot(d, aes(x = sample, y = phase, shape = stage)) +
    geom_hline(yintercept = 0, color = "#888888", linewidth = 0.35) +
    geom_point(size = 3.4, color = COL_DATA) +
    geom_text(
      aes(label = sprintf("%.3f", phase)),
      nudge_y = 0.025,
      size = 3.2,
      color = "#333333"
    ) +
    scale_shape_manual(values = c("Discovery" = 16, "Replication" = 17, "Prospective holdout" = 18)) +
    labs(
      title = "CMS phase across the replication sequence",
      subtitle = "G2 phase is the predeclared circular-mean prediction from H2 and G1, not a fitted G2 phase",
      x = NULL,
      y = "Phase [rad]",
      caption = "All phases use the same m0 = 1 GeV log-mass convention."
    ) +
    publication_theme() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))

  save_publication_plot(
    p_phase,
    file.path(out_dir, "replication_phases.png"),
    width = 8.8,
    height = 5.5
  )

  p_delta <- ggplot(d, aes(x = sample, y = delta_chi2, shape = stage)) +
    geom_point(size = 3.4, color = COL_DATA) +
    geom_text(
      aes(label = sprintf("%.2f", delta_chi2)),
      nudge_y = 5,
      size = 3.2,
      color = "#333333"
    ) +
    scale_shape_manual(values = c("Discovery" = 16, "Replication" = 17, "Prospective holdout" = 18)) +
    coord_cartesian(ylim = c(60, max(d$delta_chi2) + 15)) +
    labs(
      title = "Observed log-periodic improvement statistic across stages",
      x = NULL,
      y = expression(Delta * chi^2),
      caption = "Statistics are not combined across samples; each point is its recorded stage-specific statistic."
    ) +
    publication_theme() +
    theme(axis.text.x = element_text(angle = 18, hjust = 1))

  save_publication_plot(
    p_delta,
    file.path(out_dir, "replication_delta_chi2.png"),
    width = 8.8,
    height = 5.5
  )
}

# -----------------------------------------------------------------------------
# Execute
# -----------------------------------------------------------------------------

H1_DIR <- "results/dimuon_certified_bootstrap500"
H2_DIR <- "results/dimuon_replication_file2_frozen"
G1_DIR <- "results/dimuon_run2016g_frozen"
G2_DIR <- "results/run2016g_file2_phase_locked"
PUB_DIR <- "results/publication_figures"

cat("Rendering H1...\n")
plot_omega_scan(
  H1_DIR,
  "Run2016H H1 certified discovery frequency scan",
  "Exploratory discovery stage"
)

cat("Rendering H2...\n")
plot_omega_scan(
  H2_DIR,
  "Run2016H H2 frozen-frequency replication",
  "Primary inference is at the previously frozen frequency"
)

cat("Rendering G1...\n")
plot_omega_scan(
  G1_DIR,
  "Run2016G G1 cross-period frozen-frequency replication",
  "Primary inference is at the preregistered frozen frequency"
)
plot_spectrum_background(
  G1_DIR,
  "Run2016G opposite-sign dimuon spectrum and smooth continuum"
)
plot_residuals(
  G1_DIR,
  "Run2016G Pearson residual field"
)
plot_global_null(
  G1_DIR,
  "Run2016G exploratory global permutation null"
)

cat("Rendering G2...\n")
plot_g2_phase_locked(G2_DIR)
plot_g2_nulls(G2_DIR)

cat("Rendering cross-stage summaries...\n")
all_summaries_exist <- all(file.exists(c(
  file.path(H1_DIR, "summary.json"),
  file.path(H2_DIR, "summary.json"),
  file.path(G1_DIR, "summary.json"),
  file.path(G2_DIR, "summary.json")
)))

if (all_summaries_exist) {
  stages <- collect_stage_summary(H1_DIR, H2_DIR, G1_DIR, G2_DIR)
  plot_replication_summary(stages, PUB_DIR)
} else {
  warning("Skipping cross-stage summaries because one or more summary.json files are missing.", call. = FALSE)
}

cat("\nDone.\n")
cat("R was used only to render recorded Python outputs; no analysis parameters were refit.\n")
