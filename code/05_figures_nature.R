# =====================================================================
# 05_figures_nature.R — regenerate all manuscript figures in Nature style
#
# Nature figure guidelines applied:
#   * Arial (sans-serif) throughout, 7-9 pt at final size
#   * open axes: left + bottom spines only, ticks pointing outward
#   * no background grid, no box
#   * single-column width 89 mm, double-column 183 mm, 600 dpi
#   * Okabe-Ito colour-blind-safe palette
#   * panel-free standalone figures (each = one manuscript figure)
#
# Inputs  : output/ewas_dmps_*.csv, output/cell_estimates_GSE67530.csv,
#           output/gse199591_*.csv, data/processed/gse67530_pdata.rds
# Outputs : output/fig_volcano_*.png, output/fig_cellprops_*.png
#           (same filenames as before -> drop-in replacement)
# =====================================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
  library(data.table)
  library(tidyr)
})

# ---- project root: env var ARDS_EPI_ROOT, else walk up from script dir ----
resolve_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "R", "01_download.R")) || file.exists(file.path(d, "code", "01_download.R"))) return(d)
    d <- dirname(d)
  }
  d
}
# ---- project root: env var ARDS_EPI_ROOT, else walk up from script dir ----
resolve_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "R", "01_download.R")) || file.exists(file.path(d, "code", "01_download.R"))) return(d)
    d <- dirname(d)
  }
  d
}
BASE <- resolve_root()
OUT    <- file.path(BASE, "output")
SIG    <- c("cg03341377", "cg24310395", "cg07830557", "cg08418670")

FAM    <- "Arial"          # falls back to system sans if unavailable
COL_OK <- "#999999"        # grey   : non-significant
COL_SG <- "#0072B2"        # blue   : FDR significant
COL_HS <- "#D55E00"        # vermillion : signature CpGs

# ---------------------------------------------------------------
# Nature-style theme
# ---------------------------------------------------------------
theme_nature <- function(base_size = 8) {
  theme_bw(base_size = base_size, base_family = FAM) +
    theme(
      panel.border     = element_blank(),
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = 0.35),
      axis.ticks       = element_line(colour = "black", linewidth = 0.35),
      axis.ticks.length = unit(1.6, "pt"),
      axis.text        = element_text(colour = "black", size = base_size - 1),
      axis.title       = element_text(colour = "black", size = base_size),
      plot.title       = element_text(colour = "black", size = base_size,
                                      hjust = 0, face = "plain"),
      plot.margin      = margin(4, 6, 4, 4, "pt"),
      legend.key       = element_blank(),
      legend.background = element_rect(fill = "white", colour = NA),
      legend.box.background = element_rect(fill = "white", colour = NA),
      legend.text      = element_text(colour = "black", size = base_size - 1),
      legend.title     = element_text(colour = "black", size = base_size - 1),
      legend.margin    = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text       = element_text(colour = "black", size = base_size - 1)
    )
}

mm <- function(w, h) list(width = w, height = h)

# ---------------------------------------------------------------
# Volcano plot factory
# ---------------------------------------------------------------
volcano <- function(dmp_file, out_file, title,
                    width_mm = 89, height_mm = 82, ymax_cap = 30) {
  d <- fread(file.path(OUT, dmp_file))
  d[, y := -log10(adj.P.Val)]
  capped_n <- sum(d$y > ymax_cap)
  if (capped_n > 0) d[y > ymax_cap, y := ymax_cap]
  d[, sig := adj.P.Val < 0.05]
  d[, grp := fcase(
    probe %in% SIG, "Signature CpG",
    sig,            "FDR < 0.05",
    default =        "NS"
  )]
  d[, grp := factor(grp, levels = c("NS", "FDR < 0.05", "Signature CpG"))]
  d[, pr := probe %in% SIG]

  ymax_axis <- min(max(d$y) + 1, ymax_cap + 1.5)

  p <- ggplot(d[order(grp)], aes(logFC, y)) +
    geom_point(aes(colour = grp, size = grp, alpha = grp),
               shape = 16, stroke = 0) +
    scale_colour_manual(
      values = c("NS" = COL_OK, "FDR < 0.05" = COL_SG,
                 "Signature CpG" = COL_HS),
      breaks = c("FDR < 0.05", "Signature CpG"),
      guide  = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    scale_size_manual(values = c("NS" = 0.25, "FDR < 0.05" = 0.45,
                                 "Signature CpG" = 1.6), guide = "none") +
    scale_alpha_manual(values = c("NS" = 0.35, "FDR < 0.05" = 0.7,
                                  "Signature CpG" = 1), guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "black", linewidth = 0.25) +
    labs(title = title,
         x = expression("Methylation difference (" * Delta * "M-value)"),
         y = expression(-log[10] * " (FDR-adjusted " * italic(P) * ")"),
         colour = NULL) +
    annotate("text", x = 0, y = -log10(0.05) + ymax_axis * 0.025,
             label = "FDR = 0.05", size = 2.1, family = FAM, hjust = 0.5,
             vjust = -0.3) +
    coord_cartesian(xlim = c(-2.4, 2.4),
                    ylim = c(0, ymax_axis), clip = "off") +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    theme_nature(8) +
    theme(legend.position = c(0.02, 0.98),
          legend.justification = c(0, 1),
          legend.direction = "vertical")

  # label signature CpGs present in the dataset (repel avoids overlap)
  sub <- d[pr == TRUE]
  if (nrow(sub) > 0) {
    p <- p + geom_text_repel(
      data = sub, aes(label = probe), size = 2.1, family = FAM,
      colour = "black", segment.colour = "black", segment.size = 0.2,
      min.segment.length = 0, box.padding = 0.28, point.padding = 0.2,
      max.overlaps = 20, seed = 42)
  }
  if (capped_n > 0) {
    p <- p + annotate("text", x = -Inf, y = ymax_cap,
                      label = sprintf("%d probes beyond axis", capped_n),
                      size = 2.0, family = FAM, hjust = -0.05,
                      vjust = 1.4, fontface = "italic")
  }

  ggsave(file.path(OUT, out_file), p,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 600, bg = "white", limitsize = FALSE)
  cat("saved:", out_file, " (capped:", capped_n, ")\n")
}

# ---------------------------------------------------------------
# Cell-proportion boxplot factory
# ---------------------------------------------------------------
cellbox <- function(props, pdata_group, out_file, title,
                    width_mm = 183, height_mm = 70,
                    group_levels, group_labels, group_cols,
                    base_size = 8) {
  ct <- c("B", "NK", "CD4T", "CD8T", "Mono", "Neutro", "Eosino")
  props <- as.data.frame(props)
  props <- props[, c("geo_accession", ct), drop = FALSE]
  gdf <- data.frame(geo_accession = names(pdata_group),
                    group = unname(pdata_group), stringsAsFactors = FALSE)
  m <- merge(props, gdf, by = "geo_accession")
  stopifnot(nrow(m) == nrow(gdf))
  m$group <- factor(m$group, levels = group_levels, labels = group_labels)
  long <- pivot_longer(m, all_of(ct), names_to = "cell", values_to = "prop")
  long$cell <- factor(long$cell, levels = ct)

  p <- ggplot(long, aes(cell, prop, fill = group)) +
    geom_boxplot(outlier.size = 0.35, outlier.alpha = 0.5,
                 linewidth = 0.3, outlier.shape = 21,
                 outlier.stroke = 0.2, width = 0.75,
                 colour = "black") +
    scale_fill_manual(values = group_cols, name = NULL) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       expand = expansion(mult = c(0, 0.03))) +
    labs(title = title,
         x = NULL,
         y = "Estimated cell proportion") +
    theme_nature(base_size) +
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
          legend.position = "top",
          legend.justification = "left",
          legend.direction = "horizontal")

  ggsave(file.path(OUT, out_file), p,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 600, bg = "white", limitsize = FALSE)
  cat("saved:", out_file, "\n")
}

# ---------------------------------------------------------------
# 1-3. GSE67530 volcanoes
# ---------------------------------------------------------------
volcano("ewas_dmps_primary_paper.csv",  "fig_volcano_primary_paper.png",
        "ARDS vs healthy")
volcano("ewas_dmps_icu_vs_healthy.csv", "fig_volcano_icu_vs_healthy.png",
        "ICU controls vs healthy")
volcano("ewas_dmps_ards_vs_icu.csv",    "fig_volcano_ards_vs_icu.png",
        "ARDS vs ICU controls")

# ---------------------------------------------------------------
# 4. GSE67530 cell proportions
# ---------------------------------------------------------------
pdata <- readRDS(file.path(BASE, "data", "processed", "gse67530_pdata.rds"))
if (!"geo_accession" %in% names(pdata)) pdata$geo_accession <- rownames(pdata)
pdata$geo_accession <- as.character(pdata$geo_accession)

# rebuild the group label the same way 03_deconvolution.R does
ch <- names(pdata)[grepl("characteristics", names(pdata), ignore.case = TRUE)]
kv <- setNames(vector("list", length(ch)), ch)
for (cc in ch) {
  v   <- as.character(pdata[[cc]])
  k   <- trimws(sub(":.*$", "", v))
  val <- trimws(sub("^[^:]*:\\s*", "", v))
  for (kk in unique(k)) kv[[kk]] <- ifelse(k == kk, val, NA_character_)
}
ards_raw <- trimws(kv$ards)
grp <- ifelse(ards_raw == "1", "ARDS",
       ifelse(ards_raw == "0", "ICU_control",
       ifelse(ards_raw %in% c("NA", "", "na"), "healthy", NA_character_)))
names(grp) <- pdata$geo_accession
cat("GSE67530 groups:\n"); print(table(grp, useNA = "ifany"))

props67530 <- fread(file.path(OUT, "cell_estimates_GSE67530.csv"))
props67530$geo_accession <- as.character(props67530$geo_accession)
props67530 <- props67530[props67530$geo_accession %in% names(grp), ]

cellbox(props67530, grp, "fig_cellprops_GSE67530.png",
        "GSE67530 (450K): EpiDISH-estimated leukocyte proportions",
        width_mm = 183, height_mm = 72,
        group_levels = c("healthy", "ICU_control", "ARDS"),
        group_labels = c("Healthy (n=30)", "ICU controls (n=75)",
                         "ARDS (n=39)"),
        group_cols   = c(COL_OK, COL_SG, COL_HS))

# ---------------------------------------------------------------
# 5. GSE199591 volcano (Severe vs Mild)
# ---------------------------------------------------------------
volcano("gse199591_dmps_severe_vs_mild.csv", "fig_gse199591_volcano.png",
        "COVID-19 severe vs mild (GSE199591)", ymax_cap = 20)

# ---------------------------------------------------------------
# 6. GSE199591 cell proportions
# ---------------------------------------------------------------
pd2 <- fread(file.path(OUT, "gse199591_pdata.csv"))
grp2 <- setNames(ifelse(pd2$severe == 1, "Severe", "Mild"),
                 pd2$geo_accession)
cat("GSE199591 groups:\n"); print(table(grp2, useNA = "ifany"))

props199 <- fread(file.path(OUT, "gse199591_cellprops.csv"))
props199$geo_accession <- as.character(props199$geo_accession)
props199 <- props199[props199$geo_accession %in% names(grp2), ]

cellbox(props199, grp2, "fig_gse199591_cellprops.png",
        "GSE199591 (EPIC): EpiDISH-estimated leukocyte proportions",
        width_mm = 140, height_mm = 72,
        group_levels = c("Mild", "Severe"),
        group_labels = c("Mild (n=123)", "Severe (n=64)"),
        group_cols   = c(COL_SG, COL_HS))

cat("ALL FIGURES DONE\n")
