## 51_sens_neutro.R — 敏感性分析:仅以中性粒细胞比例(1个协变量)校正
suppressPackageStartupMessages(library(limma))
# ---- project root: env var ARDS_EPI_ROOT, else walk up from script dir ----
resolve_root <- function() {
  env <- Sys.getenv("ARDS_EPI_ROOT", "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  start <- if (length(f)) dirname(normalizePath(f, mustWork = FALSE)) else getwd()
  d <- normalizePath(start, mustWork = FALSE)
  while (!identical(d, dirname(d))) {
    if (file.exists(file.path(d, "code", "01_download.R"))) return(d)
    d <- dirname(d)
  }
  d
}
ROOT <- resolve_root()
setwd(ROOT)
beta <- readRDS("data/processed/gse67530_beta.rds")
pd   <- readRDS("data/processed/gse67530_pdata.rds")
rmat <- read.csv("output/cell_estimates_GSE67530.csv", row.names = 1, check.names = FALSE)
pf <- function(x) trimws(sub("^[^:]*:\\s*","",as.character(x)))
grp <- suppressWarnings(as.numeric(pf(pd$characteristics_ch1_4)))
stopifnot(identical(rownames(rmat), colnames(beta)))
probes <- c("cg03341377","cg24310395","cg07830557","cg08418670")
m <- log2(beta[probes, ]/(1-beta[probes, ]))
neut <- rmat$Neutro
for (nm in c("ARDS_vs_healthy","ICU_vs_healthy")) {
  if (nm=="ARDS_vs_healthy") { cs <- grp %in% 1; ct <- is.na(grp) } else { cs <- grp %in% 0; ct <- is.na(grp) }
  idx <- cs | ct; g <- as.integer(cs[idx])
  X <- model.matrix(~ g + neut[idx])
  fit <- eBayes(lmFit(m[probes, idx], X))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  cat("\n==",nm,"(仅校正中性粒细胞比例) ==\n")
  print(data.frame(probe=rownames(tt), logFC=round(tt$logFC,3),
                   P=signif(tt$P.Value,3), adjP=signif(tt$adj.P.Val,3)), row.names=FALSE)
}
## 全基因组层面:仅校正中性粒细胞后 DMP 数
for (nm in c("ARDS_vs_healthy","ICU_vs_healthy")) {
  if (nm=="ARDS_vs_healthy") { cs <- grp %in% 1; ct <- is.na(grp) } else { cs <- grp %in% 0; ct <- is.na(grp) }
  idx <- cs | ct; g <- as.integer(cs[idx])
  X <- model.matrix(~ g + neut[idx])
  fit <- eBayes(lmFit(beta[, idx], X))
  tt <- topTable(fit, coef=2, number=Inf, sort.by="none")
  cat(sprintf("%s: 仅校正中性粒 -> FDR<0.05 的 DMP = %d\n", nm, sum(tt$adj.P.Val < 0.05)))
}
