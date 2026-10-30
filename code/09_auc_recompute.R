suppressPackageStartupMessages(library(pROC))
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
pd <- readRDS("data/processed/gse67530_pdata.rds")
pf <- function(x) trimws(sub("^[^:]*:\\s*","",as.character(x)))
grp <- suppressWarnings(as.numeric(pf(pd$characteristics_ch1_4)))
probes <- c("cg03341377","cg24310395","cg07830557","cg08418670")
ctr <- list(ARDS_vs_healthy=list(grp %in% 1, is.na(grp)),
            ICU_vs_healthy=list(grp %in% 0, is.na(grp)),
            ARDS_vs_ICU=list(grp %in% 1, grp %in% 0))
for (nm in names(ctr)) {
  idx <- ctr[[nm]][[1]] | ctr[[nm]][[2]]; g <- as.integer(ctr[[nm]][[1]][idx])
  for (p in probes) {
    r <- roc(g, as.numeric(beta[p, idx]), quiet=TRUE)   # direction auto
    ci <- as.numeric(ci.auc(r))
    cat(sprintf("%-16s %-12s AUC=%.3f (%.3f-%.3f)\n", nm, p, as.numeric(auc(r)), ci[1], ci[3]))
  }
}
