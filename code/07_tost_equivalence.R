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
idx <- grp %in% 1 | grp %in% 0; g <- as.integer(grp[idx] %in% 1)
for (p in probes) {
  y <- as.numeric(beta[p, idx]); tt <- t.test(y[g==1], y[g==0])
  d <- diff(rev(tt$estimate)); se <- diff(tt$conf.int)/(2*qt(0.975, tt$parameter))
  df <- tt$parameter
  for (m in c(0.05)) {
    p1 <- 1 - pt((d + m)/se, df); p2 <- pt((d - m)/se, df)
    cat(sprintf("%-12s d=%+.4f se=%.4f df=%.1f  p_lo=%.3g p_hi=%.3g TOST_p=%.3g\n", p, d, se, df, p1, p2, max(p1,p2)))
  }
}
