## 46_equiv_v2.R — 未校正Δβ的95%CI + AUC尺度等价性(单侧检验 vs 0.70/0.65)
suppressPackageStartupMessages({library(pROC)})
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
# lock the repo dir absolutely BEFORE setwd() changes the cwd
ca <- commandArgs(FALSE)
cf <- sub("^--file=", "", ca[grep("^--file=", ca)])
REPO_DIR <- if (length(cf)) dirname(dirname(normalizePath(cf[1], mustWork = FALSE))) else ROOT
setwd(ROOT)
# output/ (pipeline) with results/ (repo snapshot) fallback
resfile <- function(p) {
  cands <- c(file.path(ROOT, "output", p),
             file.path(ROOT, "results", p),
             file.path(REPO_DIR, "results", p))
  for (x in cands) if (file.exists(x)) return(x)
  cands[1]
}

beta <- readRDS("data/processed/gse67530_beta.rds")
pd   <- readRDS("data/processed/gse67530_pdata.rds")
pf <- function(x) trimws(sub("^[^:]*:\\s*", "", as.character(x)))
grp <- suppressWarnings(as.numeric(pf(pd$characteristics_ch1_4)))
probes <- c("cg03341377", "cg24310395", "cg07830557", "cg08418670")

contr <- list(
  ARDS_vs_healthy = list(case = grp %in% 1,  ctrl = is.na(grp)),
  ICU_vs_healthy  = list(case = grp %in% 0,  ctrl = is.na(grp)),
  ARDS_vs_ICU     = list(case = grp %in% 1,  ctrl = grp %in% 0)
)
rows <- list()
for (nm in names(contr)) {
  cs <- contr[[nm]]$case; ct <- contr[[nm]]$ctrl
  idx <- cs | ct; g <- as.integer(cs[idx])
  for (p in probes) {
    y  <- as.numeric(beta[p, idx])
    tt <- t.test(y[g == 1], y[g == 0])
    rows[[length(rows) + 1]] <- data.frame(contrast = nm, probe = p,
      n_case = sum(g == 1), n_ctrl = sum(g == 0),
      delta_beta = mean(y[g == 1]) - mean(y[g == 0]),
      lo95 = tt$conf.int[1], hi95 = tt$conf.int[2], p = tt$p.value)
  }
}
res <- do.call(rbind, rows)
cat("== 未校正 Δβ(β 尺度) 与 95%CI ==\n")
print(res, digits = 4, row.names = FALSE)
write.csv(res, "output/equiv_delta_ci_unadj.csv", row.names = FALSE)

## 比值:ARDS vs ICU 的点估计与上界 相对 ARDS vs healthy
wide <- reshape(res[, c("contrast","probe","delta_beta","lo95","hi95")],
                idvar = "probe", timevar = "contrast", direction = "wide")
cat("\n== 比值 ARDS_vs_ICU / ARDS_vs_healthy(绝对值) ==\n")
print(data.frame(probe = wide$probe,
  ratio_point = abs(wide$delta_beta.ARDS_vs_ICU) / abs(wide$delta_beta.ARDS_vs_healthy),
  ratio_hi95  = abs(wide$hi95.ARDS_vs_ICU) / abs(wide$delta_beta.ARDS_vs_healthy)),
  digits = 3, row.names = FALSE)

## AUC 尺度:联合评分 0.561,DeLong SE,单侧 H0: AUC >= 0.70 / 0.65 / 0.60
sc <- function(idx) {
  b <- beta[probes, idx, drop = FALSE]
  (as.numeric(b[1, ]) + as.numeric(b[2, ]) + (1 - as.numeric(b[3, ])) +
     (1 - as.numeric(b[4, ]))) / 4
}
idx <- grp %in% 1 | grp %in% 0; g <- as.integer(grp[idx] %in% 1)
s <- sc(idx)
r <- roc(g, s, quiet = TRUE, direction = "<")
cat("\n== 联合评分 AUC ==\n"); print(ci.auc(r))
se <- (ci.auc(r)[3] - ci.auc(r)[1]) / (2 * qnorm(0.975))
cat("AUC =", round(as.numeric(auc(r)), 4), " DeLong SE =", round(se, 4), "\n")
for (th in c(0.70, 0.65, 0.60)) {
  z <- (as.numeric(auc(r)) - th) / se
  cat(sprintf("H0: AUC >= %.2f -> z = %.2f, one-sided p = %.4f\n", th, z, pnorm(z)))
}
