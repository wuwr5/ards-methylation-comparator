## 48_table3.R — 从 β 矩阵重算表3:组均值、Δβ、95%CI、AUC,确保与 logFC 自洽
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
dir_sign <- c(1, 1, -1, -1)

sc <- function(idx) {
  b <- beta[probes, idx, drop = FALSE]
  as.numeric(colMeans(sweep(b, 1, dir_sign, "*")))
}
ctr <- list(ARDS_vs_healthy = list(cs = grp %in% 1, ct = is.na(grp)),
            ICU_vs_healthy  = list(cs = grp %in% 0, ct = is.na(grp)),
            ARDS_vs_ICU     = list(cs = grp %in% 1, ct = grp %in% 0))

out <- list()
for (nm in names(ctr)) {
  cs <- ctr[[nm]]$cs; ct <- ctr[[nm]]$ct; idx <- cs | ct
  g <- as.integer(cs[idx]); lab <- ifelse(g == 1, 1, 0)
  for (i in seq_along(probes)) {
    y  <- as.numeric(beta[probes[i], idx])
    tt <- t.test(y[g == 1], y[g == 0])
    r  <- roc(g, y, quiet = TRUE, direction = "<")
    ci <- as.numeric(ci.auc(r))
    out[[length(out) + 1]] <- data.frame(contrast = nm, probe = probes[i],
      mean_case = mean(y[g == 1]), mean_ctrl = mean(y[g == 0]),
      delta = mean(y[g == 1]) - mean(y[g == 0]),
      lo95 = tt$conf.int[1], hi95 = tt$conf.int[2],
      auc = as.numeric(auc(r)), auc_lo = ci[1], auc_hi = ci[3], row.names = NULL)
  }
  s <- sc(idx); r <- roc(g, s, quiet = TRUE, direction = "<")
  ci <- as.numeric(ci.auc(r))
  out[[length(out) + 1]] <- data.frame(contrast = nm, probe = "combined",
    mean_case = NA, mean_ctrl = NA, delta = NA, lo95 = NA, hi95 = NA,
    auc = as.numeric(auc(r)), auc_lo = ci[1], auc_hi = ci[3])
}
res <- do.call(rbind, out)
cat("== 重算结果 ==\n"); print(res, digits = 4, row.names = FALSE)
write.csv(res, "output/table3_recomputed.csv", row.names = FALSE)

cat("\n== 关键量 ==\n")
p <- res[res$contrast == "ARDS_vs_healthy" & res$probe != "combined", ]
q <- res[res$contrast == "ICU_vs_healthy"  & res$probe != "combined", ]
r <- res[res$contrast == "ARDS_vs_ICU"     & res$probe != "combined", ]
cat(sprintf("ARDS vs healthy  Δβ 范围: %.3f ~ %.3f\n", min(p$delta), max(p$delta)))
cat(sprintf("ICU  vs healthy  Δβ 范围: %.3f ~ %.3f\n", min(q$delta), max(q$delta)))
cat(sprintf("ARDS vs ICU  Δβ: %s\n", paste(sprintf("%.3f(%.3f,%.3f)", r$delta, r$lo95, r$hi95), collapse = " ")))
cat(sprintf("|Δβ| 最小病例-健康效应: %.3f ; 其一半: %.3f\n", min(abs(p$delta)), min(abs(p$delta)) / 2))
cat(sprintf("ARDS vs ICU |Δβ| 上界最大: %.3f\n", max(abs(c(r$lo95, r$hi95)))))
cat(sprintf("比值 上界: %s\n", paste(sprintf("%.2f", abs(c(r$lo95, r$hi95))[c(2,4,6,8)] / abs(p$delta)), collapse = " ")))
