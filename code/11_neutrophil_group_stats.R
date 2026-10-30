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
rmat <- read.csv("output/cell_estimates_GSE67530.csv", row.names=1, check.names=FALSE)
pd <- readRDS("data/processed/gse67530_pdata.rds")
pf <- function(x) trimws(sub("^[^:]*:\\s*","",as.character(x)))
grp <- suppressWarnings(as.numeric(pf(pd$characteristics_ch1_4)))
stopifnot(identical(rownames(rmat), colnames(readRDS("data/processed/gse67530_beta.rds"))))
n <- rmat$Neutro
cat("n:", length(n), " groups:", table(grp, useNA="ifany"), "\n\n")
for (g in c("healthy(NA)","ICU(0)","ARDS(1)")) {
  idx <- switch(g, "healthy(NA)"=is.na(grp), "ICU(0)"=grp %in% 0, "ARDS(1)"=grp %in% 1)
  x <- n[idx]
  cat(sprintf("%-12s n=%3d mean=%.4f (%.2f%%) sd=%.4f median=%.4f IQR=%.4f-%.4f\n",
      g, sum(idx), mean(x), 100*mean(x), sd(x), median(x), quantile(x,.25), quantile(x,.75)))
}
a <- n[grp %in% 1]; b <- n[grp %in% 0]
cat("\n-- ARDS vs ICU --\n")
tt <- t.test(a,b); print(tt)
d <- mean(a)-mean(b)
cat(sprintf("diff mean = %+.5f (%.2f pct-point), 95%%CI %.5f to %.5f\n", d, 100*d, tt$conf.int[1], tt$conf.int[2]))
cat(sprintf("Cohen d = %+.4f ; pooled sd = %.4f\n", d/sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2)),
            sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))))
# Wilcoxon + Hodges-Lehmann
w <- wilcox.test(a,b, conf.int=TRUE)
cat("Wilcoxon p =", signif(w$p.value,3), " HL estimate =", signif(w$estimate,5),
    " 95%CI", signif(w$conf.int[1],5), signif(w$conf.int[2],5), "\n")
# Cliff's delta
cd <- function(x,y){ m<-outer(x,y,"-"); (sum(m>0)-sum(m<0))/(length(x)*length(y)) }
cat("Cliff delta =", round(cd(a,b),4), "\n")
