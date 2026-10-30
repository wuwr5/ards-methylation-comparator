# ==============================================================
# 01_download.R
# Download the ARDS methylation datasets and build the beta matrix
# + phenotype tables. Outputs to data/raw/ and data/processed/.
#
# Cohorts:
#   - GSE67530  (Zhang 2019, PMID 31623626): 39 ARDS + 30 healthy controls,
#     whole blood 450K (GPL13534) -> EWAS replication of the 4 signature CpGs.
#     NOTE: GSE32707 is an EXPRESSION array (Dolinay 2012, PMID 22461369,
#     GPL10558), NOT methylation -- never use it for the EWAS. The "99 ARDS +
#     136 healthy" in that paper is the SUM over 5 GEO series (39+18+29+13
#     ARDS; 30+106 healthy), only ONE of which is methylation (GSE67530).
#   - GSE199591 (Calzari 2023): COVID-19 severe/death, blood EPIC -> signature.
#     Series matrix = phenotype only; beta matrix is a ~1 GB suppl file.
#
# DEPENDENCIES (default path): NONE beyond base R.
#   The script downloads the GEO *series matrix* (.txt.gz) directly from
#   NCBI and parses it with base R -- GEOquery is NOT required.
#
#   Why not GEOquery: on this machine the GEOquery install is broken
#   (its dependency httr2 needs rlang >= 1.3.0, only a 1.2.0 binary is
#   available; library(GEOquery) segfaults in R-4.4.2 -> "could not find
#   function getGSE"). A guarded GEOquery fallback is kept for machines
#   where it works.
#
#   Download robustness chain (in order; each validated by gzip magic bytes):
#     1. utils::download.file(method="libcurl", extra="--ssl-no-revoke")
#        -- R >= 4.2 on Windows uses the native schannel TLS stack; the flag
#           disables the CRL revocation check that causes
#           CRYPT_E_REVOCATION_OFFLINE (0x80092013) "SSL connect error".
#     2. libcurl + --ssl-revoke-best-effort / plain libcurl / wininet / auto
#     3. external curl.exe  --ssl-revoke-best-effort -> --ssl-no-revoke
#     4. powershell Invoke-WebRequest (CheckCertificateRevocation = $false)
#     5. manual: browser (or the curl.exe one-liner printed) -> data/raw/
#   Behind a local proxy (Clash/V2Ray etc.)? Set CONFIG$proxy below.
#
#   OPTIONAL : here (falls back to file.path() if absent)
#   OPTIONAL : minfi / ChAMP -- only when CONFIG$process_idat = TRUE
#              (processing RAW IDAT files with champ.load()); the
#              default path reads GEO's processed matrix and skips them.
#
# NOTE: verify the GSE accessions against GEO/PubMed before running.
# ==============================================================

options(timeout = 900)   # GEO files are large; allow long transfers

# ---- project root: ALWAYS resolve to an absolute path ----
# (getwd() depends on the RStudio session; a wrong wd silently scatters
#  outputs. We look for the directory that actually contains R/01_download.R.)
find_project_root <- function() {
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
PROJECT_ROOT <- find_project_root()
here <- function(...) file.path(PROJECT_ROOT, ...)
cat(sprintf("   PROJECT_ROOT (absolute): %s\n", PROJECT_ROOT))
cat("   (all data/results are read from / written under this directory)\n")

CONFIG <- list(
  gse_ewas      = "GSE67530",    # ARDS EWAS (Zhang 2019 PMID 31623626):
                                 #   450K (GPL13534), 39 ARDS + 30 healthy
  gse_ewas_wrong= "GSE32707",    # DO NOT use for EWAS: expression array (GPL10558),
                                 #   Dolinay 2012 PMID 22461369 -- mRNA, not methylation
  gse_death     = "GSE199591",   # COVID-19 severe/death signature (EPIC, PUBLIC)
  # GSE199591's series matrix carries phenotypes ONLY; the processed beta
  # matrix is a ~1.0 GB supplementary file. Set TRUE to fetch it (slow!):
  gse_death_suppl       = "GSE199591_COVID_Matrix_processed.txt.gz",
  download_death_suppl  = FALSE,
  bradic_pmid  = "36371196",     # controlled-access; request from author if used
  dest_raw     = here("data", "raw"),
  dest_proc    = here("data", "processed"),
  proxy        = NULL,           # e.g. "http://127.0.0.1:7890" if you run a local proxy
  # where to look for a MANUALLY downloaded series matrix (browser default)
  extra_search_dirs = c(here("data", "processed"),
                        path.expand("~/Downloads"), path.expand("~/Desktop")),
  process_idat = FALSE,          # TRUE + idat_dir -> champ.load() (needs minfi/ChAMP)
  idat_dir     = NULL
)

dir.create(CONFIG$dest_raw,  showWarnings = FALSE, recursive = TRUE)
dir.create(CONFIG$dest_proc, showWarnings = FALSE, recursive = TRUE)

# ==============================================================
# helpers
# ==============================================================

# GEO sub-directory naming: GSE32707 -> GSE32nnn (keep all but last 3 digits)
gse_subdir <- function(acc) sub("([0-9])[0-9]{3}$", "\\1nnn", acc)

gse_series_url <- function(acc)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s_series_matrix.txt.gz",
          gse_subdir(acc), acc, acc)

# Alternative endpoint (www host) -- useful when the ftp host is blocked.
gse_series_url_alt <- function(acc)
  sprintf("https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s&format=file&file=%s_series_matrix.txt.gz",
          acc, acc)

gse_suppl_url <- function(acc)
  sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/",
          gse_subdir(acc), acc)

clean_tok <- function(x) trimws(gsub('^"|"$', "", x))

is_gzip <- function(path) {
  if (!file.exists(path)) return(FALSE)
  con <- file(path, "rb"); on.exit(close(con))
  magic <- readBin(con, "raw", n = 2)
  length(magic) == 2 && magic[1] == as.raw(0x1f) && magic[2] == as.raw(0x8b)
}

# Find an already-downloaded series matrix: project data/raw first, then
# data/processed (some users drop browser downloads there), ~/Downloads, Desktop.
# Returns an absolute path or NULL.
locate_gz <- function(acc, extra_dirs = character()) {
  fname <- sprintf("%s_series_matrix.txt.gz", acc)
  cands <- c(file.path(PROJECT_ROOT, "data", "raw", fname),
             file.path(extra_dirs, fname))
  for (p in cands) {
    if (file.exists(p) && is_gzip(p)) return(normalizePath(p))
  }
  NULL
}

set_proxy_env <- function(proxy) {
  if (!is.null(proxy) && nzchar(proxy)) {
    Sys.setenv("http_proxy"  = proxy, "HTTP_PROXY"  = proxy,
               "https_proxy" = proxy, "HTTPS_PROXY" = proxy)
  }
  invisible(NULL)
}

# External curl.exe (own TLS stack, honours http_proxy/https_proxy env).
# `ssl_flags` are tried in order: --ssl-revoke-best-effort first (only skips
# the CRL check when the revocation server is offline -- the exact
# CRYPT_E_REVOCATION_OFFLINE case), then --ssl-no-revoke (always skip),
# then no flag. Returns TRUE only for a valid non-empty gzip file.
download_via_curl_exe <- function(url, dest, proxy = NULL,
                                  ssl_flags = c("--ssl-revoke-best-effort",
                                                "--ssl-no-revoke", "")) {
  curl <- Sys.which("curl")
  if (!nzchar(curl)) {
    for (p in c("C:/Windows/System32/curl.exe", "C:/Windows/SysWOW64/curl.exe")) {
      if (file.exists(p)) { curl <- p; break }
    }
  }
  if (!nzchar(curl)) return(FALSE)
  # NCBI throttles non-browser UAs ~5x on some egress IPs; send a browser UA.
  ua <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/130.0.0.0 Safari/537.36"
  set_proxy_env(proxy)
  for (fl in unique(ssl_flags)) {
    args <- c("-sSL", "--fail", "--connect-timeout", "30", "--max-time", "1800",
              "-A", ua, "--retry", "5", "--retry-delay", "10")
    if (nzchar(fl)) args <- c(args, fl)
    if (!is.null(proxy) && nzchar(proxy)) args <- c(args, "-x", proxy)
    args <- c(args, "-o", paste0('"', dest, '"'), paste0('"', url, '"'))
    ok <- tryCatch({
      rc <- system2(curl, args)
      (rc == 0) && file.exists(dest) && is_gzip(dest)
    }, error = function(e) FALSE)
    if (isTRUE(ok)) return(TRUE)
    if (file.exists(dest)) unlink(dest)
  }
  FALSE
}

# PowerShell Invoke-WebRequest fallback (.NET TLS stack, revocation disabled).
# PS 5.1: uses ServicePointManager.CheckCertificateRevocation = $false;
# PS 7+  : -SkipCertificateCheck is supported directly.
download_via_powershell <- function(url, dest, proxy = NULL) {
  ps <- Sys.which("powershell")
  if (!nzchar(ps)) return(FALSE)
  set_proxy_env(proxy)
  script <- file.path(tempdir(), "ards_dl.ps1")
  dest_fwd <- gsub("\\\\", "/", dest)
  writeLines(c(
    "$ProgressPreference = 'SilentlyContinue'",
    "try { [Net.ServicePointManager]::CheckCertificateRevocation = $false } catch {}",
    sprintf("$u = '%s'", gsub("'", "''", url)),
    sprintf("$o = '%s'", gsub("'", "''", dest_fwd)),
    "try {",
    "    Invoke-WebRequest -Uri $u -OutFile $o -UseBasicParsing -SkipCertificateCheck -ErrorAction Stop",
    "} catch {",
    "    Invoke-WebRequest -Uri $u -OutFile $o -UseBasicParsing -ErrorAction Stop",
    "}"
  ), script)
  ok <- tryCatch({
    rc <- system2(ps, c("-NoProfile", "-ExecutionPolicy", "Bypass",
                        "-File", shQuote(script)))
    (rc == 0) && file.exists(dest) && is_gzip(dest)
  }, error = function(e) FALSE)
  if (!isTRUE(ok) && file.exists(dest)) unlink(dest)
  isTRUE(ok)
}

# Try several R download methods and an external curl, over two NCBI
# endpoints. Every attempt is validated with a gzip magic-number check,
# so HTML error pages / truncated transfers are rejected, not silently used.
download_robust <- function(url, gz, proxy = NULL, verbose = TRUE) {
  acc  <- sub("^.*/matrix/(GSE[0-9]+)_series_matrix\\.txt\\.gz$", "\\1", url)
  urls <- unique(c(url, gse_series_url_alt(acc)))
  # (method, extra) pairs. R >= 4.2 on Windows uses schannel; --ssl-no-revoke
  # disables the CRL check behind CRYPT_E_REVOCATION_OFFLINE. Unknown extra
  # options are ignored by R's libcurl module, so these are always safe.
  attempts <- list(
    list(m = "libcurl", extra = "--ssl-no-revoke"),
    list(m = "libcurl", extra = "--ssl-revoke-best-effort"),
    list(m = "libcurl", extra = NULL),
    list(m = "wininet",  extra = NULL),
    list(m = "auto",     extra = NULL)
  )
  for (u in urls) {
    for (a in attempts) {
      extra <- a$extra
      if (!is.null(proxy) && nzchar(proxy) && a$m == "libcurl")
        extra <- unique(c(extra, "-x", proxy))
      ok <- tryCatch({
        utils::download.file(u, gz, mode = "wb", method = a$m, quiet = TRUE,
                             cacheOK = FALSE, extra = extra)
        file.exists(gz) && is_gzip(gz)
      }, error = function(e) {
        if (verbose) cat(sprintf("   !! method='%s' error: %s\n", a$m, conditionMessage(e)))
        FALSE
      }, warning = function(w) {
        if (verbose) cat(sprintf("   !! method='%s' warning: %s\n", a$m, conditionMessage(w)))
        FALSE
      })
      if (isTRUE(ok)) {
        if (verbose) cat(sprintf("   downloaded (method='%s', extra='%s', host=%s)\n",
                                 a$m, if (is.null(a$extra)) "-" else a$extra,
                                 sub("^https://([^/]+)/.*$", "\\1", u)))
        return(TRUE)
      }
      if (file.exists(gz)) unlink(gz)
    }
  }
  if (verbose) cat("   -> trying external curl.exe (schannel --ssl-no-revoke) ...\n")
  if (download_via_curl_exe(url, gz, proxy = proxy)) {
    if (verbose) cat("   downloaded (external curl.exe)\n")
    return(TRUE)
  }
  if (verbose) cat("   -> trying powershell Invoke-WebRequest ...\n")
  if (download_via_powershell(url, gz, proxy = proxy)) {
    if (verbose) cat("   downloaded (powershell Invoke-WebRequest)\n")
    return(TRUE)
  }
  FALSE
}

# Phenotype table from the !Sample_* header lines (works even when the
# series matrix has NO probe-level table, e.g. GSE199591).
# NOTE: GEO repeats the key "!Sample_characteristics_ch1" for every
# characteristic line -- we must NOT overwrite; suffix them _2, _3, ...
extract_pheno <- function(lines) {
  pheno <- list()
  for (kl in grep("^!Sample_", lines, value = TRUE)) {
    parts <- strsplit(kl, "\t", fixed = TRUE)[[1]]
    key <- sub("^!Sample_", "", parts[1L])
    vals <- clean_tok(parts[-1L])
    base <- key; k <- 1L
    while (key %in% names(pheno)) { k <- k + 1L; key <- paste0(base, "_", k) }
    pheno[[key]] <- vals
  }
  if (!length(pheno)) return(NULL)
  n_smp <- length(pheno[[1L]])
  pd <- as.data.frame(lapply(pheno, function(v) { length(v) <- n_smp; v }),
                      stringsAsFactors = FALSE)
  pd
}

# Parse a GEO series matrix (already read into `lines`) into
#   beta : numeric matrix, probes x samples (rownames = probe IDs)
#   pd   : data.frame, one row per sample, from the !Sample_* header lines
parse_series_lines <- function(lines) {
  b <- grep("^!series_matrix_table_begin", lines)
  e <- grep("^!series_matrix_table_end",  lines)
  if (length(b) != 1L || length(e) != 1L || e <= b + 1L)
    stop("no probe-level table in series matrix (IDAT-only series?)")

  pd <- extract_pheno(lines)

  # ---- probe table ----
  hdr <- clean_tok(strsplit(lines[b + 1L], "\t", fixed = TRUE)[[1]])
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  writeLines(lines[(b + 2L):(e - 1L)], tmp)
  m <- utils::read.delim(tmp, header = FALSE, sep = "\t", quote = "",
                         col.names = hdr, check.names = FALSE,
                         colClasses = c("character", rep("numeric", length(hdr) - 1L)),
                         na.strings = c("NA", "null", "NaN", ""))
  # GEO writes the ID_REF values WITH literal double quotes (e.g. "cg00000029");
  # since we read with quote="", those quotes are preserved as part of the
  # string. Strip them so probe IDs match downstream CpG lookups exactly.
  ids <- clean_tok(as.character(m[[1L]]))
  rownames(m) <- make.unique(ids)
  m[[1L]] <- NULL
  list(beta = as.matrix(m), pd = pd)
}

# GEO supplementary processed beta matrix (RnBeads export, e.g.
# GSE199591_COVID_Matrix_processed.txt.gz): header row = sample IDs,
# first column = CpG ID, remaining columns = beta values.
parse_suppl_matrix <- function(gz_path, chunk = 50000L) {
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::fread(gz_path, sep = "\t", header = TRUE, showProgress = FALSE)
    m <- as.matrix(dt[, -1, with = FALSE])
    rownames(m) <- dt[[1]]
    return(m)
  }
  con <- gzfile(gz_path, "rt"); on.exit(close(con), add = TRUE)
  hdr <- clean_tok(strsplit(readLines(con, n = 1L, warn = FALSE), "\t",
                            fixed = TRUE)[[1]])
  samples <- hdr[-1L]
  probes <- character(); blocks <- list()
  repeat {
    lines <- readLines(con, n = chunk, warn = FALSE)
    if (!length(lines)) break
    parts <- strsplit(lines, "\t", fixed = TRUE)
    probes <- c(probes, vapply(parts, `[`, "", 1L))
    blocks[[length(blocks) + 1L]] <-
      do.call(rbind, lapply(parts, function(p) as.numeric(p[-1L])))
  }
  m <- do.call(rbind, blocks)
  dimnames(m) <- list(probes, samples)
  m
}

# Memory-efficient series-matrix parser for VERY large files.
# GSE67530's series matrix is ~568 MB gz (~1.33 GB uncompressed,
# 485,577 probes x 144 full-precision beta values): reading the whole
# text into R would need >3 GB RAM, so stream the probe table in chunks.
parse_series_gz <- function(gz_path, chunk_rows = 50000L) {
  con <- gzfile(gz_path, "rt")
  on.exit(if (isOpen(con)) close(con), add = TRUE)

  # 1) header lines up to and including !series_matrix_table_begin
  hdr_lines <- character()
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (!length(l)) stop("!series_matrix_table_begin not found in ", gz_path)
    hdr_lines <- c(hdr_lines, l)
    if (grepl("^!series_matrix_table_begin", l)) break
  }
  pd <- extract_pheno(hdr_lines)

  # 2) column header row ("ID_REF", GSM..., GSM...)
  col_line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(col_line)) stop("truncated series matrix (no column header)")
  hdr <- clean_tok(strsplit(col_line, "\t", fixed = TRUE)[[1]])
  nsmp <- length(hdr) - 1L

  # 3) stream data rows into a temp file, dropping the end marker
  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  tc <- file(tmp, "wt")
  done <- FALSE
  repeat {
    l <- readLines(con, n = chunk_rows, warn = FALSE)
    if (!length(l)) break
    keep <- l[!grepl("^!series_matrix_table_end", l)]
    if (length(keep) < length(l)) done <- TRUE
    if (length(keep)) writeLines(keep, tc)
    if (done) break
  }
  close(tc)

  m <- utils::read.delim(tmp, header = FALSE, sep = "\t", quote = "",
                         col.names = hdr, check.names = FALSE,
                         colClasses = c("character", rep("numeric", nsmp)),
                         na.strings = c("NA", "null", "NaN", ""))
  # GEO writes ID_REF values with literal double quotes; strip them so probe
  # IDs (cgXXXXXXX) match downstream CpG lookups exactly.
  ids <- clean_tok(as.character(m[[1L]]))
  rownames(m) <- make.unique(ids)
  m[[1L]] <- NULL
  list(beta = as.matrix(m), pd = pd)
}

# Print the candidate phenotype columns so the user can code `group`.
print_group_help <- function(acc, pd) {
  cat(sprintf("   [%s] candidate grouping columns (%d samples):\n", acc, nrow(pd)))
  ch <- grep("characteristics|title|description", names(pd),
             ignore.case = TRUE, value = TRUE)
  if (!length(ch)) ch <- names(pd)[seq_len(min(3L, ncol(pd)))]
  for (cc in ch[seq_len(min(4L, length(ch)))]) {
    tb <- table(pd[[cc]])
    cat(sprintf("     - %s : %d unique -> %s\n", cc, length(tb),
                paste(names(tb), tb, sep = ":", collapse = " | ")))
  }
  cat("     >>> MANUALLY code pd$group, e.g.\n")
  cat(sprintf("         pd <- readRDS(file.path('data','processed','%s_pdata.rds'))\n", acc))
  cat("         pd$group <- ifelse(grepl('ARDS', pd[[1]]), 'ARDS', 'control')\n")
  cat(sprintf("         saveRDS(pd, file.path('data','processed','%s_pdata.rds'))\n", acc))
}

print_manual_route <- function(acc, gz, proxy) {
  cat("   !! all automatic download methods failed. Manual route:\n")
  cat("      1a. FASTEST -- paste this ONE line into cmd.exe (Windows 10+):\n")
  cat(sprintf("         curl.exe --ssl-no-revoke -L -o \"%s\" \"%s\"\n",
              gsub("/", "\\\\", gz), gse_series_url(acc)))
  cat("         (--ssl-no-revoke skips the CRL check that fails with\n")
  cat("          CRYPT_E_REVOCATION_OFFLINE / 'SSL connect error')\n")
  cat("      1b. OR download in a browser (Edge/Chrome usually works):\n")
  cat("         ", gse_series_url(acc), "\n", sep = "")
  cat("         (alt: ", gse_series_url_alt(acc), ")\n", sep = "")
  cat("      2. Save the file EXACTLY as:\n         ", gz, "\n", sep = "")
  cat("      3. Re-run this script -- the cached file will be detected.\n")
  if (is.null(proxy) || !nzchar(proxy)) {
    cat("      (Using a local proxy (Clash/V2Ray/...)? Set CONFIG$proxy <-\n")
    cat("       'http://127.0.0.1:7890' and re-run.)\n")
  }
}

# Optional fallback for machines with a working GEOquery install.
geoquery_fallback <- function(acc, prefix, dest_raw, dest_proc) {
  cat("   -> trying GEOquery fallback ...\n")
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    cat("   GEOquery not installed. Fix network first, or install it:\n")
    cat("       BiocManager::install('GEOquery')\n")
    return(invisible(NULL))
  }
  if (!"getGSE" %in% getNamespaceExports("GEOquery")) {
    cat("   GEOquery installed but BROKEN (no getGSE export). Reinstall:\n")
    cat("       remove.packages('GEOquery'); BiocManager::install('GEOquery')\n")
    return(invisible(NULL))
  }
  g <- tryCatch(GEOquery::getGSE(acc, destdir = dest_raw, getGPL = FALSE),
                error = function(e) {
                  cat("   getGSE error:", conditionMessage(e), "\n"); NULL
                })
  if (is.null(g) || !length(g)) return(invisible(NULL))
  beta <- Biobase::exprs(g[[1L]])
  pd   <- as.data.frame(Biobase::pData(Biobase::phenoData(g[[1L]])))
  pd$geo_accession <- rownames(pd)
  cat(sprintf("   beta matrix: %d probes x %d samples\n", nrow(beta), ncol(beta)))
  saveRDS(beta, file.path(dest_proc, sprintf("%s_beta.rds", prefix)))
  saveRDS(pd,   file.path(dest_proc, sprintf("%s_pdata.rds", prefix)))
  cat(sprintf("   saved %s_beta.rds / %s_pdata.rds (via GEOquery)\n", prefix, prefix))
  invisible(list(beta = beta, pd = pd))
}

# GSE199591: phenotype is in the series matrix, but the beta matrix is a
# ~1 GB supplementary file. Fetch it only when explicitly enabled.
fetch_death_suppl <- function(acc, prefix, dest_raw, dest_proc, proxy = NULL) {
  fname <- CONFIG$gse_death_suppl
  if (!nzchar(fname)) return(invisible(NULL))
  url <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/suppl/%s",
                 gse_subdir(acc), acc, fname)
  gz  <- file.path(dest_raw, fname)
  cat(sprintf("   beta matrix lives in the supplementary file (~1 GB):\n      %s\n", url))
  if (!isTRUE(CONFIG$download_death_suppl)) {
    cat("      To fetch it automatically, set CONFIG$download_death_suppl <- TRUE\n")
    cat("      (slow!), or download it yourself with ONE cmd.exe line:\n")
    cat(sprintf("        curl.exe --ssl-no-revoke -L -o \"%s\" \"%s\"\n",
                gsub("/", "\\\\", gz), url))
    cat("      then re-run this script (the cached file is detected).\n")
    return(invisible(NULL))
  }
  if (!(file.exists(gz) && is_gzip(gz))) {
    ok <- download_robust(url, gz, proxy = proxy)
    if (!isTRUE(ok)) { cat("   !! supplementary download failed.\n"); return(invisible(NULL)) }
  }
  cat("   parsing processed beta matrix (866k probes -- may take minutes) ...\n")
  m <- parse_suppl_matrix(gz)
  cat(sprintf("   beta matrix: %d probes x %d samples\n", nrow(m), ncol(m)))
  saveRDS(m, file.path(dest_proc, sprintf("%s_beta.rds", prefix)))
  cat(sprintf("   saved %s_beta.rds -> data/processed/\n", prefix))
  invisible(m)
}

# Download + parse one cohort; returns list(beta, pd) or NULL on failure.
fetch_cohort <- function(acc, prefix, dest_raw, dest_proc, proxy = NULL) {
  cat(sprintf("\n>> [%s] series matrix ...\n", acc))
  gz  <- file.path(dest_raw, sprintf("%s_series_matrix.txt.gz", acc))
  url <- gse_series_url(acc)

  # 1) an existing valid file anywhere known (data/raw, data/processed,
  #    ~/Downloads, Desktop)
  cached <- locate_gz(acc, CONFIG$extra_search_dirs)
  if (!is.null(cached)) {
    cat("   using existing file:", cached, "\n")
    if (normalizePath(cached) != normalizePath(gz)) {
      dir.create(dest_raw, showWarnings = FALSE, recursive = TRUE)
      file.copy(cached, gz, overwrite = TRUE)
      cat("   copied into:", gz, "\n")
    }
  } else {
    if (file.exists(gz)) unlink(gz)   # stale non-gzip residue
    ok <- download_robust(url, gz, proxy = proxy)
    if (!isTRUE(ok)) {
      cat("   !! direct download failed.\n")
      print_manual_route(acc, gz, proxy)
      return(geoquery_fallback(acc, prefix, dest_raw, dest_proc))
    }
  }

  # 2) parse -- stream large files (GSE67530 is ~568 MB gz) instead of
  #    reading the whole text into RAM
  big <- file.info(gz)$size > 50e6
  cat(sprintf("   parsing series matrix (%.1f MB, %s) ...\n",
              file.info(gz)$size / 1e6, if (big) "streaming" else "in-memory"))
  res <- tryCatch({
    if (big) {
      parse_series_gz(gz)
    } else {
      con <- gzfile(gz, "rt"); ln <- readLines(con, warn = FALSE); close(con)
      parse_series_lines(ln)
    }
  }, error = function(e) {
    cat("   !! parse error:", conditionMessage(e), "\n"); NULL
  })

  if (is.null(res) || nrow(res$beta) == 0L) {
    # No probe table in the series matrix -- phenotype may still be there.
    pd <- if (!is.null(res)) res$pd else {
      con <- gzfile(gz, "rt"); ln <- readLines(con, warn = FALSE); close(con)
      tryCatch(extract_pheno(ln), error = function(e) NULL)
    }
    if (!is.null(pd) && nrow(pd) > 0L) {
      saveRDS(pd, file.path(dest_proc, sprintf("%s_pdata.rds", prefix)))
      cat(sprintf("   series matrix has no probe table; saved %s_pdata.rds ",
                  prefix))
      cat(sprintf("(%d samples x %d cols)\n", nrow(pd), ncol(pd)))
      print_group_help(acc, pd)
    }
    if (identical(acc, CONFIG$gse_death))
      return(fetch_death_suppl(acc, prefix, dest_raw, dest_proc, proxy = proxy))
    cat("      Check supplementary files at:\n        ", gse_suppl_url(acc), "\n", sep = "")
    return(invisible(NULL))
  }

  cat(sprintf("   beta matrix: %d probes x %d samples\n",
              nrow(res$beta), ncol(res$beta)))
  print_group_help(acc, res$pd)
  saveRDS(res$beta, file.path(dest_proc, sprintf("%s_beta.rds", prefix)))
  saveRDS(res$pd,   file.path(dest_proc, sprintf("%s_pdata.rds", prefix)))
  cat(sprintf("   saved %s_beta.rds / %s_pdata.rds -> data/processed/\n", prefix, prefix))
  invisible(gc(verbose = FALSE))
  invisible(res)
}

# ==============================================================
# MAIN -- set options(ards.run_download = FALSE) to skip (unit testing)
# ==============================================================
if (isTRUE(getOption("ards.run_download", TRUE))) {

  fetch_cohort(CONFIG$gse_ewas,  "gse67530",
               CONFIG$dest_raw, CONFIG$dest_proc, proxy = CONFIG$proxy)
  fetch_cohort(CONFIG$gse_death, "gse199591",
               CONFIG$dest_raw, CONFIG$dest_proc, proxy = CONFIG$proxy)

  # --------------------------------------------------------------
  # (Optional) Process RAW IDATs with ChAMP/minfi -- needs those packages.
  #   Only runs when CONFIG$process_idat = TRUE. Skipped by default.
  # --------------------------------------------------------------
  if (isTRUE(CONFIG$process_idat)) {
    if (!requireNamespace("ChAMP", quietly = TRUE) ||
        !requireNamespace("minfi", quietly = TRUE)) {
      stop("process_idat = TRUE but ChAMP/minfi not installed. Install via ",
           "BiocManager::install(c('minfi','ChAMP')), or keep process_idat = ",
           "FALSE to use the GEO processed matrix.")
    }
    suppressMessages({ library(ChAMP); library(minfi) })
    cat(">> Processing raw IDATs with champ.load() ...\n")
    # myLoad <- champ.load(dir = CONFIG$idat_dir, arraytype = "450K")
    # beta_mat <- myLoad$beta
    # pd <- as.data.frame(pData(myLoad$pd)); pd$geo_accession <- rownames(pd)
    # ... then overwrite the saved rds files above.
    cat("   >>> Uncomment the champ.load block above and supply CONFIG$idat_dir.\n")
  }

  # --------------------------------------------------------------
  # (Optional) Bradic 2022 (PMID 36371196) -- CONTROLLED ACCESS.
  #   Data NOT in GEO/public repos; request from corresponding author
  #   (Charbel Abi Khalil) with a signed data access agreement.
  #   If obtained, load below and use as an independent validation set:
  #   bradic    <- readRDS(<your_bradic_beta.rds>)
  #   bradic_pd <- readRDS(<your_bradic_pdata.rds>)  # needs age, sex, time, event
  # --------------------------------------------------------------

  cat("\n>> Done. Next steps:\n")
  cat("   1. Code pd$group for each cohort (see the printed tables above).\n")
  cat("   2. source('R/02_ewas.R')                 # GSE67530 4-CpG AUC\n")
  cat("   3. source('R/06_death_signature_gse199591.R')  # severe/death signature\n")
  cat("      (needs the ~1 GB GSE199591 suppl beta matrix -- see notes above)\n")

  # ---- final verification with ABSOLUTE paths (no guessing) ----
  cat("\n>> Verification:\n")
  expected <- c(file.path(CONFIG$dest_proc, "gse67530_beta.rds"),
                file.path(CONFIG$dest_proc, "gse67530_pdata.rds"),
                file.path(CONFIG$dest_proc, "gse199591_pdata.rds"),
                file.path(CONFIG$dest_proc, "gse199591_beta.rds"))
  n_ok <- 0L
  for (f in expected) {
    ok <- file.exists(f)
    n_ok <- n_ok + ok
    note <- ""
    if (!ok && grepl("gse199591_beta", f))
      note <- "   (optional for now: needs the ~1 GB suppl file)"
    cat(sprintf("   %s  %s%s\n", if (ok) "OK     " else "MISSING", f, note))
  }
  if (!file.exists(file.path(CONFIG$dest_proc, "gse67530_beta.rds"))) {
    cat("   !! EWAS inputs missing -> the GSE67530 download did NOT succeed.\n")
    cat("      Use the manual route printed above (browser or the curl.exe\n")
    cat("      one-liner), then re-run.\n")
  }
}
