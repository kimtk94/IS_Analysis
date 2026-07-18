#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) { i <- which(args == flag); if (!length(i) || i == length(args)) return(default); args[i + 1] }
has_flag <- function(flag) flag %in% args

# ---------- CLI ----------
gene_file <- get_arg("--gene-file")
source_file_list <- get_arg("--source-file-list")
batch_id  <- get_arg("--batch-id", "batch_001")
outdir    <- get_arg("--outdir", "results/exposure_batches")
batch_output <- get_arg("--batch-output")
rawdir    <- get_arg("--rawdir", "data/rawdata/pqtl/selected_targets")
tmp_root  <- get_arg("--tmpdir", "/content/ukbppp_tmp")
ancestries <- trimws(unlist(strsplit(get_arg("--ancestries", "EUR,EAS"), ",")))
p_threshold <- as.numeric(get_arg("--p-threshold", "5e-8"))
cis_window  <- as.integer(get_arg("--cis-window", "1000000"))

eur_first <- has_flag("--eur-first")
copy_to_local <- has_flag("--copy-to-local")
keep_all <- has_flag("--keep-all")
no_cis_filter <- has_flag("--no-cis-filter")
no_p_filter <- has_flag("--no-p-filter")
no_f_filter <- has_flag("--no-f-filter")
force <- has_flag("--force")
test_mode <- has_flag("--test")
max_file_lines <- as.integer(get_arg("--max-file-lines", if (test_mode) "200" else "0"))
if (is.na(max_file_lines) || max_file_lines < 0) stop("[ERROR] --max-file-lines must be a non-negative integer")

if (is.null(gene_file)) stop("[ERROR] --gene-file is required")
if (!file.exists(gene_file)) stop("[ERROR] gene file does not exist: ", gene_file)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_root, recursive = TRUE, showWarnings = FALSE)
per_gene_dir <- file.path(outdir, "per_gene", batch_id)
log_dir <- file.path(outdir, "logs")
dir.create(per_gene_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

batch_out <- if (is.null(batch_output)) file.path(outdir, paste0("exposure_", batch_id, ".tsv")) else batch_output
dir.create(dirname(batch_out), recursive = TRUE, showWarnings = FALSE)
gene_status_out <- file.path(log_dir, paste0(batch_id, "_gene_status.tsv"))

msg <- function(...) message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " ", paste0(..., collapse = ""))
safe_upper <- function(x) toupper(trimws(as.character(x)))
safe_label <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

target_genes <- unique(safe_upper(readLines(gene_file, warn = FALSE)))
target_genes <- target_genes[target_genes != ""]
if (!length(target_genes)) stop("[ERROR] gene file has no genes")

allowed_source_files <- NULL
if (!is.null(source_file_list)) {
  if (!file.exists(source_file_list)) stop("[ERROR] source file list does not exist: ", source_file_list)
  allowed_source_files <- unique(trimws(readLines(source_file_list, warn = FALSE)))
  allowed_source_files <- allowed_source_files[allowed_source_files != ""]
  if (!length(allowed_source_files)) stop("[ERROR] source file list is empty: ", source_file_list)
}

msg("[INFO] batch_id=", batch_id, " genes=", length(target_genes), " eur_first=", eur_first, " copy_to_local=", copy_to_local, " test=", test_mode, " max_file_lines=", max_file_lines)
msg("[INFO] rawdir=", rawdir, " outdir=", outdir, " tmpdir=", tmp_root)
if (!is.null(allowed_source_files)) msg("[INFO] exact source_file filter enabled: ", length(allowed_source_files), " files")

# ---------- column mapping ----------
pick_col <- function(nms, candidates) {
  lower_map <- setNames(nms, tolower(nms))
  for (cand in candidates) {
    key <- tolower(cand)
    if (key %in% names(lower_map)) return(lower_map[[key]])
  }
  NA_character_
}

COLS <- list(
  chr = c("Chrom", "chr", "chromosome", "CHR", "#CHROM"),
  pos = c("Pos(hg38)", "pos_hg38", "position_hg38", "Pos", "POS", "position", "bp"),
  snp = c("rsids", "rsid", "SNP", "variant", "Name", "MarkerName", "ID"),
  ea = c("effectAllele", "effect_allele", "EA", "A1", "ALT", "testedAllele"),
  oa = c("otherAllele", "other_allele", "OA", "A2", "REF", "nonEffectAllele"),
  beta = c("Beta(SD)", "Beta", "beta", "BETA", "effect", "Effect"),
  se = c("SE", "se", "standard_error", "StdErr"),
  p = c("Pval", "pval", "p_value", "P", "p", "Pvalue", "P-value"),
  mlogp = c("min_log10_pval", "minus_log10_pval", "neg_log10_pval", "log10p"),
  n = c("N", "n", "samplesize", "sample_size", "OBS_CT"),
  eaf = c("ImpMAF", "EAF", "eaf", "effect_allele_frequency", "effectAlleleFrequency", "MAF", "maf")
)

choose_needed_columns <- function(nms) unique(na.omit(c(
  pick_col(nms, COLS$chr), pick_col(nms, COLS$pos), pick_col(nms, COLS$snp),
  pick_col(nms, COLS$ea), pick_col(nms, COLS$oa), pick_col(nms, COLS$beta),
  pick_col(nms, COLS$se), pick_col(nms, COLS$p), pick_col(nms, COLS$mlogp),
  pick_col(nms, COLS$n), pick_col(nms, COLS$eaf)
)))

empty_exposure_dt <- function() {
  data.table(
    gene_symbol=character(), ancestry=character(), source_file=character(),
    chr=character(), pos_hg38=integer(), effect_allele=character(), other_allele=character(),
    beta=numeric(), se=numeric(), SNP=character(), pval=numeric(), samplesize=numeric(), eaf=numeric(),
    exposure=character(), id.exposure=character(), beta.exposure=numeric(), se.exposure=numeric(),
    effect_allele.exposure=character(), other_allele.exposure=character(),
    eaf.exposure=numeric(), pval.exposure=numeric(), samplesize.exposure=numeric(),
    is_cis=logical(), F_stat=numeric()
  )
}

write_empty_exposure <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(empty_exposure_dt(), path, sep = "\t")
}

parse_variant_id <- function(x) {
  # Supports chr1:12345:A:G, 1:12345:A:G, chr1_12345_A_G, and chr1:12345_A_G.
  # Allele1 is interpreted as other/ref and allele2 as effect/alt.
  x <- as.character(x)
  clean <- gsub("^chr", "", x, ignore.case = TRUE)
  clean <- gsub("[|/]", ":", clean)
  clean <- gsub("_", ":", clean)
  parts <- tstrsplit(clean, ":", fixed = TRUE, fill = NA)

  n <- length(parts)
  chr <- if (n >= 1) parts[[1]] else rep(NA_character_, length(x))
  pos <- if (n >= 2) suppressWarnings(as.integer(parts[[2]])) else rep(NA_integer_, length(x))
  a1 <- if (n >= 3) safe_upper(parts[[3]]) else rep(NA_character_, length(x))
  a2 <- if (n >= 4) safe_upper(parts[[4]]) else rep(NA_character_, length(x))

  data.table(chr_from_id = chr, pos_from_id = pos, other_from_id = a1, effect_from_id = a2)
}

standardize_pqtl <- function(dt, gene, ancestry, source_file) {
  nms <- names(dt)
  c_chr <- pick_col(nms, COLS$chr); c_pos <- pick_col(nms, COLS$pos)
  c_snp <- pick_col(nms, COLS$snp); c_ea <- pick_col(nms, COLS$ea); c_oa <- pick_col(nms, COLS$oa)
  c_beta <- pick_col(nms, COLS$beta); c_se <- pick_col(nms, COLS$se)
  c_p <- pick_col(nms, COLS$p); c_mlogp <- pick_col(nms, COLS$mlogp)
  c_n <- pick_col(nms, COLS$n); c_eaf <- pick_col(nms, COLS$eaf)

  if (is.na(c_beta) || is.na(c_se)) stop("Missing required beta/se columns. Found: ", paste(nms, collapse = ", "))
  if (is.na(c_snp) && (is.na(c_pos) || is.na(c_ea) || is.na(c_oa))) {
    stop("Missing SNP/ID column and cannot infer variant fields. Found: ", paste(nms, collapse = ", "))
  }

  variant_id <- if (!is.na(c_snp)) as.character(dt[[c_snp]]) else NA_character_
  parsed <- if (!is.na(c_snp)) parse_variant_id(variant_id) else data.table(
    chr_from_id = rep(NA_character_, nrow(dt)),
    pos_from_id = rep(NA_integer_, nrow(dt)),
    other_from_id = rep(NA_character_, nrow(dt)),
    effect_from_id = rep(NA_character_, nrow(dt))
  )

  chr_vec <- if (!is.na(c_chr)) as.character(dt[[c_chr]]) else parsed$chr_from_id
  pos_vec <- if (!is.na(c_pos)) suppressWarnings(as.integer(dt[[c_pos]])) else parsed$pos_from_id
  ea_vec <- if (!is.na(c_ea)) safe_upper(dt[[c_ea]]) else parsed$effect_from_id
  oa_vec <- if (!is.na(c_oa)) safe_upper(dt[[c_oa]]) else parsed$other_from_id

  out <- data.table(
    gene_symbol = gene, ancestry = ancestry, source_file = source_file,
    chr = chr_vec, pos_hg38 = pos_vec,
    effect_allele = ea_vec, other_allele = oa_vec,
    beta = suppressWarnings(as.numeric(dt[[c_beta]])),
    se = suppressWarnings(as.numeric(dt[[c_se]]))
  )

  out[, SNP := if (!is.na(c_snp)) variant_id else paste0("chr", chr, ":", pos_hg38, ":", other_allele, ":", effect_allele)]
  if (!is.na(c_p)) {
    out[, pval := suppressWarnings(as.numeric(dt[[c_p]]))]
  } else if (!is.na(c_mlogp)) {
    out[, pval := 10^(-suppressWarnings(as.numeric(dt[[c_mlogp]])))]
  } else {
    out[, pval := NA_real_]
  }

  out[, samplesize := if (!is.na(c_n)) suppressWarnings(as.numeric(dt[[c_n]])) else NA_real_]
  out[, eaf := if (!is.na(c_eaf)) suppressWarnings(as.numeric(dt[[c_eaf]])) else NA_real_]
  out[, exposure := paste(gene_symbol, ancestry, sep = "__")]
  out[, id.exposure := exposure]
  out[, beta.exposure := beta]
  out[, se.exposure := se]
  out[, effect_allele.exposure := effect_allele]
  out[, other_allele.exposure := other_allele]
  out[, eaf.exposure := eaf]
  out[, pval.exposure := pval]
  out[, samplesize.exposure := samplesize]

  out <- out[!is.na(chr) & !is.na(pos_hg38) & !is.na(beta) & !is.na(se)]
  out <- out[effect_allele %in% c("A", "C", "G", "T") & other_allele %in% c("A", "C", "G", "T")]
  if (nrow(out)) { setorder(out, SNP, pval); out <- unique(out, by = c("gene_symbol", "ancestry", "SNP")) }
  out
}

# ---------- optional cis coordinate ----------
load_gene_coords <- function() {
  for (p in c("results/qc/gene_coordinates_hg38.tsv", "data/reference/gene_coordinates_hg38.tsv", "reference/gene_coordinates_hg38.tsv")) {
    if (file.exists(p)) {
      x <- fread(p); names(x) <- tolower(names(x))
      if (all(c("gene_symbol", "chr", "start", "end") %in% names(x))) {
        x[, gene_symbol := safe_upper(gene_symbol)]
        x[, chr := as.character(chr)]
        x[, start := as.integer(start)]
        x[, end := as.integer(end)]
        msg("[INFO] Loaded gene coordinates: ", p)
        return(x)
      }
    }
  }
  msg("[WARN] No gene coordinate table found. Cis filter will be skipped.")
  NULL
}
gene_coords <- load_gene_coords()

apply_filters <- function(dt, gene) {
  if (!nrow(dt)) return(dt)
  dt[, is_cis := NA]
  if (!no_cis_filter && !is.null(gene_coords)) {
    coord <- gene_coords[gene_symbol == gene]
    if (nrow(coord)) {
      coord <- coord[1]
      gene_chr <- gsub("^chr", "", as.character(coord$chr), ignore.case = TRUE)
      vchr <- gsub("^chr", "", as.character(dt$chr), ignore.case = TRUE)
      dt[, is_cis := vchr == gene_chr & pos_hg38 >= (coord$start - cis_window) & pos_hg38 <= (coord$end + cis_window)]
      dt <- dt[is_cis == TRUE]
    } else {
      msg("[WARN] No coordinate for ", gene, "; cis filter skipped for this gene.")
    }
  }
  if (!no_p_filter && !keep_all) dt <- dt[!is.na(pval) & pval < p_threshold]
  if (!nrow(dt)) return(dt)
  dt[, F_stat := (beta / se)^2]
  if (!no_f_filter && !keep_all) dt <- dt[!is.na(F_stat) & F_stat > 10]
  dt
}

# ---------- tar helpers: stream read, no full untar ----------
find_gene_tar <- function(gene, ancestry) {
  d <- file.path(rawdir, ancestry)
  if (!dir.exists(d)) return(character(0))
  fs <- list.files(d, pattern = paste0("^", gene, "_.*\\.tar$"), full.names = TRUE)
  fs <- fs[!grepl("\\.synapse_download_", fs)]
  if (!is.null(allowed_source_files)) fs <- fs[basename(fs) %in% allowed_source_files]
  if (!length(fs)) return(character(0))
  info <- file.info(fs)
  fs[!is.na(info$size) & info$size > 0]
}

stage_tar_to_local <- function(tar_file, ancestry) {
  if (!copy_to_local) return(tar_file)
  d <- file.path(tmp_root, "staged_tar", batch_id, ancestry)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  local <- file.path(d, basename(tar_file))
  if (!file.exists(local) || file.info(local)$size != file.info(tar_file)$size) {
    msg("[STAGE] ", basename(tar_file))
    ok <- file.copy(tar_file, local, overwrite = TRUE)
    if (!ok) { msg("[WARN] local copy failed; using Drive file"); return(tar_file) }
  }
  local
}

list_tar_files <- function(tar_file) {
  x <- tryCatch(system2("tar", args = c("-tf", shQuote(tar_file)), stdout = TRUE, stderr = TRUE), error = function(e) character(0))
  x <- x[!grepl("^tar:", x)]
  x[x != ""]
}

detect_summary_file_in_tar <- function(tar_file) {
  fs <- list_tar_files(tar_file)
  if (!length(fs)) stop("No files listed in tar: ", tar_file)
  cand <- fs[grepl("\\.(tsv|txt|csv|gz|bgz)$", fs, ignore.case = TRUE)]
  if (!length(cand)) cand <- fs
  bad <- grepl("readme|manifest|metadata|index|html|json|yaml|yml", basename(cand), ignore.case = TRUE)
  if (any(!bad)) cand <- cand[!bad]
  pref <- cand[grepl("sumstat|summary|pqtl|gwas|assoc|variant|txt|tsv", basename(cand), ignore.case = TRUE)]
  if (length(pref)) cand <- pref
  cand[1]
}

build_tar_cmd <- function(tar_file, inner) {
  cmd <- paste("tar -xOf", shQuote(tar_file), shQuote(inner))
  if (grepl("\\.(gz|bgz)$", inner, ignore.case = TRUE)) cmd <- paste(cmd, "| gzip -dc")
  cmd
}

read_pqtl_from_tar_fast <- function(tar_file) {
  inner <- detect_summary_file_in_tar(tar_file)
  cmd <- build_tar_cmd(tar_file, inner)
  header <- fread(cmd = cmd, nrows = 0, showProgress = FALSE, data.table = TRUE)
  select_cols <- choose_needed_columns(names(header))
  if (!length(select_cols)) stop("Could not identify useful columns in ", basename(tar_file))
  if (max_file_lines > 0) {
    data_rows <- max(0L, max_file_lines - 1L)
    return(fread(cmd = cmd, select = select_cols, nrows = data_rows, showProgress = FALSE, data.table = TRUE))
  }
  fread(cmd = cmd, select = select_cols, showProgress = FALSE, data.table = TRUE)
}

process_one <- function(gene, ancestry) {
  status0 <- data.table(batch_id = batch_id, gene_symbol = gene, ancestry = ancestry, status = "pending", n_raw_rows = NA_integer_, n_filtered_rows = NA_integer_, tar_file = NA_character_, output_file = NA_character_, message = NA_character_)
  tars <- find_gene_tar(gene, ancestry)
  if (!length(tars)) { status0[, `:=`(status = "missing_tar", message = "No valid .tar file found")]; return(list(data = NULL, status = status0)) }
  res_list <- list(); st_list <- list()
  for (tar0 in tars) {
    panel <- safe_label(tools::file_path_sans_ext(basename(tar0)))
    out <- file.path(per_gene_dir, paste0(gene, "__", ancestry, "__", panel, ".tsv"))
    st <- copy(status0); st[, `:=`(tar_file = tar0, output_file = out)]
    if (!force && file.exists(out) && file.info(out)$size > 0) {
      existing <- tryCatch(fread(out, showProgress = FALSE), error = function(e) NULL)
      if (!is.null(existing)) {
        st[, `:=`(status = "completed_existing", n_filtered_rows = nrow(existing), message = "Existing output loaded")]
        res_list[[length(res_list) + 1]] <- existing; st_list[[length(st_list) + 1]] <- st; next
      }
    }
    msg("[PROCESS] ", batch_id, " | ", gene, " | ", ancestry, " | ", basename(tar0))
    one <- tryCatch({
      tar <- stage_tar_to_local(tar0, ancestry)
      raw <- read_pqtl_from_tar_fast(tar); st[, n_raw_rows := nrow(raw)]
      std <- standardize_pqtl(raw, gene, ancestry, basename(tar0)); rm(raw); gc(FALSE)
      filt <- apply_filters(std, gene); rm(std); gc(FALSE)
      st[, n_filtered_rows := nrow(filt)]
      if (!nrow(filt)) { st[, `:=`(status = "no_variants_after_filter", message = "No variants after filters")]; write_empty_exposure(out); NULL }
      else { fwrite(filt, out, sep = "\t"); st[, `:=`(status = "completed", message = "Saved")]; filt }
    }, error = function(e) {
      st[, `:=`(status = "failed", message = conditionMessage(e))]
      msg("[ERROR] ", gene, " | ", ancestry, " | ", basename(tar0), " | ", conditionMessage(e)); NULL
    })
    if (!is.null(one)) res_list[[length(res_list) + 1]] <- one
    st_list[[length(st_list) + 1]] <- st
    gc(FALSE)
  }
  list(data = if (length(res_list)) rbindlist(res_list, fill = TRUE) else NULL, status = rbindlist(st_list, fill = TRUE))
}

# ---------- main ----------
all_status <- list(); batch_result_files <- character(0)
for (gene in target_genes) {
  msg("[GENE] ", gene)
  gene_results <- list()
  if (eur_first && all(c("EUR", "EAS") %in% ancestries)) {
    eur <- process_one(gene, "EUR")
    if (!is.null(eur$data) && nrow(eur$data)) gene_results[[length(gene_results) + 1]] <- eur$data
    all_status[[length(all_status) + 1]] <- eur$status; fwrite(rbindlist(all_status, fill = TRUE), gene_status_out, sep = "\t")
    if (is.null(eur$data) || !nrow(eur$data)) {
      msg("[SKIP EAS] ", gene, " has no EUR instrument after filters")
      all_status[[length(all_status) + 1]] <- data.table(batch_id = batch_id, gene_symbol = gene, ancestry = "EAS", status = "skipped_no_eur_instrument", n_raw_rows = NA_integer_, n_filtered_rows = 0L, tar_file = NA_character_, output_file = NA_character_, message = "EUR-first mode: EAS skipped because EUR had no instruments")
      fwrite(rbindlist(all_status, fill = TRUE), gene_status_out, sep = "\t")
    } else {
      eas <- process_one(gene, "EAS")
      if (!is.null(eas$data) && nrow(eas$data)) gene_results[[length(gene_results) + 1]] <- eas$data
      all_status[[length(all_status) + 1]] <- eas$status; fwrite(rbindlist(all_status, fill = TRUE), gene_status_out, sep = "\t")
    }
  } else {
    for (anc in ancestries) {
      x <- process_one(gene, anc)
      if (!is.null(x$data) && nrow(x$data)) gene_results[[length(gene_results) + 1]] <- x$data
      all_status[[length(all_status) + 1]] <- x$status; fwrite(rbindlist(all_status, fill = TRUE), gene_status_out, sep = "\t")
    }
  }
  if (length(gene_results)) {
    gdt <- rbindlist(gene_results, fill = TRUE)
    gout <- file.path(per_gene_dir, paste0(gene, "__combined.tsv"))
    fwrite(gdt, gout, sep = "\t"); batch_result_files <- c(batch_result_files, gout)
    rm(gdt)
  }
  rm(gene_results); gc(FALSE)
}

msg("[INFO] Combining batch outputs")
valid_files <- batch_result_files[file.exists(batch_result_files) & file.info(batch_result_files)$size > 0]
if (!length(valid_files)) {
  write_empty_exposure(batch_out)
} else {
  lst <- list()
  for (f in valid_files) { x <- fread(f, showProgress = FALSE); if (nrow(x)) lst[[length(lst) + 1]] <- x; rm(x); gc(FALSE) }
  if (!length(lst)) write_empty_exposure(batch_out) else {
    bdt <- rbindlist(lst, fill = TRUE)
    ord <- intersect(c("gene_symbol", "ancestry", "pval", "SNP"), names(bdt)); if (length(ord)) setorderv(bdt, ord)
    fwrite(bdt, batch_out, sep = "\t")
    msg("[INFO] Final batch rows: ", nrow(bdt))
    rm(bdt)
  }
}

fwrite(rbindlist(all_status, fill = TRUE), gene_status_out, sep = "\t")
staged <- file.path(tmp_root, "staged_tar", batch_id)
if (dir.exists(staged)) unlink(staged, recursive = TRUE, force = TRUE)
msg("[INFO] Completed batch: ", batch_id)
msg("[INFO] Output: ", batch_out)
