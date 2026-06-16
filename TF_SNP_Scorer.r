#!/usr/bin/env Rscript

# ----------------------- EDIT THESE SETTINGS -----------------------
#this could be made into a config file eventually, for now just leaving like this
cfg <- list(
  input_tsv = "",
  jaspar_pfm = ,        #this is your input 
  fasta = "",  #genome build used.
  out_dir = ".../jaspar_pwm_ref_alt",    
  output_prefix = "",
  rsid = "",                     ##########CHANGE THIS TO THE VARIANT YOU WANT TO PLOT.
  filter_col = "",
  filter_value = "",
  filter_regex = FALSE,
  max_variants = 0,
  rank_tsv = "",
  rank_col = "exact_h5_shap",
  rank_n = 0,
  motif_regex = "",
  flank = 30,
  top_n_plot = 20,
  top_per_variant = 5,
  min_report_pwm_score = 0,
  overlap_variant_only = TRUE,
  chr_col = "chrom",
  pos_col = "pos",
  ref_col = "ref",
  alt_col = "alt",
  rsid_col = "rsid",
  variant_id_col = "variant_id",
  plot_width = 11,
  png_dpi = 300,
  save_pdf = TRUE
)
# -------------------------------------------------------------------

parse_bool <- function(x) {
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

read_yaml_config <- function(path) {
  if (!nzchar(path)) return(list())
  if (!requireNamespace("yaml", quietly = TRUE)) {
    stop("Install the yaml R package to use --config.", call. = FALSE)
  }
  if (!file.exists(path)) stop("Missing config file: ", path, call. = FALSE)
  config <- yaml::read_yaml(path)
  if (is.null(config)) list() else config
}

list_to_vector <- function(x) {
  if (is.list(x) && is.null(names(x))) return(unlist(x, use.names = FALSE))
  x
}

merge_config <- function(cfg, config) {
  for (key in intersect(names(config), names(cfg))) {
    cfg[[key]] <- list_to_vector(config[[key]])
  }
  cfg
}

find_config_arg <- function(args) {
  eq <- grep("^--config=", args, value = TRUE)
  if (length(eq) > 0) return(sub("^--config=", "", eq[[1]]))

  idx <- match("--config", args)
  if (!is.na(idx)) {
    if (idx == length(args)) stop("--config requires a file path.", call. = FALSE)
    return(args[[idx + 1L]])
  }

  ""
}

drop_config_args <- function(args) {
  out <- character()
  skip_next <- FALSE
  for (arg in args) {
    if (skip_next) {
      skip_next <- FALSE
      next
    }
    if (identical(arg, "--config")) {
      skip_next <- TRUE
      next
    }
    if (startsWith(arg, "--config=")) next
    out <- c(out, arg)
  }
  out
}

parse_args <- function(cfg) {
  args <- commandArgs(trailingOnly = TRUE)
  cfg <- merge_config(cfg, read_yaml_config(find_config_arg(args)))
  args <- drop_config_args(args)

  for (arg in args) {
    if (!startsWith(arg, "--")) next
    pieces <- base::strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- gsub("-", "_", pieces[1])
    value <- if (length(pieces) > 1) paste(pieces[-1], collapse = "=") else "TRUE"

    if (key %in% c("flank", "top_n_plot", "max_variants", "rank_n", "top_per_variant", "png_dpi")) {
      cfg[[key]] <- suppressWarnings(as.integer(value))
    } else if (key %in% c("min_report_pwm_score", "plot_width")) {
      cfg[[key]] <- suppressWarnings(as.numeric(value))
    } else if (key %in% c("overlap_variant_only", "filter_regex", "save_pdf")) {
      cfg[[key]] <- parse_bool(value)
    } else if (key %in% names(cfg)) {
      cfg[[key]] <- value
    }
  }
  cfg
}

required_packages <- c(
  "ggplot2",
  "dplyr",
  "readr",
  "stringr",
  "tidyr",
  "tibble",
  "Rsamtools",
  "Biostrings",
  "GenomicRanges",
  "IRanges"
)

missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install missing packages first: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Rsamtools)
  library(Biostrings)
  library(GenomicRanges)
  library(IRanges)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tibble)
})

cfg <- parse_args(cfg)

if (!file.exists(cfg$input_tsv)) stop("Missing input TSV: ", cfg$input_tsv, call. = FALSE)
if (!nzchar(cfg$jaspar_pfm)) stop("Set cfg$jaspar_pfm to a JASPAR PFM/TRANSFAC file.", call. = FALSE)
if (!file.exists(cfg$jaspar_pfm)) stop("Missing JASPAR PFM: ", cfg$jaspar_pfm, call. = FALSE)
if (nzchar(cfg$rank_tsv) && !file.exists(cfg$rank_tsv)) stop("Missing rank TSV: ", cfg$rank_tsv, call. = FALSE)
if (!toupper(cfg$fasta) %in% c("BSGENOME", "BSGENOME.HG38") && !file.exists(cfg$fasta)) {
  stop("Missing FASTA: ", cfg$fasta, call. = FALSE)
}
if (toupper(cfg$fasta) %in% c("BSGENOME", "BSGENOME.HG38") && !requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
  stop("cfg$fasta uses BSGENOME.HG38, but package BSgenome.Hsapiens.UCSC.hg38 is not installed.", call. = FALSE)
}
if (is.na(cfg$png_dpi) || cfg$png_dpi < 72) stop("png_dpi must be at least 72.", call. = FALSE)
if (!is.finite(cfg$plot_width) || cfg$plot_width <= 0) stop("plot_width must be positive.", call. = FALSE)
if (cfg$flank < 1) stop("flank must be positive.", call. = FALSE)
if (cfg$top_n_plot < 1) stop("top_n_plot must be positive.", call. = FALSE)
if (cfg$top_per_variant < 1) stop("top_per_variant must be positive.", call. = FALSE)

dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x

save_pub_plot <- function(plot, png_path, width, height) {
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = cfg$png_dpi, bg = "white", limitsize = FALSE)
  paths <- c(png = png_path)
  if (isTRUE(cfg$save_pdf)) {
    pdf_path <- sub("\\.png$", ".pdf", png_path)
    ggplot2::ggsave(pdf_path, plot, width = width, height = height, device = "pdf", bg = "white", limitsize = FALSE)
    paths <- c(paths, pdf = pdf_path)
  }
  paths
}

write_tsv_atomic <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  if (file.exists(tmp)) unlink(tmp)
  readr::write_tsv(x, tmp)
  if (file.exists(path)) unlink(path)
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop("Cannot replace output TSV: ", path, call. = FALSE)
  }
  invisible(path)
}

save_empty_plot <- function(path, title, subtitle, height = 4.5) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No motif hits passed the current filters.", size = 4) +
    labs(title = title, subtitle = subtitle) +
    xlim(-1, 1) +
    ylim(-1, 1) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey35")
    )
  save_pub_plot(p, path, width = cfg$plot_width, height = height)
}

read_jaspar_simple_pfm <- function(path) {
  lines <- readLines(path, warn = FALSE)
  motifs <- list()
  current_id <- NULL
  current_name <- NULL
  current <- list()

  flush_current <- function() {
    if (is.null(current_id)) return(NULL)
    bases <- c("A", "C", "G", "T")
    if (!all(bases %in% names(current))) return(NULL)
    mat <- do.call(rbind, current[bases])
    rownames(mat) <- bases
    list(id = current_id, name = current_name, counts = mat)
  }

  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next

    if (startsWith(line, ">")) {
      old <- flush_current()
      if (!is.null(old)) motifs[[length(motifs) + 1L]] <- old

      header <- sub("^>", "", line)
      pieces <- base::strsplit(header, "\\s+")[[1]]
      current_id <- pieces[1]
      current_name <- if (length(pieces) > 1) paste(pieces[-1], collapse = " ") else pieces[1]
      current <- list()
      next
    }

    if (grepl("^[ACGT]\\s+", line)) {
      base <- substr(line, 1, 1)
      nums <- gregexpr("-?[0-9]+\\.?[0-9]*", line, perl = TRUE)
      vals <- regmatches(line, nums)[[1]]
      vals <- as.numeric(vals)
      current[[base]] <- vals
    }
  }

  old <- flush_current()
  if (!is.null(old)) motifs[[length(motifs) + 1L]] <- old

  if (length(motifs) == 0) {
    stop("No motifs parsed from simple JASPAR PFM format.", call. = FALSE)
  }

  motifs
}

read_jaspar_transfac <- function(path) {
  lines <- readLines(path, warn = FALSE)
  motifs <- list()
  current_id <- NULL
  current_name <- NULL
  rows <- list()
  in_matrix <- FALSE

  flush_current <- function() {
    if (is.null(current_id) || length(rows) == 0) return(NULL)
    mat_pos <- do.call(rbind, rows)
    mat <- t(mat_pos)
    rownames(mat) <- c("A", "C", "G", "T")
    colnames(mat) <- seq_len(ncol(mat))
    list(id = current_id, name = current_name %||% current_id, counts = mat)
  }

  for (line in lines) {
    line <- trimws(line)
    if (!nzchar(line)) next

    if (startsWith(line, "AC ")) {
      current_id <- sub("^AC\\s+", "", line)
      current_name <- current_id
      rows <- list()
      in_matrix <- FALSE
      next
    }

    if (startsWith(line, "ID ")) {
      current_name <- sub("^ID\\s+", "", line)
      next
    }

    if (startsWith(line, "DE ")) {
      de <- sub("^DE\\s+", "", line)
      pieces <- base::strsplit(de, "\\s+")[[1]]
      if (length(pieces) >= 2) current_name <- pieces[2]
      next
    }

    if (grepl("^P[O0]\\s+", line)) {
      in_matrix <- TRUE
      next
    }

    if (line == "//") {
      old <- flush_current()
      if (!is.null(old)) motifs[[length(motifs) + 1L]] <- old
      current_id <- NULL
      current_name <- NULL
      rows <- list()
      in_matrix <- FALSE
      next
    }

    if (in_matrix && grepl("^[0-9]+\\s+", line)) {
      parts <- base::strsplit(line, "\\s+")[[1]]
      vals <- suppressWarnings(as.numeric(parts[2:5]))
      if (length(vals) == 4 && all(is.finite(vals))) {
        rows[[length(rows) + 1L]] <- vals
      }
    }
  }

  old <- flush_current()
  if (!is.null(old)) motifs[[length(motifs) + 1L]] <- old

  if (length(motifs) == 0) {
    stop("No motifs parsed from TRANSFAC format.", call. = FALSE)
  }

  motifs
}

read_motif_file <- function(path) {
  first_lines <- readLines(path, n = 50, warn = FALSE)
  if (any(startsWith(trimws(first_lines), "AC "))) {
    return(read_jaspar_transfac(path))
  }
  read_jaspar_simple_pfm(path)
}

counts_to_pwm <- function(counts, pseudocount = 0.8, bg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25)) {
  mat <- counts + pseudocount
  probs <- sweep(mat, 2, colSums(mat), "/")
  log2(sweep(probs, 1, bg[rownames(probs)], "/"))
}

score_kmer <- function(kmer, pwm) {
  codes <- utf8ToInt(kmer)
  idx <- match(codes, c(65L, 67L, 71L, 84L))
  if (length(idx) != ncol(pwm) || anyNA(idx)) return(NA_real_)
  sum(pwm[cbind(idx, seq_along(idx))])
}

scan_pwm <- function(seq_string, pwm, center_i, overlap_variant_only = TRUE) {
  seq_string <- toupper(seq_string)
  motif_width <- ncol(pwm)
  seq_width <- nchar(seq_string)
  if (motif_width > seq_width) {
    return(tibble(
      score = NA_real_,
      strand = NA_character_,
      start = NA_integer_,
      end = NA_integer_,
      rel_start = NA_integer_,
      rel_end = NA_integer_,
      variant_motif_pos = NA_integer_,
      kmer = NA_character_
    ))
  }

  rc_seq <- as.character(reverseComplement(DNAString(seq_string)))
  rows <- list()

  if (overlap_variant_only) {
    start_min <- max(1L, center_i - motif_width + 1L)
    start_max <- min(center_i, seq_width - motif_width + 1L)
    if (start_min > start_max) {
      return(tibble(
        score = NA_real_,
        strand = NA_character_,
        start = NA_integer_,
        end = NA_integer_,
        rel_start = NA_integer_,
        rel_end = NA_integer_,
        variant_motif_pos = NA_integer_,
        kmer = NA_character_
      ))
    }
    start_values <- seq.int(start_min, start_max)
  } else {
    start_values <- seq_len(seq_width - motif_width + 1L)
  }

  for (strand in c("+", "-")) {
    scan_seq <- if (strand == "+") seq_string else rc_seq
    for (start_i in start_values) {
      end_i <- start_i + motif_width - 1L
      kmer <- substr(scan_seq, start_i, end_i)
      score <- score_kmer(kmer, pwm)

      if (strand == "+") {
        genomic_start_i <- start_i
        genomic_end_i <- end_i
      } else {
        genomic_start_i <- seq_width - end_i + 1L
        genomic_end_i <- seq_width - start_i + 1L
      }

      variant_motif_pos <- if (center_i >= genomic_start_i && center_i <= genomic_end_i) {
        if (strand == "+") center_i - genomic_start_i + 1L else genomic_end_i - center_i + 1L
      } else {
        NA_integer_
      }

      rows[[length(rows) + 1L]] <- tibble(
        score = score,
        strand = strand,
        start = genomic_start_i,
        end = genomic_end_i,
        rel_start = genomic_start_i - center_i,
        rel_end = genomic_end_i - center_i,
        variant_motif_pos = variant_motif_pos,
        kmer = kmer
      )
    }
  }

  bind_rows(rows) %>%
    arrange(desc(score)) %>%
    slice_head(n = 1)
}

fetch_ref_alt_sequences <- function(df, cfg) {
  fasta <- NULL
  bsgenome <- NULL

  if (!toupper(cfg$fasta) %in% c("BSGENOME", "BSGENOME.HG38") && file.exists(cfg$fasta)) {
    fasta <- FaFile(cfg$fasta)
    open(fasta)
    on.exit(close(fasta), add = TRUE)
  }

  if (requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    bsgenome <- getExportedValue("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Hsapiens.UCSC.hg38")
  }

  out <- vector("list", nrow(df))
  for (i in seq_len(nrow(df))) {
    chr <- as.character(df[[cfg$chr_col]][i])
    pos <- as.integer(df[[cfg$pos_col]][i])
    ref <- toupper(as.character(df[[cfg$ref_col]][i]))
    alt <- toupper(as.character(df[[cfg$alt_col]][i]))

    if (nchar(ref) != 1 || nchar(alt) != 1 || !ref %in% c("A", "C", "G", "T") || !alt %in% c("A", "C", "G", "T")) {
      stop("Only single-base A/C/G/T variants are supported. Problem row: ", i, " (", chr, ":", pos, " ", ref, ">", alt, ")", call. = FALSE)
    }

    gr <- GRanges(
      seqnames = chr,
      ranges = IRanges(start = pos - cfg$flank, end = pos + cfg$flank)
    )

    ref_seq <- NULL

    if (!is.null(fasta)) {
      ref_seq <- tryCatch(
        toupper(as.character(getSeq(fasta, gr))),
        error = function(e) NULL
      )
    }

    if (is.null(ref_seq) && !is.null(bsgenome)) {
      ref_seq <- toupper(as.character(getSeq(bsgenome, gr)))
    }

    if (is.null(ref_seq)) {
      stop("Could not fetch sequence for ", chr, ":", pos, " from FASTA or BSgenome.", call. = FALSE)
    }
    center_i <- cfg$flank + 1L
    fasta_ref <- substr(ref_seq, center_i, center_i)

    if (!identical(fasta_ref, ref)) {
      warning("FASTA REF mismatch for ", chr, ":", pos, " expected ", ref, " observed ", fasta_ref)
    }

    alt_seq <- ref_seq
    substr(alt_seq, center_i, center_i) <- alt

    out[[i]] <- tibble(
      row_index = i,
      variant_id = if (cfg$variant_id_col %in% names(df)) as.character(df[[cfg$variant_id_col]][i]) else paste(chr, pos, ref, alt, sep = ":"),
      rsid = if (cfg$rsid_col %in% names(df)) as.character(df[[cfg$rsid_col]][i]) else NA_character_,
      chrom = chr,
      pos = pos,
      ref = ref,
      alt = alt,
      ref_seq = ref_seq,
      alt_seq = alt_seq,
      center_i = center_i
    )
  }

  bind_rows(out)
}

plot_pwm_scores <- function(score_df, path, top_n = 20) {
  plot_df <- score_df %>%
    mutate(max_pwm_score = pmax(ref_pwm_score, alt_pwm_score, na.rm = TRUE)) %>%
    filter(max_pwm_score >= cfg$min_report_pwm_score) %>%
    mutate(abs_delta = abs(delta_alt_minus_ref)) %>%
    arrange(desc(abs_delta), desc(max_pwm_score)) %>%
    slice_head(n = top_n) %>%
    mutate(
      variant_label = ifelse(is.na(rsid) | !nzchar(rsid), variant_id, rsid),
      motif_label = paste0(variant_label, " | ", motif_name, " (", motif_id, ")")
    ) %>%
    mutate(motif_label = stringr::str_wrap(motif_label, width = 82)) %>%
    mutate(motif_label = factor(motif_label, levels = rev(unique(motif_label))))

  if (nrow(plot_df) == 0) {
    return(save_empty_plot(
      path,
      "REF vs ALT JASPAR PWM scores",
      paste0("Motif filter: ", cfg$motif_regex, " | max allele PWM >= ", cfg$min_report_pwm_score)
    ))
  }

  long_df <- plot_df %>%
    select(variant_id, rsid, motif_label, ref_pwm_score, alt_pwm_score, delta_alt_minus_ref) %>%
    tidyr::pivot_longer(
      cols = c(ref_pwm_score, alt_pwm_score, delta_alt_minus_ref),
      names_to = "score_type",
      values_to = "score"
    ) %>%
    mutate(
      score_type = recode(
        score_type,
        ref_pwm_score = "REF PWM score",
        alt_pwm_score = "ALT PWM score",
        delta_alt_minus_ref = "ALT - REF"
      ),
      score_type = factor(score_type, levels = c("REF PWM score", "ALT PWM score", "ALT - REF"))
    )

  p <- ggplot(long_df, aes(x = score, y = motif_label, fill = score_type)) +
    geom_col(position = position_dodge2(width = 0.78, preserve = "single"), width = 0.68) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey35") +
    scale_fill_manual(values = c("REF PWM score" = "grey35", "ALT PWM score" = "#C7362F", "ALT - REF" = "#2F6FB3")) +
    labs(
      title = "REF vs ALT JASPAR PWM scores",
      subtitle = paste0(
        "Top ", nrow(plot_df), " variant-motif changes with max allele PWM >= ",
        cfg$min_report_pwm_score, " | motif filter: ", cfg$motif_regex
      ),
      x = "PWM log-odds score",
      y = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      plot.title.position = "plot"
    )

  save_pub_plot(p, path, width = cfg$plot_width, height = max(4.5, 0.32 * nrow(plot_df) + 1.6))
}

plot_best_hit_letters <- function(score_df, path) {
  best <- score_df %>%
    filter(is.finite(delta_alt_minus_ref), is.finite(ref_pwm_score), is.finite(alt_pwm_score)) %>%
    mutate(abs_delta = abs(delta_alt_minus_ref)) %>%
    arrange(desc(abs_delta)) %>%
    slice_head(n = 1)

  if (nrow(best) == 0) {
    return(save_empty_plot(path, "Best changed JASPAR motif letter match", "No finite REF/ALT PWM score pairs were available.", height = 3.2))
  }

  ref_seq <- best$ref_seq[[1]]
  alt_seq <- best$alt_seq[[1]]
  center_i <- best$center_i[[1]]
  ref_motif_start <- best$ref_rel_start[[1]] + center_i
  ref_motif_end <- best$ref_rel_end[[1]] + center_i
  alt_motif_start <- best$alt_rel_start[[1]] + center_i
  alt_motif_end <- best$alt_rel_end[[1]] + center_i

  seq_df <- bind_rows(
    tibble(allele = "REF", pos_i = seq_len(nchar(ref_seq)), base = base::strsplit(ref_seq, "", fixed = TRUE)[[1]]),
    tibble(allele = "ALT", pos_i = seq_len(nchar(alt_seq)), base = base::strsplit(alt_seq, "", fixed = TRUE)[[1]])
  ) %>%
    mutate(
      rel_pos = pos_i - center_i,
      is_variant = pos_i == center_i,
      in_pwm_hit = case_when(
        allele == "REF" ~ pos_i >= ref_motif_start & pos_i <= ref_motif_end,
        allele == "ALT" ~ pos_i >= alt_motif_start & pos_i <= alt_motif_end,
        TRUE ~ FALSE
      )
    )

  p <- ggplot(seq_df, aes(x = rel_pos, y = allele)) +
    geom_tile(aes(fill = in_pwm_hit), color = "white", linewidth = 0.25, height = 0.72) +
    geom_tile(data = filter(seq_df, is_variant), fill = "#F6D04D", color = "black", linewidth = 0.45, height = 0.72) +
    geom_text(aes(label = base), size = 3.2) +
    scale_fill_manual(values = c("FALSE" = "grey92", "TRUE" = "#DCE9F8"), guide = "none") +
    labs(
      title = paste0("Best changed JASPAR motif letter match: ", best$motif_name[[1]], " (", best$motif_id[[1]], ")"),
      subtitle = paste0(
        best$rsid[[1]], " ", best$chrom[[1]], ":", best$pos[[1]], " ",
        best$ref[[1]], ">", best$alt[[1]],
        " | REF=", round(best$ref_pwm_score[[1]], 3),
        " ALT=", round(best$alt_pwm_score[[1]], 3),
        " delta=", round(best$delta_alt_minus_ref[[1]], 3)
      ),
      x = "bp relative to variant",
      y = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title.position = "plot",
      axis.ticks.y = element_blank()
    )

  save_pub_plot(p, path, width = max(cfg$plot_width, 12), height = 3.2)
}

top_per_variant_scores <- function(score_df, n = 5) {
  score_df %>%
    filter(max_pwm_score >= cfg$min_report_pwm_score) %>%
    group_by(chrom, pos, ref, alt, rsid, variant_id) %>%
    arrange(desc(abs(delta_alt_minus_ref)), desc(max_pwm_score), .by_group = TRUE) %>%
    slice_head(n = n) %>%
    ungroup()
}

plot_top_per_variant_scores <- function(score_df, path, n = 5) {
  plot_df <- top_per_variant_scores(score_df, n = n) %>%
    mutate(
      variant_label = ifelse(is.na(rsid) | !nzchar(rsid), variant_id, rsid),
      motif_label = paste0(motif_name, " (", motif_id, ")"),
      row_label = paste0(variant_label, " | ", motif_label),
      row_label = stringr::str_wrap(row_label, width = 82),
      row_label = factor(row_label, levels = rev(unique(row_label)))
    )

  if (nrow(plot_df) == 0) {
    return(save_empty_plot(
      path,
      paste0("Top ", n, " JASPAR PWM changes per variant"),
      paste0("Motif filter: ", cfg$motif_regex, " | max allele PWM >= ", cfg$min_report_pwm_score)
    ))
  }

  long_df <- plot_df %>%
    select(variant_label, row_label, ref_pwm_score, alt_pwm_score, delta_alt_minus_ref) %>%
    tidyr::pivot_longer(
      cols = c(ref_pwm_score, alt_pwm_score, delta_alt_minus_ref),
      names_to = "score_type",
      values_to = "score"
    ) %>%
    mutate(
      score_type = recode(
        score_type,
        ref_pwm_score = "REF PWM score",
        alt_pwm_score = "ALT PWM score",
        delta_alt_minus_ref = "ALT - REF"
      ),
      score_type = factor(score_type, levels = c("REF PWM score", "ALT PWM score", "ALT - REF"))
    )

  p <- ggplot(long_df, aes(x = score, y = row_label, fill = score_type)) +
    geom_col(position = position_dodge2(width = 0.78, preserve = "single"), width = 0.68) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey35") +
    scale_fill_manual(values = c("REF PWM score" = "grey35", "ALT PWM score" = "#C7362F", "ALT - REF" = "#2F6FB3")) +
    labs(
      title = paste0("Top ", n, " JASPAR PWM changes per variant"),
      subtitle = paste0("Motif filter: ", cfg$motif_regex, " | max allele PWM >= ", cfg$min_report_pwm_score),
      x = "PWM log-odds score",
      y = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      legend.position = "top",
      plot.title.position = "plot"
    )

  save_pub_plot(p, path, width = max(cfg$plot_width, 13), height = max(5, 0.24 * nrow(long_df) + 1.8))
}

input_df <- readr::read_tsv(cfg$input_tsv, show_col_types = FALSE, progress = FALSE)
for (nm in names(input_df)) {
  input_df[[nm]] <- type.convert(input_df[[nm]], as.is = TRUE)
}

required_cols <- c(cfg$chr_col, cfg$pos_col, cfg$ref_col, cfg$alt_col)
missing_cols <- setdiff(required_cols, names(input_df))
if (length(missing_cols) > 0) {
  stop("Input TSV missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

if (nzchar(cfg$filter_col) || nzchar(cfg$filter_value)) {
  filter_cols <- base::strsplit(cfg$filter_col, ",", fixed = TRUE)[[1]]
  filter_values <- base::strsplit(cfg$filter_value, ",", fixed = TRUE)[[1]]
  filter_cols <- trimws(filter_cols)
  filter_values <- trimws(filter_values)

  if (length(filter_cols) != length(filter_values)) {
    stop("--filter-col and --filter-value must have the same comma-separated length.", call. = FALSE)
  }

  for (i in seq_along(filter_cols)) {
    col <- filter_cols[[i]]
    value <- filter_values[[i]]
    if (!col %in% names(input_df)) stop("Filter column missing from input TSV: ", col, call. = FALSE)

    if (cfg$filter_regex) {
      input_df <- input_df %>% filter(grepl(value, as.character(.data[[col]]), ignore.case = TRUE))
    } else {
      input_df <- input_df %>% filter(as.character(.data[[col]]) == value)
    }
  }

  if (nrow(input_df) == 0) {
    stop("No input rows remained after --filter-col/--filter-value filters.", call. = FALSE)
  }
}

if (nzchar(cfg$rsid) && !toupper(cfg$rsid) %in% c("ALL", "*", "NONE")) {
  if (!cfg$rsid_col %in% names(input_df)) stop("Cannot filter by rsid because column is missing: ", cfg$rsid_col, call. = FALSE)
  input_df <- input_df %>% filter(.data[[cfg$rsid_col]] == cfg$rsid)
  if (nrow(input_df) == 0) stop("No rows found for --rsid=", cfg$rsid, call. = FALSE)
}

if (nzchar(cfg$rank_tsv) && cfg$rank_n > 0) {
  if (!cfg$rsid_col %in% names(input_df)) stop("Cannot rank-filter because input rsid column is missing: ", cfg$rsid_col, call. = FALSE)
  rank_df <- readr::read_tsv(cfg$rank_tsv, show_col_types = FALSE, progress = FALSE)
  if (!cfg$rsid_col %in% names(rank_df)) stop("Rank TSV missing rsid column: ", cfg$rsid_col, call. = FALSE)
  if (!cfg$rank_col %in% names(rank_df)) stop("Rank TSV missing rank column: ", cfg$rank_col, call. = FALSE)

  top_rank <- rank_df %>%
    mutate(.rank_value = suppressWarnings(as.numeric(.data[[cfg$rank_col]]))) %>%
    filter(is.finite(.rank_value)) %>%
    arrange(desc(abs(.rank_value))) %>%
    distinct(.data[[cfg$rsid_col]], .keep_all = TRUE) %>%
    slice_head(n = cfg$rank_n) %>%
    pull(.data[[cfg$rsid_col]]) %>%
    as.character()

  input_df <- input_df %>% filter(as.character(.data[[cfg$rsid_col]]) %in% top_rank)
  if (nrow(input_df) == 0) stop("No input variants matched top ranked rsids from rank TSV.", call. = FALSE)
}

input_df <- input_df %>%
  distinct(.data[[cfg$chr_col]], .data[[cfg$pos_col]], .data[[cfg$ref_col]], .data[[cfg$alt_col]], .keep_all = TRUE)

if (cfg$max_variants > 0 && nrow(input_df) > cfg$max_variants) {
  input_df <- input_df %>% slice_head(n = cfg$max_variants)
}

motifs <- read_motif_file(cfg$jaspar_pfm)
motif_meta <- tibble(
  motif_index = seq_along(motifs),
  motif_id = vapply(motifs, `[[`, character(1), "id"),
  motif_name = vapply(motifs, `[[`, character(1), "name")
)

if (nzchar(cfg$motif_regex) && !toupper(cfg$motif_regex) %in% c("ALL", "*", "NONE")) {
  keep <- grepl(cfg$motif_regex, motif_meta$motif_name, ignore.case = TRUE) |
    grepl(cfg$motif_regex, motif_meta$motif_id, ignore.case = TRUE)
  motifs <- motifs[keep]
  motif_meta <- motif_meta[keep, , drop = FALSE]
}

if (length(motifs) == 0) {
  stop("No motifs matched --motif-regex=", cfg$motif_regex, call. = FALSE)
}

seq_df <- fetch_ref_alt_sequences(input_df, cfg)

score_rows <- list()
for (variant_i in seq_len(nrow(seq_df))) {
  v <- seq_df[variant_i, ]

  for (motif_i in seq_along(motifs)) {
    motif <- motifs[[motif_i]]
    pwm <- counts_to_pwm(motif$counts)

    ref_hit <- scan_pwm(v$ref_seq[[1]], pwm, v$center_i[[1]], cfg$overlap_variant_only)
    alt_hit <- scan_pwm(v$alt_seq[[1]], pwm, v$center_i[[1]], cfg$overlap_variant_only)

    score_rows[[length(score_rows) + 1L]] <- tibble(
      variant_id = v$variant_id[[1]],
      rsid = v$rsid[[1]],
      chrom = v$chrom[[1]],
      pos = v$pos[[1]],
      ref = v$ref[[1]],
      alt = v$alt[[1]],
      motif_id = motif$id,
      motif_name = motif$name,
      motif_width = ncol(pwm),
      ref_pwm_score = ref_hit$score[[1]],
      alt_pwm_score = alt_hit$score[[1]],
      delta_alt_minus_ref = alt_hit$score[[1]] - ref_hit$score[[1]],
      ref_strand = ref_hit$strand[[1]],
      alt_strand = alt_hit$strand[[1]],
      ref_rel_start = ref_hit$rel_start[[1]],
      ref_rel_end = ref_hit$rel_end[[1]],
      alt_rel_start = alt_hit$rel_start[[1]],
      alt_rel_end = alt_hit$rel_end[[1]],
      ref_variant_motif_pos = ref_hit$variant_motif_pos[[1]],
      alt_variant_motif_pos = alt_hit$variant_motif_pos[[1]],
      ref_best_kmer = ref_hit$kmer[[1]],
      alt_best_kmer = alt_hit$kmer[[1]],
      ref_seq = v$ref_seq[[1]],
      alt_seq = v$alt_seq[[1]],
      center_i = v$center_i[[1]]
    )
  }
}

score_df <- bind_rows(score_rows) %>%
  mutate(max_pwm_score = pmax(ref_pwm_score, alt_pwm_score, na.rm = TRUE)) %>%
  arrange(desc(abs(delta_alt_minus_ref)), desc(max_pwm_score))

prefix <- if (nzchar(cfg$output_prefix)) {
  cfg$output_prefix
} else if (nzchar(cfg$rsid) && !toupper(cfg$rsid) %in% c("ALL", "*", "NONE")) {
  cfg$rsid
} else {
  "all_variants"
}
prefix <- safe_name(prefix)

score_path <- file.path(cfg$out_dir, paste0(prefix, "_jaspar_pwm_ref_alt_scores.tsv"))
top_path <- file.path(cfg$out_dir, paste0(prefix, "_top_jaspar_pwm_ref_alt_changes.tsv"))
top_per_variant_path <- file.path(cfg$out_dir, paste0(prefix, "_top", cfg$top_per_variant, "_motifs_per_variant.tsv"))
motif_path <- file.path(cfg$out_dir, paste0(prefix, "_matched_jaspar_motifs.tsv"))
barplot_path <- file.path(cfg$out_dir, paste0(prefix, "_jaspar_pwm_ref_alt_scores.png"))
top_per_variant_plot_path <- file.path(cfg$out_dir, paste0(prefix, "_top", cfg$top_per_variant, "_motifs_per_variant.png"))
letters_path <- file.path(cfg$out_dir, paste0(prefix, "_top_changed_pwm_hit_letters.png"))

write_tsv_atomic(score_df, score_path)
write_tsv_atomic(
  score_df %>%
    filter(max_pwm_score >= cfg$min_report_pwm_score) %>%
    slice_head(n = max(cfg$top_n_plot, 100)),
  top_path
)
write_tsv_atomic(top_per_variant_scores(score_df, n = cfg$top_per_variant), top_per_variant_path)
write_tsv_atomic(motif_meta, motif_path)
pwm_plot_paths <- plot_pwm_scores(score_df, barplot_path, top_n = cfg$top_n_plot)
top_per_variant_plot_paths <- plot_top_per_variant_scores(score_df, top_per_variant_plot_path, n = cfg$top_per_variant)
letters_plot_paths <- plot_best_hit_letters(score_df, letters_path)

message("Wrote PWM score table: ", score_path)
message("Wrote top PWM changes table: ", top_path)
message("Wrote top motifs per variant table: ", top_per_variant_path)
message("Wrote matched motif table: ", motif_path)
message("Wrote PWM score plot(s): ", paste(pwm_plot_paths, collapse = ", "))
message("Wrote top motifs per variant plot(s): ", paste(top_per_variant_plot_paths, collapse = ", "))
message("Wrote top-hit letter plot(s): ", paste(letters_plot_paths, collapse = ", "))


