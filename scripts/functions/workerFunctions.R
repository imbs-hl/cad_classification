
#' Change compression of Oxford format files from bz2 to gz
#'
#' This function changes the compression of Oxford format files from bz2 to gz.
#' Additionally, the chromosome ID is extracted and written in first column.
#'
#' @param bz.file  [\code{string}]\cr
#'                 Input file name of bz2 compressed gen file.
#' @param gz.file  [\code{string}]\cr
#'                 Output file name of gz compressed gen file.
#'
#' @return The captured system2 output.
#'
bz2gz <- function(bz.file, gz.file) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(bz.file, add = collection)
  checkmate::assertDirectory(dirname(gz.file), add = collection)

  checkmate::reportAssertions(collection)

  dir.create(dirname(gz.file), recursive = TRUE)
  sys_out <- system2(command = "bunzip2",
                     args = c("-c", "<", bz.file,
                              # Write in first column the correct chromosome number
                              "|", "awk '{split($1, a, \":\"); $1=a[1]; print}'",
                              # gzip the output
                              "|", "gzip", "-c", ">", gz.file),
                     stdout = TRUE,
                     stderr = TRUE)
  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }
  return(sys_out)
}

#' Convert Oxford IMPUTE2 files into \code{PLINK} binary files
#'
#' @param gz.file          [\code{string}]\cr
#'                         Filepath of the gen file.
#' @param sample.file      [\code{string}]\cr
#'                         Filepath of the sample file.
#' @param pheno            [\code{string}]\cr
#'                         Name of phenotype column.
#' @param call.threshold   [\code{sting} or \code{numeric}]\cr
#'                         Maxiumum distance of genotype probabilities to allow
#'                         a hard call. If "random" the genotype is called
#'                         according to the given genotype probabilities.
#' @param exec             [\code{string}]\cr
#'                         Name of the \code{PLINK} executable.
#'
#' @return A list with the dirname of the converted files and the captured output
#' from \code{\link{system2}}.
#'
gz2bed <- function(gz.file, sample.file, pheno, call.threshold, exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(gz.file, add = collection)
  checkmate::assertFile(sample.file, add = collection)
  checkmate::assertString(pheno,
                          na.ok = FALSE, null.ok = FALSE,
                          add = collection)
  checkmate::assertNumber(call.threshold,
                          lower = 0, upper = 0.499,
                          null.ok = FALSE, na.ok = FALSE, finite = TRUE,
                          add = collection)

  checkmate::reportAssertions(collection)

  sys_out <- system2(command = exec,
                     args = c("--gen", gz.file,
                              "--sample", sample.file,
                              "--oxford-pheno-name", pheno,
                              "--hard-call-threshold", call.threshold,
                              "--make-bed",
                              "--seed", as.integer(runif(1, min = 0, max = .Machine$integer.max)),
                              "--out", out <- sub("^([^.]*).*", "\\1", gz.file)),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  return(list(dir = dirname(gz.file),
              out = out))
}

#' Copy sample files and fix coding of phenotype
#'
#' @param file [\code{string}]\cr
#'             Filepath of sample file to copy.
#'
#' @return Filepath of new sample file
#'
copy_sample_files <- function(file) {

  checkmate::assertFile(file)

  cohort <- basename(dirname(file))
  new_file <- file.path(proc_dir, cohort, "SampleInfos-cohort.txt")
  sys_out <- system2(command = "awk",
                     args = c(sprintf("'{if(FNR>2){$1=$1\"-%s\"; $2=$2\"-%s\"; $6=($6-1)}; print $0}'", cohort, cohort),
                              file, ">", new_file),
                     stdout = TRUE,
                     stderr = TRUE)
  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }
  return(new_file)
}

#' Merge multiple \code{PLINK} files from a list
#'
#' @param mergelist [\code{string}]\cr
#'                  A filepath to a list of \code{PLINK} files to merge.
#' @param output    [\code{string}]\cr
#'                  Filepath prefix of \code{PLINK} output.
#' @param exec      [\code{string}]\cr
#'                  Name of the \code{PLINK} executable.
#'
#' @return A filepath as \code{string} to the merged \code{PLINK} output.
#'
merge_plink <- function(mergelist, output, exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(mergelist, add = collection)
  checkmate::assertDirectory(dirname(output), add = collection)

  checkmate::reportAssertions(collection)

  sys_out <- system2(command = exec,
                     args = c("--merge-list", mergelist,
                              "--make-bed",
                              "--out", output),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  return(output)
}

#' Calculate SNP and sample statistics in IMPUTE2 files
#'
#' @param gen.file     [\code{string}]\cr
#'                     Filepath of IMPUTE2 gen file in Oxford format.
#' @param sample.file  [\code{string}]\cr
#'                     Filepath of sample file in Oxford format.
#' @param exec         [\code{string}]\cr
#'                     \code{qctool} executable.
#'
#' @return A list with the system output, the file paths of newly created files,
#' a \code{\link{data.table}} with SNP statistics and a \code{\link{data.table}}
#' with sample statistics
#'
gen_stats <- function(gen.file, sample.file, exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(gen.file, add = collection)
  checkmate::assertFile(sample.file, add = collection)

  checkmate::reportAssertions(collection)

  snp_stats_file <- paste(tools::file_path_sans_ext(gen.file, compression = TRUE),
                          "snp-stats", sep = ".")
  sample_stats_file <- paste(tools::file_path_sans_ext(gen.file, compression = TRUE),
                             "sample-stats", sep = ".")
  sample_file <- paste(tools::file_path_sans_ext(gen.file, compression = TRUE),
                       "sample", sep = ".")

  sys_out <- system2(command = exec,
                     args = c("-g", gen.file,
                              "-s", sample.file,
                              "-snp-stats", snp_stats_file,
                              "-sample-stats", sample_stats_file,
                              "-os", sample_file),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  snp_stats <- data.table::fread(snp_stats_file,
                                 drop = c("chromosome",
                                          "A_allele",
                                          "B_allele",
                                          "AA",
                                          "AB",
                                          "BB",
                                          "AA_calls",
                                          "AB_calls",
                                          "BB_calls",
                                          "missing"))
  snp_stats[, c("STUDY", "CHR") := list(basename(dirname(snp_stats_file)),
                                        gsub(pattern = "(\\d+):.*",
                                             replacement = "\\1",
                                             x = SNPID))]
  data.table::setnames(snp_stats,
                       c("position", "minor_allele", "major_allele", "missing_calls", "information"),
                       c("POS", "MINORALLELE", "MAJORALLELE", "MISSINGCALLS", "INFO"))
  snp_stats[, SNPID := NULL]

  sample_stats <- data.table::fread(sample_stats_file)
  data.table::setnames(sample_stats,
                       c("ID_1", "ID_2", "missing", "heterozygosity"),
                       c("ID1", "ID2", "MISSING", "HETEROZYGOSITY"))

  return(list(sys_out = sys_out,
              SNPSTATSFILE = snp_stats_file,
              SAMPLESTATFILE = sample_stats_file,
              SAMPLEFILE = sample_file,
              SNPSTATS = snp_stats,
              SAMPLESTATS = sample_stats))
}

#' Find samples to exclude based on heterozygosity
#'
#' @param plink.file.prefix      [\code{string}]\cr
#'                      Prefix path of \code{PLINK} binary files.
#' @param het.fac       [\code{number}]\cr
#'                      Maximum factor of standard deviations exceeding the
#'                      mean heterozygosity of study population for a sample
#'                      heterozygosity.
#' @param keep          [\code{string} or \code{character}]\cr
#'                      Filepath to list of sample IDs or character vector
#'                      with sample IDs to keep.
#' @param remove        [\code{string} or \code{character}]\cr
#'                      Filepath to list of sample IDs or character vector
#'                      with sample IDs to remove.
#' @param extract       [\code{string} or \code{character}]\cr
#'                      Filepath to list of variants or character vector
#'                      with variants to keep.
#' @param exclude       [\code{string} or \code{character}]\cr
#'                      Filepath to list of variants or character vector
#'                      with variants to remove.
#' @param exec          [\code{string}]\cr
#'                      Name of the \code{PLINK} executable.
#'
#' @return A \code{\link{list}} with the name of the input BED file, the filename
#' with sample identifiers to exclude and the heterozygosity statistics as a
#' \code{\link{data.table}}.
#'
plink_heterozygosity <- function(plink.file.prefix, het.fac,
                                 keep, remove,
                                 extract, exclude,
                                 exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  checkmate::assertNumber(het.fac,
                          lower = 1, upper = Inf,
                          finite = TRUE, null.ok = FALSE, na.ok = FALSE,
                          add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }

  checkmate::reportAssertions(collection)

  het_stats_file_prefix <- sprintf("%s_stats", plink.file.prefix)
  het_stats_file <- sprintf("%s.het", het_stats_file_prefix)
  het_fail_file <- sprintf("%s.het.fail", het_stats_file_prefix)

  sys_out <- system2(command = exec,
                     args = c("--bfile", plink.file.prefix,
                              "--het",
                              keep,
                              remove,
                              exclude,
                              extract,
                              "--out", het_stats_file_prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  het_stats <- data.table::fread(het_stats_file, head = TRUE)
  setnames(het_stats, c("FID", "IID", "OBS_HOM", "EXP_HOM", "N_SNPS", "F"))

  # Calculate heterozygosity rate
  het_stats[, HET_RATE := (N_SNPS - OBS_HOM)/N_SNPS]

  # Find samples with high deviation from mean heterozygosity
  write.table(het_stats[HET_RATE < mean(HET_RATE) - het.fac*sd(HET_RATE) |
                          HET_RATE > mean(HET_RATE) + het.fac*sd(HET_RATE),
                        .SD,
                        .SDcols = c("FID", "IID")],
              file = het_fail_file,
              quote = FALSE,
              col.names = FALSE,
              row.names = FALSE)

  return(list(BEDFILE = plink.file.prefix,
              HETFAILFILE = het_fail_file,
              HETEROZYGOSITYSTATS = het_stats))
}

#' Call to \code{PLINK} QC functions
#'
#' @param plink.file.prefix [\code{string}]\cr
#'                          Prefix path of \code{PLINK} binary files.
#' @param keep              [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to keep.
#' @param remove            [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to remove.
#' @param extract           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to keep.
#' @param exclude           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to remove.
#' @param hwe.pvalue        [\code{number}]\cr
#'                          P value threshold for HWE test with \code{midp}
#'                          correction. Variants with a P value below this
#'                          threshold are removed.
#' @param geno              [\code{number}]\cr
#'                          Maximum missing fraction in variants.
#' @param mind              [\code{number}]\cr
#'                          Maximum missing fraction in samples.
#' @param maf               [\code{number}]\cr
#'                          Minimium minor allele frequency of variants.
#' @param qc.iteration      [\code{int}]\cr
#'                          Number of QC iteration for output file name.
#' @param exec              [\code{string}]\cr
#'                          Name of the \code{PLINK} executable.
#'
#' @return A \code{\link{list}} with the prefix path of the newly created
#' \code{PLINK} binary files and from log extracted numbers of samples and variants
#' removed or kept during QC.
#'
plink_qc <- function(plink.file.prefix,
                     keep, remove,
                     extract, exclude,
                     hwe.pvalue = 0, geno = 1, mind = 1, maf = 0,
                     qc.iteration = 1,
                     exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }
  checkmate::assertNumber(hwe.pvalue,
                          lower = -Inf, upper = 0.5,
                          finite = TRUE, null.ok = FALSE, na.ok = FALSE,
                          add = collection)
  checkmate::assertNumber(geno,
                          lower = 0, upper = 1,
                          finite = TRUE, null.ok = FALSE, na.ok = FALSE,
                          add = collection)
  checkmate::assertNumber(mind,
                          lower = 0, upper = 1,
                          finite = TRUE, null.ok = FALSE, na.ok = FALSE,
                          add = collection)
  checkmate::assertNumber(maf,
                          lower = 0, upper = 1,
                          finite = TRUE, null.ok = FALSE, na.ok = FALSE,
                          add = collection)
  checkmate::assertInt(qc.iteration,
                       lower = 1,
                       null.ok = FALSE, na.ok = FALSE,
                       add = collection)

  checkmate::reportAssertions(collection)

  out_file <- sprintf("%s_QC%s", gsub(sprintf("_QC%d", qc.iteration - 1), "", plink.file.prefix), qc.iteration)

  # Apply all filters
  sys_out <- system2(command = exec,
                     args = c("--bfile", plink.file.prefix,
                              "--make-bed",
                              "--geno", geno,
                              "--mind", mind,
                              ifelse(maf > 0, sprintf("--maf %f", maf), ""),
                              "--snps-only", "just-acgt",
                              "--hwe", hwe.pvalue, "midp",
                              keep,
                              remove,
                              exclude,
                              extract,
                              "--out", out_file),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  return(list(BEDFILE = out_file,
              HWE = gsub(pattern = "[^0-9]",
                         x = sys_out[max(grep("--hwe", sys_out))],
                         replacement = ""),
              KEEP = gsub(pattern = "[^0-9]",
                          x = sys_out[max(grep("--keep", sys_out))],
                          replacement = ""),
              REMOVE = gsub(pattern = "[^0-9]",
                            x = sys_out[max(grep("--remove", sys_out))],
                            replacement = ""),
              EXTRACT = gsub(pattern = "[^0-9]",
                             x = sys_out[max(grep("--keep", sys_out))],
                             replacement = ""),
              EXCLUDE = gsub(pattern = "[^0-9]",
                             x = sys_out[max(grep("--remove", sys_out))],
                             replacement = ""),
              MIND = gsub(pattern = "[^0-9]",
                          x = sys_out[max(grep("--mind", sys_out))],
                          replacement = ""),
              GENO = gsub(pattern = "[^0-9]",
                          x = sys_out[max(grep("--geno", sys_out))],
                          replacement = ""),
              MAF = gsub(pattern = "[^0-9]",
                         x = sys_out[max(grep("minor allele", sys_out))],
                         replacement = "")))
}

#' LD based variant pruning on a single chromosome using \code{PLINK}
#'
#' @param plink.file.prefix [\code{string}]\cr
#'                          Prefix path of \code{PLINK} binary files.
#' @param chr               [\code{int}]\cr
#'                          Specifies the chromosome to analyse.
#' @param keep              [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to keep.
#' @param remove            [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to remove.
#' @param extract           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to keep.
#' @param exclude           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to remove.
#' @param ld.window.size    [\code{int}]\cr
#'                          Specifies the window size.
#' @param ld.step.size      [\code{int}]\cr
#'                          Specifies the step size.
#' @param ld.threshold      [\code{number}]\cr
#'                          Specifies the LD threshold. One of two variants
#'                          exceeding this threshold is removed.
#' @param exec              [\code{string}]\cr
#'                          Name of the \code{PLINK} executable.
#'
#' @return A \code{\link{list}} with the file names of the list of variants to
#' keep or remove.
#'
plink_ld_pruning <- function(plink.file.prefix, chr,
                             keep, remove,
                             extract, exclude,
                             ld.window.size, ld.step.size, ld.threshold, exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  checkmate::assertInt(chr,
                       lower = 1, upper = 25,
                       na.ok = FALSE, null.ok = FALSE,
                       add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }
  checkmate::assertInt(ld.window.size, lower = 1, upper = Inf,
                       na.ok = FALSE, null.ok = FALSE,
                       add = collection)
  checkmate::assertInt(ld.step.size, lower = 1, upper = Inf,
                       null.ok = FALSE, na.ok = FALSE,
                       add = collection)
  checkmate::assertNumber(ld.threshold, lower = 0, upper = 1,
                          na.ok = FALSE, finite = TRUE, null.ok = FALSE,
                          add = collection)

  checkmate::reportAssertions(collection)

  ld_file_prefix <- sprintf("%s.chr%02d", plink.file.prefix, chr)
  ld_prune_in_file <- sprintf("%s.prune.in", ld_file_prefix)
  ld_prune_out_file <- sprintf("%s.prune.out", ld_file_prefix)

  sys_out <- system2(command = exec,
                     args = c("--bfile", plink.file.prefix,
                              "--chr", chr,
                              keep,
                              remove,
                              exclude,
                              extract,
                              sprintf("--indep-pairwise %d %d %f", ld.window.size, ld.step.size, ld.threshold),
                              "--out", ld_file_prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  return(list(PLINKPREFIX = plink.file.prefix,
              PRUNEIN = ld_prune_in_file,
              PRUNEOUT = ld_prune_out_file))
}

#' Calculate IBD statistics using \code{PLINK}
#'
#' @param plink.file.prefix [\code{string}]\cr
#'                          Prefix path of \code{PLINK} binary files.
#' @param pi.hat            [\code{number}]\cr
#'                          Exclude sample pairs below this threshold.
#' @param keep              [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to keep.
#' @param remove            [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to remove.
#' @param extract           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to keep.
#' @param exclude           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to remove.
#' @param exec              [\code{string}]\cr
#'                          Name of the \code{PLINK} executable.
#'
#' @return A \code{\link{data.table}} with IBD statistics.
#'
plink_ibd <- function(plink.file.prefix, pi.hat,
                      keep, remove,
                      extract, exclude,
                      exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  checkmate::assertNumber(pi.hat, lower = 0, upper = 1,
                          na.ok = FALSE, finite = TRUE, null.ok = FALSE,
                          add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }

  checkmate::reportAssertions(collection)

  ibd_stats_file_prefix <- sprintf("%s_stats", plink.file.prefix)
  ibd_stats_file <- sprintf("%s.genome", ibd_stats_file_prefix)
  ibd_fail_file <- sprintf("%s.ibd.fail", ibd_stats_file_prefix)

  sys_out <- system2(command = exec,
                     args = c("--bfile", plink.file.prefix,
                              keep,
                              remove,
                              exclude,
                              extract,
                              "--genome",
                              "--min", pi.hat,
                              "--out", ibd_stats_file_prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  ibd_stats <- data.table::fread(ibd_stats_file)

  data.table::fwrite(ibd_stats[, .(FID2, IID2)], file = ibd_fail_file,
                     row.names = FALSE, col.names = FALSE, sep = "\t")

  return(list(IBDSTATS = ibd_stats,
              IBDFAILFILE = ibd_fail_file))
}

#' Run logistic association test with \code{PLINK}
#'
#' This function runs a sex adjusted logistic association test using \code{PLINK}.
#'
#' @param plink.file.prefix [\code{string}]\cr
#'                          Prefix path of \code{PLINK} binary files.
#' @param chr               [\code{int}]\cr
#'                          Specifies the chromosome to analyse.
#' @param keep              [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to keep.
#' @param remove            [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to remove.
#' @param extract           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to keep.
#' @param exclude           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to remove.
#' @param exec              [\code{string}]\cr
#'                          Name of the \code{PLINK} executable.
#'
#' @return A list with the output file path of \code{PLINK} and a \code{\link{data.table}} containing the results.
#'
plink_gwa <- function(plink.file.prefix, chr,
                      keep, remove,
                      extract, exclude,
                      out.dir, exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  checkmate::assertInt(chr,
                       lower = 1, upper = 25,
                       na.ok = FALSE, null.ok = FALSE,
                       add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }
  if (missing(out.dir)) {
    out.dir <- dirname(plink.file.prefix)
  } else {
    checkmate::assertDirectory(out.dir, add = collection)
  }

  checkmate::reportAssertions(collection)

  gwa_file_prefix <- file.path(out.dir,
                               sprintf("%s.chr%02d",
                                       basename(plink.file.prefix),
                                       chr))
  gwa_file <- sprintf("%s.assoc.logistic", gwa_file_prefix)

  sys_out <- system2(command = exec,
                     args = c("--bfile", plink.file.prefix,
                              "--logistic", "sex",
                              "--chr", chr,
                              keep,
                              remove,
                              exclude,
                              extract,
                              "--out", gwa_file_prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  gwa <- data.table::fread(gwa_file)

  return(list(GWAFILE = gwa_file,
              GWA = gwa))

}

#' Create list of intersecting SNPs from multiple \code{PLINK} \code{bim} files
#'
#' @param bim.files [\code{character}]\cr
#'                  Character vector of \code{PLINK} \code{bim} files.
#' @param out.file  [\code{string}]\cr
#'                  Filepath to list of intersecting SNPs.
#'
#' @return Nothing.
#'
find_snp_intersection <- function(bim.files, out.file) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertCharacter(bim.files,
                             any.missing = FALSE, all.missing = FALSE, null.ok = FALSE,
                             min.len = 2, unique = TRUE,
                             add = collection)
  for (bim.file in bim.files) {
    checkmate::assertFile(bim.file, add = collection)
  }
  checkmate::assertDirectory(dirname(out.file), add = collection)

  checkmate::reportAssertions(collection)

  n <- length(bim.files)

  sys_out <- system2(command = "grep",
                     args = c("-Fx",
                              "-f", bim.files[1], bim.files[2],
                              ifelse(n > 2, paste(sprintf("| grep -Fx -f - %s", bim.files[3:n]), collapse = " "), ""),
                              "|",
                              "cut", "-f2",
                              ">", out.file),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

}

#' LD based variant pruning on a single chromosome using \code{PLINK}
#'
#' @param plink.file.prefix  [\code{string}]\cr
#'                           Prefix path of \code{PLINK} binary files.
#' @param chr                [\code{int}]\cr
#'                           Specifies the chromosome to analyse.
#' @param keep               [\code{string} or \code{character}]\cr
#'                           Filepath to list of sample IDs or character vector
#'                           with sample IDs to keep.
#' @param remove             [\code{string} or \code{character}]\cr
#'                           Filepath to list of sample IDs or character vector
#'                           with sample IDs to remove.
#' @param extract            [\code{string} or \code{character}]\cr
#'                           Filepath to list of variants or character vector
#'                           with variants to keep.
#' @param exclude            [\code{string} or \code{character}]\cr
#'                           Filepath to list of variants or character vector
#'                           with variants to remove.
#' @param ld.window.size     [\code{int}]\cr
#'                           Specifies the window size.
#' @param ld.step.size       [\code{int}]\cr
#'                           Specifies the step size.
#' @param ld.threshold       [\code{number}]\cr
#'                           Specifies the LD threshold. One of two variants
#'                           exceeding this threshold is removed.
#' @param exec               [\code{string}]\cr
#'                           Name of the \code{PLINK} executable.
#'
#' @return A \code{\link{list}} with the file names of the list of variants to
#' keep or remove.
#'
plink_ld_pruning <- function(plink.file.prefix, chr,
                             keep, remove,
                             extract, exclude,
                             ld.window.size, ld.step.size, ld.threshold, exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  checkmate::assertInt(chr,
                       lower = 1, upper = 25,
                       na.ok = FALSE, null.ok = FALSE,
                       add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }
  checkmate::assertInt(ld.window.size, lower = 1, upper = Inf,
                       na.ok = FALSE, null.ok = FALSE,
                       add = collection)
  checkmate::assertInt(ld.step.size, lower = 1, upper = Inf,
                       null.ok = FALSE, na.ok = FALSE,
                       add = collection)
  checkmate::assertNumber(ld.threshold, lower = 0, upper = 1,
                          na.ok = FALSE, finite = TRUE, null.ok = FALSE,
                          add = collection)

  checkmate::reportAssertions(collection)

  ld_file_prefix <- sprintf("%s.chr%02d", plink.file.prefix, chr)
  ld_prune_in_file <- sprintf("%s.prune.in", ld_file_prefix)
  ld_prune_out_file <- sprintf("%s.prune.out", ld_file_prefix)

  sys_out <- system2(command = exec,
                     args = c("--bfile", plink.file.prefix,
                              "--chr", chr,
                              extract,
                              exclude,
                              remove,
                              keep,
                              sprintf("--indep-pairwise %d %d %f", ld.window.size, ld.step.size, ld.threshold),
                              "--out", ld_file_prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  ld_prune_in <- data.table::fread(ld_prune_in_file)

  return(list(PLINKPREFIX = plink.file.prefix,
              PRUNEINFILE = ld_prune_in_file,
              PRUNEOUTFILE = ld_prune_out_file,
              PRUNEIN = ld_prune_in))
}

#' LD-based result clumping
#'
#' @param plink.file.prefix [\code{string}]\cr
#'                          Prefix path of \code{PLINK} binary files.
#' @param chr               [\code{int}]\cr
#'                          Specifies the chromosome to analyse.
#' @param keep              [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to keep.
#' @param remove            [\code{string} or \code{character}]\cr
#'                          Filepath to list of sample IDs or character vector
#'                          with sample IDs to remove.
#' @param extract           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to keep.
#' @param exclude           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or character vector
#'                          with variants to remove.
#' @param gwa.result.file   [\code{string}]\cr
#'                          File path of GWA result file in \code{PLINK} format.
#' @param p1                [\code{number}]\cr
#'                          Index variants are chosen greedily starting with the
#'                          lowest p-value. Variants with p-values which meet
#'                          the \code{p1} threshold, but have already been
#'                          assigned to another clump, do not start their own
#'                          clumps.
#' @param p2                [\code{number}]\cr
#'                          Variants that have association p-value smaller than
#'                          0.01 are assigned to some index variant's clump
#'                          (unless they have been previously been assigned to
#'                          another clump.
#' @param kb                [\code{number}]\cr
#'                          Sites which are less than \code{kb} away from an
#'                          index variant are assigned to some index variant's
#'                          clump (unless they have been previously been
#'                          assigned to another clump.
#' @param r2                [\code{number}]\cr
#'                          Sites which have LD larger than \code{r2} are
#'                          assigned to some index variant's clump (unless they
#'                          have been previously been assigned to another clump.
#' @param exec              [\code{string}]\cr
#'                          Name of the \code{PLINK} executable.
#'
#' @return A \code{\link[data.table]{data.table}} object containing the clumping
#' results.
#'
plink_clumping <- function(plink.file.prefix, gwa.result.file, chr,
                           keep, remove,
                           extract, exclude,
                           p1, p2, kb, r2, exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  checkmate::assertFile(gwa.result.file, add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }
  checkmate::assertInt(chr,
                       lower = 1, upper = 22,
                       na.ok = FALSE, null.ok = FALSE,
                       add = collection)
  checkmate::assertNumber(p1,
                          lower = 0, upper = 1,
                          na.ok = FALSE, null.ok = FALSE, finite = TRUE,
                          add = collection)
  checkmate::assertNumber(p2,
                          lower = 0, upper = 1,
                          na.ok = FALSE, null.ok = FALSE, finite = TRUE,
                          add = collection)
  checkmate::assertNumber(kb,
                          lower = 0, upper = Inf,
                          na.ok = FALSE, null.ok = FALSE, finite = TRUE,
                          add = collection)
  checkmate::assertNumber(r2,
                          lower = 0, upper = 1,
                          na.ok = FALSE, null.ok = FALSE, finite = TRUE,
                          add = collection)

  checkmate::reportAssertions(collection)

  clump_file_prefix <- file.path(sprintf("%s.chr%02d",
                                         plink.file.prefix,
                                         chr))
  clump_file <- sprintf("%s.clumped", clump_file_prefix)

  sys_out <- system2(command = exec,
                     args = c("--bfile", plink.file.prefix,
                              "--chr", chr,
                              extract,
                              exclude,
                              remove,
                              keep,
                              "--clump", gwa.result.file,
                              "--clump-p1", p1,
                              "--clump-p2", p2,
                              "--clump-kb", kb,
                              "--clump-r2", r2,
                              "--out", clump_file_prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  res <- data.table::fread(clump_file)

  return(res)
}

#' Convert \code{IMPUTE2} files to dosages
#'
#' @param gen.file          [\code{string}]\cr
#'                          File path of \code{IMPUTE2} gen file.
#' @param plink.file.prefix [\code{string}]\cr
#'                          Prefix path of \code{PLINK} binary files.
#' @param out.prefix        [\code{string}]\cr
#'                          Prefix path of output files.
#' @param chr               [\code{int}]\cr
#'                          Specifies the chromosome to analyse.
#' @param keep              [\code{string} or \code{data.frame}]\cr
#'                          Filepath to list of sample IDs or \code{data.frame}
#'                          with sample IDs to keep.
#' @param remove            [\code{string} or \code{data.frame}]\cr
#'                          Filepath to list of sample IDs or \code{data.frame}
#'                          with sample IDs to remove.
#' @param extract           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or \code{character} vector
#'                          with variants to keep.
#' @param exclude           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or \code{character} vector
#'                          with variants to remove.
#' @param qctool.exec       [\code{string}]\cr
#'                          Name of the \code{QCTOOL} executable.
#' @param fcgene.exec       [\code{string}]\cr
#'                          Name of the \code{fcGENE} executable.
#'
#' @return A \code{\link[data.table]{data.table}} object containing the sample \code{ID}, its \code{STATUS} and the dosage of the alternative allele.
#'
gen2dosage <- function(gen.file, sample.file, plink.file.prefix,
                       out.prefix = plink.file.prefix,
                       chr,
                       keep, remove,
                       extract, exclude,
                       qctool.exec, fcgene.exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(gen.file, add = collection)
  checkmate::assertString(plink.file.prefix, add = collection)
  checkmate::assertString(out.prefix, add = collection)
  checkmate::assertFile(sample.file, add = collection)
  checkmate::assertFile(sprintf("%s.fam", plink.file.prefix), add = collection)
  checkmate::assertInt(chr,
                       lower = 1, upper = 22,
                       na.ok = FALSE, null.ok = FALSE,
                       add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("-incl-samples %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("-excl-samples %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("-incl-rsids %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("-excl-rsids %s", exclude)
  }

  checkmate::reportAssertions(collection)

  temp_gen_file <- sprintf("%s.chr%d.gen.gz", out.prefix, chr)
  fam_file <- sprintf("%s.fam", plink.file.prefix)
  pedinfo_file <- sprintf("%s.chr%d.pedinfo", out.prefix, chr)
  dosage_file <- sprintf("%s_dose.txt", out.prefix)
  affection_file <- sprintf("%s_affstat.txt", out.prefix)

  sys_out <- system2(command = qctool.exec,
                     args = c("-g", gen.file,
                              "-s", sample.file,
                              extract,
                              exclude,
                              remove,
                              keep,
                              "-og", temp_gen_file,
                              "-omit-chromosome"),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  num_extracted_snps <- as.numeric(gsub(".*\\(total (\\d+) snps\\).*", "\\1",
                                        grep("(total \\d+ snps)",
                                             sys_out,
                                             value = TRUE)))

  # Create pedinfo file
  fam <- data.table::fread(fam_file)
  data.table::setnames(fam, c("famid", "indid", "matid", "patid", "sex", "phenotype"))

  data.table::fwrite(fam, pedinfo_file, sep = "\t")

  # Read pedinfo file
  pedinfo <- data.table::fread(pedinfo_file, select = c("indid", "sex", "phenotype"))

  if (num_extracted_snps > 0) {
    # SNPs to create dosage data for are available
    sys_out <- system2(command = fcgene.exec,
                       args = c("--gens", temp_gen_file,
                                "--oformat", "r-dose",
                                "--out", out.prefix,
                                "--pedinfo", pedinfo_file),
                       stdout = TRUE,
                       stderr = TRUE)

    if (!is.null(attr(sys_out, "status"))) {
      stop(paste(sys_out, collapse = "\n"))
    }

    # Read dosage file
    dosage <- data.table::fread(dosage_file)
  } else {
    # No SNPs to create dosage data for are available
    # Create dummy result object
    dosage <- data.table::data.table(SMAPLE_ID = pedinfo$indid)
    data.table::fwrite(dosage, dosage_file, sep = " ")
  }

  res <- pedinfo[dosage, on = c("indid" = "SMAPLE_ID")]

  data.table::setnames(res, c("indid", "sex", "phenotype"), c("ID", "SEX", "STATUS"))

  # Change column type of SNP columns to numeric
  for (col in setdiff(colnames(res), c("ID", "STATUS", "SEX"))) {
    data.table::set(res, j = col, value = as.numeric(res[[col]]))
  }
  data.table::set(res, j = "SEX", value = factor(res[["SEX"]], levels = c(1, 2), labels = c("male", "female")))
  data.table::set(res, j = "STATUS", value = factor(res[["STATUS"]], levels = c(1, 2), labels = c("control", "case")))

  return(res)

}

#' Convert \code{PLINK} files to \code{GenABEL}
#'
#' @param plink.file.prefix [\code{string}]\cr
#'                          Prefix path of \code{PLINK} binary files.
#' @param keep              [\code{string} or \code{data.frame}]\cr
#'                          Filepath to list of sample IDs or \code{data.frame}
#'                          with sample IDs to keep.
#' @param remove            [\code{string} or \code{data.frame}]\cr
#'                          Filepath to list of sample IDs or \code{data.frame}
#'                          with sample IDs to remove.
#' @param extract           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or \code{character} vector
#'                          with variants to keep.
#' @param exclude           [\code{string} or \code{character}]\cr
#'                          Filepath to list of variants or \code{character} vector
#'                          with variants to remove.
#' @param fcgene.exec       [\code{string}]\cr
#'                          Name of the \code{fcGENE} executable.
#'
#' @return A \code{\link[data.table]{data.table}} object containing the sample \code{ID}, its \code{STATUS} and the dosage of the alternative allele.
#'
plink2gwaa <- function(plink.file.prefix,
                       out.prefix = plink.file.prefix,
                       keep, remove,
                       extract, exclude,
                       plink.exec,
                       fcgene.exec) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertString(plink.file.prefix, add = collection)
  checkmate::assertString(out.prefix, add = collection)
  checkmate::assertFile(sprintf("%s.fam", plink.file.prefix), add = collection)
  checkmate::assertFile(sprintf("%s.bim", plink.file.prefix), add = collection)
  checkmate::assertFile(sprintf("%s.bed", plink.file.prefix), add = collection)
  if (missing(keep)) {
    keep <- ""
  } else {
    if (checkmate::testDataFrame(keep)) {
      keep_file <- tempfile()
      data.table::fwrite(keep, keep_file, sep = "\t")
      keep <- keep_file
    }
    checkmate::assertFile(keep, add = collection)
    keep <- sprintf("--keep %s", keep)
  }
  if (missing(remove)) {
    remove <- ""
  } else {
    if (checkmate::testDataFrame(remove)) {
      remove_file <- tempfile()
      data.table::fwrite(remove, remove_file, sep = "\t")
      remove <- remove_file
    }
    checkmate::assertFile(remove, add = collection)
    remove <- sprintf("--remove %s", remove)
  }
  if (missing(extract)) {
    extract <- ""
  } else {
    if (checkmate::testCharacter(extract, min.len = 2)) {
      extract_file <- tempfile()
      writeLines(extract, extract_file)
      extract <- extract_file
    }
    checkmate::assertFile(extract, add = collection)
    extract <- sprintf("--extract %s", extract)
  }
  if (missing(exclude)) {
    exclude <- ""
  } else {
    if (checkmate::testCharacter(exclude, min.len = 2)) {
      exclude_file <- tempfile()
      writeLines(exclude, exclude_file)
      exclude <- exclude_file
    }
    checkmate::assertFile(exclude, add = collection)
    exclude <- sprintf("--exclude %s", exclude)
  }

  checkmate::reportAssertions(collection)

  temp_plink_prefix <- tempfile()
  raw_file <- sprintf("%s.raw", out.prefix)
  phe_file <- sprintf("%s_phe0.dat", out.prefix)

  sys_out <- system2(command = plink.exec,
                     args = c("--bfile", plink.file.prefix,
                              extract,
                              exclude,
                              remove,
                              keep,
                              "--make-bed",
                              "--out", temp_plink_prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  # Convert to GenABEL data
  sys_out <- system2(command = fcgene.exec,
                     args = c("--bed", sprintf("%s.bed", temp_plink_prefix),
                              "--bim", sprintf("%s.bim", temp_plink_prefix),
                              "--fam", sprintf("%s.fam", temp_plink_prefix),
                              "--oformat", "genabel",
                              "--out", out.prefix),
                     stdout = TRUE,
                     stderr = TRUE)

  if (!is.null(attr(sys_out, "status"))) {
    stop(paste(sys_out, collapse = "\n"))
  }

  # Load GenABEL object
  res <- GenABEL::load.gwaa.data(phenofile = phe_file, genofile = raw_file, force = TRUE)
  GenABEL::phdata(res)$disease <- factor(GenABEL::phdata(res)$disease)

  return(res)

}

#' Calculate \code{ranger} importance using a gwaa rds file
#'
#' @param gwaa.rds.file [\code{string}]\cr
#'                      File path of rds file of a gwaa object to load.
#' @param importance    [\code{string}]\cr
#'                      Importance to compute.
#' @param threads       [\code{integer}]\cr
#'                      How many CPUs to use.
#'
#' @return A \code{ranger} object.
#'
calc_ranger_importance <- function(gwaa.rds.file, mtry.perc, ...) {

  collection <- checkmate::makeAssertCollection()

  checkmate::assertFile(gwaa.rds.file, add = collection)

  checkmate::reportAssertions(collection)

  gwaa <- readRDS(gwaa.rds.file)

  num_snps <- nsnps(gwaa)

  mtry <- ceiling(sqrt(num_snps))
  if (!missing(mtry.perc)) {
    mtry <- ceiling(num_snps * mtry.perc)
  }

  ranger <- ranger::ranger(disease~0, data = gwaa, mtry = mtry, ...)

  return(ranger)

}

bcftools_query <- function(vcf.file, query, regions.file, output.file, print.header = TRUE, bcftools.exe = "bcftools") {

  sys_out <- system2(
    command = bcftools.exe,
    args = c(
      "query",
      ifelse(print.header, "-H", ""),
      "-f", query,
      "-o", output.file,
      ifelse(missing(regions.file), "", sprintf("-R %s", regions.file)),
      vcf.file
    ), stdout = TRUE, stderr = TRUE)

  result <- data.table::fread(output.file)

  return(result)

}

bcftools_index <- function(vcf.file, output.file, threads = 1, bcftools.exe = "bcftools") {

  sys_out <- system2(
    command = bcftools.exe,
    args = c(
      "index",
      "-t",
      "--threads", threads,
      "-o", output.file,
      vcf.file
    ), stdout = TRUE, stderr = TRUE)

  return(sys_out)

}
