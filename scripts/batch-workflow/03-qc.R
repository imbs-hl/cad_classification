
source("init.R", chdir = TRUE)

# Load rds files from previous steps ----
sample_files <- readRDS(sample_files_rds)
gen_files <- readRDS(gen_files_rds)
plink_files <- readRDS(plink_files_rds)

# Create SNP statistics ----
gen_stats_reg <- load_or_create_registry("genStats",
                                         packages = c("data.table", "tools"),
                                         conf.file = slurm_btconf_file)

params <- gen_files[sample_files, on = "STUDY"]
if (!nrow(ids <- batchtools::findJobs(reg = gen_stats_reg))) {
  ids <- batchtools::batchMap(reg = gen_stats_reg,
                              fun = gen_stats,
                              gen.file = params$GZFILE,
                              sample.file = params$SAMPLEFILE,
                              more.args = list(exec = qctool_exec))
}

missing_files <- which(sapply(paste(tools::file_path_sans_ext(params$GZFILE, compression = TRUE),
                                    "sample-stats", sep = "."), checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = gen_stats_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = gen_stats_reg))) {
  ids <- batchtools::findNotDone(reg = gen_stats_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = gen_stats_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_gen_stats"))
}

if (!batchtools::waitForJobs(reg = gen_stats_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(gen_stats_reg$file.dir)))
}

snp_stats <- data.table::rbindlist(batchtools::reduceResultsList(reg = gen_stats_reg,
                                                                 fun = function(job, res, ...) {
                                                                   res$SNPSTATS
                                                                 }))

saveRDS(snp_stats, snp_stats_rds)

# Write info filtered SNP lists ----
n_studies <- gen_files[, length(unique(STUDY))]
snp_stats[INFO > info_threshold, # select SNPs with minimum imputation quality
          .(count = .N),
          by = RSID][count == n_studies, # select SNPs present in all studies
                     writeLines(RSID, overlapping_snps_file)]

# QC1 only SNP level ----
qc <- 1
qc1_reg <- load_or_create_registry(sprintf("QC%d", qc), conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = qc1_reg))) {
  ids <- batchtools::batchMap(reg = qc1_reg,
                              fun = plink_qc,
                              plink.file.prefix = c(german_file_prefix,
                                                    wtccc_file_prefix,
                                                    cg_file_prefix),
                              more.args = list(geno = geno,
                                               mind = 1,
                                               hwe.pvalue = hwe_pvalue,
                                               maf = maf,
                                               extract = overlapping_snps_file,
                                               qc.iteration = qc,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s.bed", c(german_qc1_file_prefix,
                                                  wtccc_qc1_file_prefix,
                                                  cg_qc1_file_prefix)), checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = qc1_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = qc1_reg))) {
  ids <- batchtools::findNotDone(reg = qc1_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = qc1_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 10,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_qc1"))
}

if (!batchtools::waitForJobs(reg = qc1_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(qc1_reg$file.dir)))
}

# heterozygosity ----
heterozygosity_reg <- load_or_create_registry("heterozygosity",
                                              packages = c("data.table"),
                                              conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = heterozygosity_reg))) {
  ids <- batchtools::batchMap(reg = heterozygosity_reg,
                              fun = plink_heterozygosity,
                              plink.file.prefix = c(german_qc1_file_prefix,
                                                    wtccc_qc1_file_prefix,
                                                    cg_qc1_file_prefix),
                              more.args = list(het.fac = het_fac,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s_stats.het.fail", c(german_qc1_file_prefix,
                                                             wtccc_qc1_file_prefix,
                                                             cg_qc1_file_prefix)),
                              checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = heterozygosity_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = heterozygosity_reg))) {
  ids <- batchtools::findNotDone(reg = heterozygosity_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = heterozygosity_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 10,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_heterozygosity"))
}

if (!batchtools::waitForJobs(reg = heterozygosity_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(heterozygosity_reg$file.dir)))
}

het_fail_files <- unlist(batchtools::reduceResultsList(reg = heterozygosity_reg,
                                                       fun = function(job, res, ...) res$HETFAILFILE))

# QC2 sample and SNP level ----
qc <- 2
qc2_reg <- load_or_create_registry(sprintf("QC%d", qc), conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = qc2_reg))) {
  ids <- batchtools::batchMap(reg = qc2_reg,
                              fun = plink_qc,
                              plink.file.prefix = c(german_qc1_file_prefix,
                                                    wtccc_qc1_file_prefix,
                                                    cg_qc1_file_prefix),
                              remove = het_fail_files,
                              more.args = list(geno = geno,
                                               mind = mind,
                                               hwe.pvalue = hwe_pvalue,
                                               maf = maf,
                                               extract = overlapping_snps_file,
                                               qc.iteration = qc,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s.bed", c(german_qc2_file_prefix,
                                                  wtccc_qc2_file_prefix,
                                                  cg_qc2_file_prefix)), checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = qc2_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = qc2_reg))) {
  ids <- batchtools::findNotDone(reg = qc2_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = qc2_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 10,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_qc2"))
}

if (!batchtools::waitForJobs(reg = qc2_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(qc2_reg$file.dir)))
}

# IBD ----
# IBD calculation requires approximately independent SNPs
ibd_ld_reg <- load_or_create_registry("IBDLD", conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = ibd_ld_reg))) {
  ids <- batchtools::batchMap(reg = ibd_ld_reg,
                              fun = plink_ld_pruning,
                              plink.file.prefix = c(rep(german_qc2_file_prefix, 22),
                                                    rep(wtccc_qc2_file_prefix, 22),
                                                    rep(cg_qc2_file_prefix, 22)),
                              chr = rep(1:22, 3),
                              more.args = list(ld.window.size = ibd_ld_window_size,
                                               ld.step.size = ibd_ld_step_size,
                                               ld.threshold = ibd_ld_threshold,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s.chr%02d.prune.in",
                                      files = c(rep(german_qc2_file_prefix, 22),
                                                rep(wtccc_qc2_file_prefix, 22),
                                                rep(cg_qc2_file_prefix, 22)),
                                      chromosomes = 1:22),
                              checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = ibd_ld_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = ibd_ld_reg))) {
  ids <- batchtools::findNotDone(reg = ibd_ld_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = ibd_ld_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 20000,
                                          partition = "fast", walltime = 30,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_ibd_ld"))
}

if (!batchtools::waitForJobs(reg = ibd_ld_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(ibd_ld_reg$file.dir)))
}

# Combine prune in files
ibd_prune_in_files_per_chromosome <- batchtools::reduceResultsDataTable(reg = ibd_ld_reg)
ibd_prune_in_files_per_chromosome[, DATASET := basename(dirname(PLINKPREFIX))]
ibd_prune_in_files <- data.table::data.table(DATASET = c("GERMAN",
                                                         "WTCCC",
                                                         "CG"),
                                             PRUNEINFILE = c(german_ibd_prune_in_file,
                                                             wtccc_ibd_prune_in_file,
                                                             cg_ibd_prune_in_file))

ibd_prune_in_files[, file.remove(PRUNEINFILE), by = DATASET]
ibd_prune_in_files_per_chromosome[ibd_prune_in_files,
                                  file.append(PRUNEINFILE, PRUNEIN),
                                  on = .(DATASET)]

# IBD calculation
ibd_reg <- load_or_create_registry("IBD", conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = ibd_reg))) {
  ids <- batchtools::batchMap(reg = ibd_reg,
                              fun = plink_ibd,
                              plink.file.prefix = c(german_qc2_file_prefix,
                                                    wtccc_qc2_file_prefix,
                                                    cg_qc2_file_prefix),
                              extract = c(german_ibd_prune_in_file,
                                               wtccc_ibd_prune_in_file,
                                               cg_ibd_prune_in_file),
                              more.args = list(pi.hat = pi_hat,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s_stats.ibd.fail", c(german_qc2_file_prefix,
                                                             wtccc_qc2_file_prefix,
                                                             cg_qc2_file_prefix)),
                              checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = ibd_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = ibd_reg))) {
  ids <- batchtools::findNotDone(reg = ibd_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = ibd_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 10,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_ibd"))
}

if (!batchtools::waitForJobs(reg = ibd_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(ibd_reg$file.dir)))
}

ibd_fail_files <- unlist(batchtools::reduceResultsList(reg = ibd_reg,
                                                       fun = function(job, res, ...) res$IBDFAILFILE))

# Remove IBD failing samples
qc <- 3
qc3_reg <- load_or_create_registry(sprintf("QC%d", qc), conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = qc3_reg))) {
  ids <- batchtools::batchMap(reg = qc3_reg,
                              fun = plink_qc,
                              plink.file.prefix = c(german_qc2_file_prefix,
                                                    wtccc_qc2_file_prefix,
                                                    cg_qc2_file_prefix),
                              remove = ibd_fail_files,
                              more.args = list(geno = geno,
                                               mind = mind,
                                               hwe.pvalue = hwe_pvalue,
                                               maf = maf,
                                               qc.iteration = qc,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s.bed", c(german_qc3_file_prefix,
                                                  wtccc_qc3_file_prefix,
                                                  cg_qc3_file_prefix)), checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = qc3_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = qc3_reg))) {
  ids <- batchtools::findNotDone(reg = qc3_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = qc3_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_qc3"))
}

if (!batchtools::waitForJobs(reg = qc3_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(qc3_reg$file.dir)))
}

# Overlapping SNPs after QC ----
overlap_reg <- load_or_create_registry("PostQCOverlap", conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = overlap_reg))) {
  ids <- batchtools::batchMap(reg = overlap_reg,
                              fun = find_snp_intersection,
                              out.file = post_qc_overlapping_snps_file,
                              more.args = list(bim.files = sprintf("%s.bim",
                                                                   c(german_qc3_file_prefix,
                                                                     wtccc_qc3_file_prefix,
                                                                     cg_qc3_file_prefix))))
}

missing_files <- which(sapply(sprintf("%s.bed", c(german_qc3_file_prefix,
                                                  wtccc_qc3_file_prefix,
                                                  cg_qc3_file_prefix)), checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = overlap_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = overlap_reg))) {
  ids <- batchtools::findNotDone(reg = overlap_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = overlap_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 10,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_post_qc_overlap"))
}

if (!batchtools::waitForJobs(reg = overlap_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(overlap_reg$file.dir)))
}
