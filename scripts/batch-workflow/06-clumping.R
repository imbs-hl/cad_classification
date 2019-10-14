
source("init.R", chdir = TRUE)

## Load training fam ----
training_fam <- readRDS(training_fam_rds_file)

# LD clumping ----
clumping_reg <- load_or_create_registry("clumping", conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = clumping_reg))) {
  ids <- batchtools::batchMap(reg = clumping_reg,
                              fun = plink_clumping,
                              chr = 1:22,
                              more.args = list(plink.file.prefix = german_qc3_file_prefix,
                                               gwa.result.file = german_gwa_results_file,
                                               keep = training_fam[, .(FID, IID)],
                                               p1 = clumping_p1,
                                               p2 = clumping_p2,
                                               kb = clumping_kb,
                                               r2 = clumping_r2,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s.chr%02d.clumped",
                                      german_qc3_file_prefix, 1:22),
                              checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = clumping_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = clumping_reg))) {
  ids <- batchtools::findNotDone(reg = clumping_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = clumping_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_clumping"))
}

if (!batchtools::waitForJobs(reg = clumping_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(clumping_reg$file.dir)))
}

clumped_results <- rbindlist(batchtools::reduceResultsList(reg = clumping_reg))
saveRDS(clumped_results, clumped_results_rds)

writeLines(clumped_results$SNP, clumped_snps_file)
