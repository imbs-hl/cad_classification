
source("init.R", chdir = TRUE)

## Load training fam ----
training_fam <- readRDS(training_fam_rds_file)

# LD pruning ----
pruning_reg <- load_or_create_registry("pruning", conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = pruning_reg))) {
  ids <- batchtools::batchMap(reg = pruning_reg,
                              fun = plink_ld_pruning,
                              chr = 1:22,
                              more.args = list(plink.file.prefix = german_qc3_file_prefix,
                                               keep = training_fam[, .(FID, IID)],
                                               extract = post_qc_overlapping_snps_file,
                                               ld.window.size = pruning_windowsize,
                                               ld.step.size = pruning_stepsize,
                                               ld.threshold = pruning_r2,
                                               exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s.chr%02d.prune.in",
                                      german_qc3_file_prefix, 1:22),
                              checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = pruning_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = pruning_reg))) {
  ids <- batchtools::findNotDone(reg = pruning_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = pruning_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_pruning"))
}

if (!batchtools::waitForJobs(reg = pruning_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(pruning_reg$file.dir)))
}

pruning_results <- batchtools::reduceResultsList(reg = pruning_reg)
saveRDS(pruning_results, pruning_results_rds)

writeLines(unlist(sapply(pruning_results, function(x) x$PRUNEIN, simplify = TRUE)),
           prune_in_snps_file)
