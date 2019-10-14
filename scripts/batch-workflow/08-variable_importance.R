
source("init.R", chdir = TRUE)

# Create dosage data ----
importance_reg <- load_or_create_registry(
  "importance",
  conf.file = slurm_btconf_file,
  packages = c("data.table", "GenABEL", "ranger")
)

cpus <- 16

if (!nrow(ids <- batchtools::findJobs(reg = importance_reg))) {
  ids <- batchtools::batchMap(reg = importance_reg,
                              fun = calc_ranger_importance,
                              gwaa.rds.file = c(german_clumped_gwaa_rds_file,
                                                german_pruned_gwaa_rds_file),
                              more.args = list(importance = "impurity_corrected",
                                               mtry.perc = 0.1,
                                               num.trees = 50000,
                                               num.threads = cpus))
}

if (nrow(batchtools::findNotDone(reg = importance_reg))) {
  ids <- batchtools::findNotDone(reg = importance_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = importance_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = cpus,
                                          memory = "100G",
                                          partition = "batch", walltime = 0,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_importance"))
}

if (!batchtools::waitForJobs(reg = importance_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(importance_reg$file.dir)))
}

importances <- reduceResultsList(
  fun = function(res, job) {
    importance_pvalues(res,
                       method = "janitza")
  },
  reg = importance_reg
)

importances <- rbindlist(
  lapply(
    importances, as.data.table, keep.rownames = "SNP"),
  fill = TRUE, idcol = "job.id"
)

saveRDS(importances, importance_rds_file)

clumped_importance <- importances[job.id == 1]
pruned_importance <- importances[job.id == 2]
clumped_importance[, job.id := NULL]
pruned_importance[, job.id := NULL]

saveRDS(clumped_importance, clumped_importance_rds_file)
saveRDS(pruned_importance, pruned_importance_rds_file)
