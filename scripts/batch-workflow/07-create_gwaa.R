
source("init.R", chdir = TRUE)

training_fam <- readRDS(training_fam_rds_file)

# Create GenABEL data ----
gwaa_reg <- load_or_create_registry("gwaa",
                                    conf.file = slurm_btconf_file,
                                    packages = c("data.table", "GenABEL"))

if (!nrow(ids <- batchtools::findJobs(reg = gwaa_reg))) {
  ids <- batchtools::batchMap(reg = gwaa_reg,
                              fun = plink2gwaa,
                              extract = c(clumped_snps_file, prune_in_snps_file),
                              out.prefix = sprintf(c("%s_clumped", "%s_pruned"), german_qc3_file_prefix),
                              more.args = list(plink.exec = plink_exec,
                                               plink.file.prefix = german_qc3_file_prefix,
                                               keep = training_fam[, .(FID, IID)],
                                               fcgene.exec = fcgene_exec))
}

if (nrow(batchtools::findNotDone(reg = gwaa_reg))) {
  ids <- batchtools::findNotDone(reg = gwaa_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = gwaa_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 10000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_gwaa"))
}

if (!batchtools::waitForJobs(reg = gwaa_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(gwaa_reg$file.dir)))
}

job_pars <- batchtools::getJobPars(reg = gwaa_reg)
job_pars[, REDUCTION := ifelse(grepl("clumped", out.prefix), "clumped", "pruned")]

results <- reduceResultsList(reg = gwaa_reg, ids = job_pars$job.id)

saveRDS(results[[1]], get(sprintf("german_%s_gwaa_rds_file", job_pars[1, REDUCTION])))
saveRDS(results[[2]], get(sprintf("german_%s_gwaa_rds_file", job_pars[2, REDUCTION])))
