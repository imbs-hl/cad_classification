
# Consumes a lot of memory!
# Use interactive session with --mem=100G
source("init.R", chdir = TRUE)

# Load rds files from previous steps ----
plink_files <- readRDS(plink_files_rds)
gen_files <- readRDS(gen_files_rds)
sample_files <- readRDS(sample_files_rds)

# Create dosage data ----
dosage_reg <- load_or_create_registry("dosage",
                                      conf.file = slurm_btconf_file,
                                      packages = c("data.table"))

params <- gen_files[plink_files[sample_files, on = .(STUDY)], on = .(STUDY, CHR)]
params <- rbind(params, params)
params[, EXTRACT := rep(c(clumped_snps_file, prune_in_snps_file), each = .N/2)]
params[, REDUCTION := rep(c("clumping", "pruning"), each = .N/2)]
params[, PREFIX := sprintf("%s_%s", PLINK, REDUCTION)]
if (!nrow(ids <- batchtools::findJobs(reg = dosage_reg))) {
  ids <- batchtools::batchMap(reg = dosage_reg,
                              fun = gen2dosage,
                              chr = params$CHR,
                              gen.file = params$GZFILE,
                              sample.file = params$SAMPLEFILE,
                              plink.file.prefix = params$PLINK,
                              out.prefix = params$PREFIX,
                              extract = params$EXTRACT,
                              keep = sprintf("%s.fam", sapply(sprintf("%s_qc3_file_prefix", tolower(params$DATASET)), get)),
                              more.args = list(qctool.exec = qctool_exec,
                                               fcgene.exec = fcgene_exec))
}

missing_files <- which(sapply(sprintf("%s_dose.txt",
                                      params$PREFIX),
                              checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = dosage_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = dosage_reg))) {
  ids <- batchtools::findNotDone(reg = dosage_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = dosage_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_dosage"))
}

if (!batchtools::waitForJobs(reg = dosage_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(dosage_reg$file.dir)))
}

job_pars <- batchtools::getJobPars(reg = dosage_reg)
job_pars[, DATASET := basename(dirname(gen.file))]
job_pars[, REDUCTION := ifelse(grepl("clumping", out.prefix), "clumped", "pruned")]

# Merge dosages of datasets
# Consumes a lot of memory!
# Use interactive session with --mem=100G
clumped_dosage_data <- job_pars[REDUCTION == "clumped", batchtools::reduceResults(
  fun = function(aggr, res) {
    setkey(aggr, ID, SEX, STATUS)
    setkey(res, ID, SEX, STATUS)
    merge(aggr, res, all = TRUE, by = c("ID", "SEX", "STATUS"))
  },
  ids = job.id),
  by = c("DATASET")]
clumped_dosage_data[DATASET %in% german_datasets, DATASET := "GERMAN"]
clumped_dosage_data[, saveRDS(.SD, get(sprintf("%s_clumped_dosage_rds_file", tolower(unique(DATASET))))),
                    by = c("DATASET")]
rm(clumped_dosage_data)
pruned_dosage_data <- job_pars[REDUCTION == "pruned", batchtools::reduceResults(
  fun = function(aggr, res) {
    setkey(aggr, ID, SEX, STATUS)
    setkey(res, ID, SEX, STATUS)
    merge(aggr, res, all = TRUE, by = c("ID", "SEX", "STATUS"))
  },
  ids = job.id),
  by = c("DATASET")]
pruned_dosage_data[DATASET %in% german_datasets, DATASET := "GERMAN"]
pruned_dosage_data[, saveRDS(.SD, get(sprintf("%s_pruned_dosage_rds_file", tolower(unique(DATASET))))),
                    by = c("DATASET")]
rm(pruned_dosage_data)
