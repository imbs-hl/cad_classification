
source("init.R", chdir = TRUE)

# Load rds files from previous steps ----
sample_files <- readRDS(sample_files_rds)
gen_files <- readRDS(gen_files_rds)

# Write merge list ----
gz2bed_reg <- load_or_create_registry("gz2bed", conf.file = slurm_btconf_file)

plink_files <- batchtools::reduceResultsDataTable(reg = gz2bed_reg, fun = function(job, res) list(PLINK = res$out))
plink_files[, STUDY := basename(dirname(PLINK))]
plink_files[, CHR := as.integer(gsub(".*_(\\d+)", "\\1", PLINK))]
plink_files[, DATASET := STUDY]
plink_files[STUDY %in% german_datasets, DATASET := "GERMAN"]

plink_files[, writeLines(PLINK, get(sprintf("%s_merge_list_file",
                                            tolower(unique(DATASET))))),
            by = DATASET]

# Use PLINK to merge datasets
merge_reg <- load_or_create_registry("mergeDatasets", conf.file = slurm_btconf_file)

if(!nrow(ids <- batchtools::findJobs(reg = merge_reg))) {
  ids <- batchtools::batchMap(reg = merge_reg,
                             fun = merge_plink,
                             mergelist = c(german_merge_list_file, cg_merge_list_file, wtccc_merge_list_file),
                             output = c(german_file_prefix, cg_file_prefix, wtccc_file_prefix),
                             more.args = list(exec = plink_exec))
}

missing_files <- which(sapply(sprintf("%s.bed", c(german_file_prefix, cg_file_prefix, wtccc_file_prefix)),
                              checkmate::testFile) == FALSE)

if(length(missing_files)) {
  batchtools::resetJobs(reg = merge_reg, missing_files)
}

if(nrow(batchtools::findNotDone(reg = merge_reg))) {
  ids <- batchtools::findNotDone(reg = merge_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = merge_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = "90G",
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_merge"))
}

if(!batchtools::waitForJobs(reg = merge_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(merge_reg$file.dir)))
}

saveRDS(plink_files, plink_files_rds)
