
source("init.R", chdir = TRUE)

## Load training fam ----
training_fam <- readRDS(training_fam_rds_file)

## GWAS on quality controlled data ----
gwa_reg <- load_or_create_registry("GWAS", conf.file = slurm_btconf_file)

if (!nrow(ids <- batchtools::findJobs(reg = gwa_reg))) {
  ids <- batchtools::batchMap(reg = gwa_reg,
                              fun = plink_gwa,
                              chr = 1:22,
                              more.args = list(out.dir = gwa_dir,
                                               plink.file.prefix = german_qc3_file_prefix,
                                               keep = training_fam[, .(FID, IID)],
                                               extract = post_qc_overlapping_snps_file,
                                               exec = plink_exec))
}

missing_files <- which(sapply(file.path(gwa_dir, sprintf("%s.chr%02d.assoc.logistic",
                                      basename(german_qc3_file_prefix), 1:22)),
                              checkmate::testFile) == FALSE)

if (length(missing_files)) {
  batchtools::resetJobs(reg = gwa_reg, missing_files)
}

if (nrow(batchtools::findNotDone(reg = gwa_reg))) {
  ids <- batchtools::findNotDone(reg = gwa_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = gwa_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 10,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_gwa"))
}

if (!batchtools::waitForJobs(reg = gwa_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(gwa_reg$file.dir)))
}

gwas_results <- rbindlist(
  batchtools::reduceResultsList(reg = gwa_reg,
                                fun = function(res, job) {

                                  gwa <- res$GWA[TEST != "SEX"]
                                  gwa[, DATASET := basename(dirname(job$pars$plink.file.prefix))]

                                  return(gwa)
                                })
  )

# Merge with p values from MB-MDR analysis
german_regional_mbmdr_reg <- batchtools::loadRegistry(
  file.dir = "/imbs/projects/grifam/data/proc/registries/german_regional_mbmdr_analyses",
  make.default = FALSE
)
german_regional_mbmdr_results <- rbindlist(
  lapply(
    reduceResultsList(reg = german_regional_mbmdr_reg),
    function(x) x$regional_mbmdr
  )
)
german_regional_mbmdr_results <- rbind(
  german_regional_mbmdr_results[, .(SNP = Marker1, pValue)],
  german_regional_mbmdr_results[, .(SNP = Marker2, pValue)]
)
setorder(german_regional_mbmdr_results, pValue)
gwas_results[german_regional_mbmdr_results[, .SD[1], by = SNP],
             P_MBMDR := i.pValue, on = "SNP"]
gwas_results[, P_GWAS := P]
gwas_results[, P := pmin(P_MBMDR, P_GWAS, na.rm = TRUE)]

saveRDS(gwas_results, gwas_results_rds)

# Write GWA result files
gwas_results[, data.table::fwrite(.SD, get(sprintf("%s_gwa_results_file",
                                                   tolower(unique(DATASET)))),
                                  sep = "\t"),
             by = .(DATASET),
             .SDcols = c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "STAT", "P")]
