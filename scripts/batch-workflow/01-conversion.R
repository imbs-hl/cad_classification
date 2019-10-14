
source("init.R", chdir = TRUE)

# Read sample files ----
sample_files <- data.table(SAMPLEFILE = list.files(input_dir,
                                                   recursive = TRUE,
                                                   pattern = "SampleInfos.txt",
                                                   full.names = TRUE))
sample_files[, STUDY := basename(dirname(SAMPLEFILE))]

samples <- rbindlist(sapply(sample_files$SAMPLEFILE,
                            function(x) {
                              s <- data.table::fread(x, skip = 2)[, study := basename(dirname(x))]
                              setnames(s, grep("V.*", colnames(s), value = TRUE), unlist(strsplit(readLines(x, 1), " ")))
                            },
                            simplify = FALSE))

saveRDS(samples, samples_rds)

# Convert compression to gzip ----
gen_files <- data.table::data.table(BZFILE = list.files(input_dir,
                                                        recursive = TRUE,
                                                        pattern = "*.bz2",
                                                        full.names = TRUE))
gen_files[, STUDY := basename(dirname(BZFILE))]
gen_files[, CHR := as.integer(gsub(".*_(\\d+)\\..*", "\\1", BZFILE))]

gen_files[, GZFILE := file.path(proc_dir,
                                basename(dirname(BZFILE)),
                                gsub(pattern = "(.*)bz2",
                                     replacement = "\\1gz",
                                     x = basename(BZFILE)))]
gen_files[, DATASET := STUDY]
gen_files[STUDY %in% german_datasets, DATASET := "GERMAN"]

bz2gz_reg <- load_or_create_registry("bz2gz", conf.file = slurm_btconf_file)

if(!nrow(ids <- batchtools::findJobs(reg = bz2gz_reg))) {
  ids <- batchtools::batchMap(reg = bz2gz_reg,
                              fun = bz2gz,
                              bz.file = gen_files$BZFILE,
                              gz.file = gen_files$GZFILE)
}

missing_files <- gen_files[, which(!file.exists(GZFILE))]

if(length(missing_files)) {
  batchtools::resetJobs(reg = bz2gz_reg, ids = missing_files)
}

if(nrow(batchtools::findNotDone(reg = bz2gz_reg))) {
  ids <- batchtools::findNotDone(reg = bz2gz_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = bz2gz_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_bz2gz"))
}

if(!batchtools::waitForJobs(reg = bz2gz_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(bz2gz_reg$file.dir)))
}

saveRDS(gen_files, gen_files_rds)

# Copy the sample files to new location and correct the phenotype coding ----
sample_files[, ORIGSAMPLEFILE := SAMPLEFILE]
sample_files[, SAMPLEFILE := sapply(X = SAMPLEFILE,
                                    FUN = copy_sample_files)]

# Convert GEN to PLINK ----
sample_files[, PHENOCOLNAME := sapply(X = SAMPLEFILE,
                                      FUN = function(x) {
                                        unlist(strsplit(readLines(x, n = 1), " "))[6]
                                      })]

saveRDS(sample_files, sample_files_rds)

gz2bed_reg <- load_or_create_registry("gz2bed", conf.file = slurm_btconf_file)

if(!nrow(ids <- batchtools::findJobs(reg = gz2bed_reg))) {
  params <- gen_files[sample_files, on = "STUDY"]
  ids <- batchtools::batchMap(reg = gz2bed_reg,
                              fun = gz2bed,
                              gz.file = params$GZFILE,
                              sample.file = params$SAMPLEFILE,
                              pheno = params$PHENOCOLNAME,
                              more.args = list(exec = plink_exec,
                                               call.threshold = call_threshold))
}

missing_files <- which((sapply(paste(sub("^([^.]*).*", "\\1", gen_files$GZFILE), "bed", sep = "."),
                               checkmate::testFile) &
                          sapply(paste(sub("^([^.]*).*", "\\1", gen_files$GZFILE), "bim", sep = "."),
                                 checkmate::testFile) &
                          sapply(paste(sub("^([^.]*).*", "\\1", gen_files$GZFILE), "fam", sep = "."),
                                 checkmate::testFile)) == FALSE)

if(length(missing_files)) {
  batchtools::resetJobs(reg = gz2bed_reg, missing_files)
}

if(nrow(batchtools::findNotDone(reg = gz2bed_reg))) {
  ids <- batchtools::findNotDone(reg = gz2bed_reg, ids = ids)
  ids[, chunk := 1]
  batchtools::submitJobs(reg = gz2bed_reg,
                         ids = ids,
                         resources = list(ntasks = 1, ncpus = 1, memory = 6000,
                                          partition = "fast", walltime = 60,
                                          chunks.as.arrayjobs = TRUE,
                                          name = "cad_classification_gz2bed"))
}

if(!batchtools::waitForJobs(reg = gz2bed_reg)) {
  stop(sprintf("Jobs of registry '%s' failed!", basename(gz2bed_reg$file.dir)))
}
