
source("init.R", chdir = TRUE)

# Load ag_cardio/classification training task ----
training_task <- readRDS(training_task_rds_file)

# Extract feature names ----
training_features <- mlr::getTaskFeatureNames(training_task)

writeLines(
  text = gsub("^(X)(\\d+)\\.(\\d+.*)", "\\2:\\3", training_features),
  con = file.path(proc_dir, "training_features.list")
)

# Get positions of training features ----
system2(
  command = "grep",
  args = c("-wF",
           "-f", file.path(proc_dir, "training_features.list"),
           "/imbs/projects/ag_cardio/proc/classification/GERMAN/german.bim",
           "|", "cut", "-f1,4",
           ">", file.path(proc_dir, "training_features.positions"))
)

# Create list of feature names with positions ----
system2(
  command = "grep",
  args = c("-wF",
           "-f", file.path(proc_dir, "training_features.list"),
           "/imbs/projects/ag_cardio/proc/classification/GERMAN/german.bim",
           "|", "cut", "-f2,1,4,5,6",
           ">", file.path(proc_dir, "training_features_positions.list"))
)

# Extract dosage data of training features ----
vcf_files <- c(
  list.files(
    path = file.path(input_dir, "G6"),
    pattern = ".*vcf.gz$",
    recursive = TRUE,
    full.names = TRUE
  ),
  list.files(
    path = file.path(input_dir, "G7"),
    pattern = ".*vcf.gz$",
    recursive = TRUE,
    full.names = TRUE
  )
)

sapply(
  unique(
    dirname(
      dosage_files <- file.path(proc_dir, gsub(".*\\/(G[67])\\/.*", "\\1", vcf_files), sprintf("%s_%s.dosage", gsub(".*\\/(G[67])\\/.*", "\\1", vcf_files), gsub("^(\\d+).*", "\\1", basename(vcf_files))))
    )
  ), dir.create, recursive = TRUE, showWarnings = FALSE
)

sapply(
  unique(
    dirname(
      index_files_to <- file.path(proc_dir, gsub(".*\\/(G[67])\\/.*", "\\1", vcf_files), sprintf("%s.tbi", basename(vcf_files)))
    )
  ), dir.create, recursive = TRUE, showWarnings = FALSE
)

reg_name <- "cardio_classification-create_index"

reg <- imbs::load_or_create_registry(
  file.dir = file.path(reg_dir, reg_name),
  work.dir = proc_dir,
  writeable = TRUE,
  overwrite = FALSE,
  conf.file = slurm_btconf_file
)

ids <- batchtools::batchMap(
  fun = bcftools_index,
  vcf.file = vcf_files,
  output.file = index_files,
  more.args = list(bcftools.exe = bcftools_exec)
)

ids[, chunk := 1]

batchtools::submitJobs(
  ids = ids,
  resources = list(ntasks = 1, ncpus = 1, memory = "5G",
                   walltime = 59,
                   partition = "fast",
                   name = reg_name,
                   chunks.as.arrayjobs = TRUE)
)

batchtools::waitForJobs()

# Indices copied to external_data manually on head

reg_name <- "cardio_classification-extract_features_dosage"

reg <- imbs::load_or_create_registry(
  file.dir = file.path(reg_dir, reg_name),
  work.dir = proc_dir,
  writeable = TRUE,
  overwrite = FALSE,
  conf.file = slurm_btconf_file
)

ids <- batchtools::batchMap(
  fun = bcftools_query,
  vcf.file = vcf_files,
  output.file = dosage_files,
  more.args = list(query = "%ID\\\\t%CHROM\\\\t%POS\\\\t%REF\\\\t%ALT[\\\\t%DS]\\\\n",
                   regions.file = file.path(proc_dir, "training_features.positions"),
                   print.header = TRUE,
                   bcftools.exe = bcftools_exec)
)

ids[, chunk := 1]

setJobNames(ids, names = basename(dosage_files))

batchtools::submitJobs(
  ids = ids,
  resources = list(ntasks = 1, ncpus = 1, memory = "5G",
                   partition = "batch",
                   name = reg_name,
                   chunks.as.arrayjobs = TRUE)
)

batchtools::waitForJobs()

dosage_data <- reduceResultsList()
names(dosage_data) <- basename(dosage_files)

dosage_data_G6 <- rbindlist(dosage_data[which(grepl("G6", names(dosage_data)))])
dosage_data_G7 <- rbindlist(dosage_data[which(grepl("G7", names(dosage_data)))])

dosage_data <- list(G6 = dosage_data_G6, G7 = dosage_data_G7)

training_features_positions <- data.table::fread(
  input = file.path(proc_dir, "training_features_positions.list"),
  col.names = c("CHR", "RSID", "POS", "REF", "ALT")
)

dosage_data_molten <- list()

# Add information from training and flip mismatched alleles ----
for (i in 1:length(dosage_data)) {
  dosage_data[[i]][training_features_positions,
                   c("RSID", "REF", "ALT") := list(i.RSID, i.REF, i.ALT),
                   on = c("[2]CHROM" = "CHR", "[3]POS" = "POS")]
  dosage_data_molten[[i]] <- melt(
    data = dosage_data[[i]][(`[4]REF` == REF & `[5]ALT` == ALT) |
                              (`[5]ALT` == REF & `[4]REF` == ALT)],
    id.vars = c("# [1]ID", "[2]CHROM", "[3]POS", "[4]REF", "[5]ALT", "REF", "ALT", "RSID"),
    variable.name = "SAMPLE",
    value.name = "DOSAGE"
  )
  dosage_data_molten[[i]][`[4]REF` != REF, DOSAGE := 2 - DOSAGE]
}
names(dosage_data_molten) <- names(dosage_data)

# Merge datasets ----
dosage_data_molten <- rbindlist(dosage_data_molten, idcol = "DATASET")

# Transpose dosage data ----
dosage <- dcast(
  data = dosage_data_molten,
  formula = DATASET + SAMPLE ~ RSID,
  value.var = "DOSAGE"
)

# Define status ----
dosage[, SAMPLE := gsub("\\[\\d+\\](.*):DS", "\\1", SAMPLE)]

samples_G6 <- data.table::fread(file.path(input_dir, "G6/HRC_imputation_G6/sample_file/vcf2gen_SNPpostHRCimputeQC_chr1.sample"), skip = 1, col.names = c("ID1", "ID2", "missing", "sex", sprintf("PC%02d", 1:10), "STATUS"))
samples_G7 <- data.table::fread(file.path(input_dir, "G7/HRC_imputation_G7/sample_file/g7.sample+new"), skip = 2, col.names = c("ID1", "ID2", "missing", "sex", "STATUS", sprintf("PC%02d", 1:6)))

dosage[samples_G6, STATUS := i.STATUS, on = c("SAMPLE" = "ID2")]
dosage[samples_G7, STATUS := i.STATUS, on = c("SAMPLE" = "ID2")]
dosage[, STATUS := factor(STATUS, levels = c(0, 1), labels = c("control", "case"))]

# Fix names ----
data.table::setnames(dosage, make.names(names(dosage)))

# Define 1000GP3 based task ----
task_G6_data <- as.data.frame(
  dosage[grepl("G6", DATASET)]
)
rownames(task_G6_data) <- task_G6_data$SAMPLE

task_G6 <- mlr::makeClassifTask(
  id = "GerMIFSVI",
  data = task_G6_data[, !names(task_G6_data) %in% c("DATASET", "SAMPLE")],
  target = "STATUS",
  positive = "control"
)
saveRDS(task_G6, task_G6_rds_file)

# Define HRC based task ----
task_G7_data <- as.data.frame(
  dosage[grepl("G7", DATASET)]
)
rownames(task_G7_data) <- task_G7_data$SAMPLE

task_G7 <- mlr::makeClassifTask(
  id = "GerMIFSVII",
  data = task_G7_data[, !names(task_G7_data) %in% c("DATASET", "SAMPLE")],
  target = "STATUS",
  positive = "control"
)
saveRDS(task_G7, task_G7_rds_file)

# Define measures ----
measures <- list(mlr::auc, # use AUC as performance measure during tuning
                 mlr::mmce, mlr::npv, mlr::fpr, mlr::f1, mlr::fnr, mlr::ssr,
                 mlr::tp, mlr::tn, mlr::gpr, mlr::lsr, mlr::acc,
                 mlr::wkappa, mlr::ppv, mlr::logloss, mlr::ber, mlr::tpr,
                 mlr::brier, mlr::gmean, mlr::fdr, mlr::tnr, mlr::qsr, mlr::bac,
                 mlr::brier.scaled, mlr::fp, mlr::fn, mlr::kappa)

# Load models ----
rs_model <- readRDS(rs_model_rds_file)
nb_model <- readRDS(nb_model_rds_file)
rf_model <- readRDS(rf_model_rds_file)
gb_model <- readRDS(gb_model_rds_file)
svm_model <- readRDS(svm_model_rds_file)

# Dirty hotfix to adjust for min and max score
rs_model$learner.model$next.model$learner.model$next.model$learner.model$min_score <- rs_model$learner.model$next.model$learner.model$next.model$learner.model$risk_model[weight < 0, sum(2*weight)]
rs_model$learner.model$next.model$learner.model$next.model$learner.model$max_score <- rs_model$learner.model$next.model$learner.model$next.model$learner.model$risk_model[weight > 0, sum(2*weight)]

# Predict ----
# GerMIFSVI
rs_g6_prediction <- predict(rs_model, task = task_G6)
saveRDS(rs_g6_prediction, rs_g6_prediction_rds_file)

nb_g6_prediction <- predict(nb_model, task = task_G6)
saveRDS(nb_g6_prediction, nb_g6_prediction_rds_file)

rf_g6_prediction <- predict(rf_model, task = task_G6)
saveRDS(rf_g6_prediction, rf_g6_prediction_rds_file)

gb_g6_prediction <- predict(gb_model, task = task_G6)
saveRDS(gb_g6_prediction, gb_g6_prediction_rds_file)

svm_g6_prediction <- predict(svm_model, task = task_G6)
saveRDS(svm_g6_prediction, svm_g6_prediction_rds_file)

# GerMIFSVII
rs_g7_prediction <- predict(rs_model, task = task_G7)
saveRDS(rs_g7_prediction, rs_g7_prediction_rds_file)

nb_g7_prediction <- predict(nb_model, task = task_G7)
saveRDS(nb_g7_prediction, nb_g7_prediction_rds_file)

rf_g7_prediction <- predict(rf_model, task = task_G7)
saveRDS(rf_g7_prediction, rf_g7_prediction_rds_file)

gb_g7_prediction <- predict(gb_model, task = task_G7)
saveRDS(gb_g7_prediction, gb_g7_prediction_rds_file)

svm_g7_prediction <- predict(svm_model, task = task_G7)
saveRDS(svm_g7_prediction, svm_g7_prediction_rds_file)

# Performance ----
# g6
rs_g6_performance <- performance(rs_g6_prediction,
                                 measures = measures)
nb_g6_performance <- performance(nb_g6_prediction,
                                 measures = measures)
rf_g6_performance <- performance(rf_g6_prediction,
                                 measures = measures)
gb_g6_performance <- performance(gb_g6_prediction,
                                 measures = measures)
svm_g6_performance <- performance(svm_g6_prediction,
                                  measures = measures)
# G7
rs_g7_performance <- performance(rs_g7_prediction,
                                 measures = measures)
nb_g7_performance <- performance(nb_g7_prediction,
                                 measures = measures)
rf_g7_performance <- performance(rf_g7_prediction,
                                 measures = measures)
gb_g7_performance <- performance(gb_g7_prediction,
                                 measures = measures)
svm_g7_performance <- performance(svm_g7_prediction,
                                  measures = measures)
