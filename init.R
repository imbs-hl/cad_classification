
##============================================================================+
## directories setup ----
##============================================================================+
user <- Sys.info()["user"]
switch (user,
        gola = main_dir <- path.expand("~/Documents/Work/ag_cardio-classification/")
)

dir.create(scripts_dir <- file.path(main_dir, "scripts"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir <- file.path(main_dir, "report"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(paper_dir <- file.path(main_dir, "paper"),
           recursive = TRUE, showWarnings = FALSE)

input_dir <- "/imbs/external_data/ag_cardio/imputed_munich"
dir.create(proc_dir <- "/imbs/projects/ag_cardio/proc/classification",
           recursive = TRUE, showWarnings = FALSE)

dir.create(reg_dir <- file.path(proc_dir, "registries"),
           recursive = TRUE, showWarnings = FALSE)

# Analyses
dir.create(analyses_dir <- file.path(proc_dir, "analyses"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(gwa_dir <- file.path(analyses_dir, "gwa"),
           recursive = TRUE, showWarnings = FALSE)

##============================================================================+
## rds files ----
##============================================================================+
sample_files_rds <- file.path(proc_dir, "sample_files.rds")
samples_rds <- file.path(proc_dir, "samples.rds")
gen_files_rds <- file.path(proc_dir, "gen_files.rds")
plink_files_rds <- file.path(proc_dir, "plink_files.rds")
snp_stats_rds <- file.path(proc_dir, "snp_stats.rds")
gwas_results_rds <- file.path(gwa_dir, "gwas_results.rds")
pruning_results_rds <- file.path(proc_dir, "pruning.results.rds")
clumped_results_rds <- file.path(proc_dir, "clumped.results.rds")

##============================================================================+
## parameters ----
##============================================================================+
info_threshold <- 0.9 # Imputed SNPs with lower info score are filtered
call_threshold <- 0.1 # Sample genotype based on imputation probabilities with certainty
# QC
het_fac <- 3
hwe_pvalue <- 1e-5
geno <- 0.02
mind <- 0.02
maf <- 0.05
# IBD
ibd_ld_window_size <- 2e6 # 2 MB
ibd_ld_step_size <- 0.1*ibd_ld_window_size
ibd_ld_threshold <- 0.2
pi_hat <- 0.2
# Pruning
pruning_windowsize <- 10000
pruning_r2 <- 0.5
pruning_stepsize <- 0.1 * pruning_windowsize
# Clumping
clumping_p1 <- 0.25
clumping_p2 <- 1
clumping_kb <- 10000
clumping_r2 <- 0.5

##============================================================================+
## Training/Validation ----
##============================================================================+
validation_size <- 1000

##============================================================================+
## benchmark parameters ----
##============================================================================+
bm_stratify <- TRUE
bm_dw <- 1/5
bm_inner_cv_iters <- 5
bm_outer_cv_iters <- 10
bm_tuning_iters <- 100
bm_time_budget <- 4*7*24*60*60 # 4 weeks

##============================================================================+
## tuning parameters ----
##============================================================================+
tuning_stratify <- TRUE
tuning_cv_iters <- 10
tuning_iters <- 100
tuning_time_budget <- 4*7*24*60*60 # 4 weeks

##============================================================================+
## external software ----
##============================================================================+
plink_exec <- "plink1.9b4.4"
qctool_exec <- "qctool1.4"
fcgene_exec <- "fcgene"
bcftools_exec <- "bcftools1.9"

##============================================================================+
## tuned models ----
##============================================================================+
rs_model_rds_file <- file.path(proc_dir, "rs_model.rds") # riskScoreR model
nb_model_rds_file <- file.path(proc_dir, "nb_model.rds") # Naive Bayes model
rf_model_rds_file <- file.path(proc_dir, "rf_model.rds") # RandomForest model
svm_model_rds_file <- file.path(proc_dir, "svm_model.rds") # SVM model
gb_model_rds_file <- file.path(proc_dir, "gb_model.rds") # Gradient boosting model

##============================================================================+
## predictions ----
##============================================================================+
rs_validation_prediction_rds_file <- file.path(proc_dir, "rs_validation_prediction.rds")
nb_validation_prediction_rds_file <- file.path(proc_dir, "nb_validation_prediction.rds")
rf_validation_prediction_rds_file <- file.path(proc_dir, "rf_validation_prediction.rds")
gb_validation_prediction_rds_file <- file.path(proc_dir, "gb_validation_prediction.rds")
svm_validation_prediction_rds_file <- file.path(proc_dir, "svm_validation_prediction.rds")
rs_wtccc_prediction_rds_file <- file.path(proc_dir, "rs_wtccc_prediction.rds")
nb_wtccc_prediction_rds_file <- file.path(proc_dir, "nb_wtccc_prediction.rds")
rf_wtccc_prediction_rds_file <- file.path(proc_dir, "rf_wtccc_prediction.rds")
gb_wtccc_prediction_rds_file <- file.path(proc_dir, "gb_wtccc_prediction.rds")
svm_wtccc_prediction_rds_file <- file.path(proc_dir, "svm_wtccc_prediction.rds")
rs_cg_prediction_rds_file <- file.path(proc_dir, "rs_cg_prediction.rds")
nb_cg_prediction_rds_file <- file.path(proc_dir, "nb_cg_prediction.rds")
rf_cg_prediction_rds_file <- file.path(proc_dir, "rf_cg_prediction.rds")
gb_cg_prediction_rds_file <- file.path(proc_dir, "gb_cg_prediction.rds")
svm_cg_prediction_rds_file <- file.path(proc_dir, "svm_cg_prediction.rds")
rs_g6_prediction_rds_file <- file.path(proc_dir, "rs_g6_prediction")
nb_g6_prediction_rds_file <- file.path(proc_dir, "nb_g6_prediction")
rf_g6_prediction_rds_file <- file.path(proc_dir, "rf_g6_prediction")
gb_g6_prediction_rds_file <- file.path(proc_dir, "gb_g6_prediction")
svm_g6_prediction_rds_file <- file.path(proc_dir, "svm_g6_prediction")
rs_g7_prediction_rds_file <- file.path(proc_dir, "rs_g7_prediction")
nb_g7_prediction_rds_file <- file.path(proc_dir, "nb_g7_prediction")
rf_g7_prediction_rds_file <- file.path(proc_dir, "rf_g7_prediction")
gb_g7_prediction_rds_file <- file.path(proc_dir, "gb_g7_prediction")
svm_g7_prediction_rds_file <- file.path(proc_dir, "svm_g7_prediction")
