
source("../init.R", chdir = TRUE)

##============================================================================+
## load R packages ----                                                       |
##============================================================================+
if(!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(data.table)
pacman::p_load(mlr)
pacman::p_load(snpStats)
pacman::p_load(plyr)
pacman::p_load(parallelMap)
pacman::p_load(ggplot2)

##============================================================================+
## mlr setup ----
##============================================================================+
mlr::configureMlr(on.learner.error = "warn",
                  on.par.without.desc = "warn")

##============================================================================+
## seed ----
##============================================================================+
seed <- 20170106
set.seed(seed)

##============================================================================+
## files & dirs ----
##============================================================================+
overlapping_snps_file <- file.path(proc_dir, "overlapping.snps")
post_qc_overlapping_snps_file <- file.path(proc_dir, "post_qc_overlapping.snps")
prune_in_snps_file <- file.path(proc_dir, "prune_in.snps")
clumped_snps_file <- file.path(proc_dir, "clumped.snps")
gwas_importance_rds_file <- file.path(proc_dir, "gwas_importance.rds")

# Define German datasets
german_datasets <- c("G1", "G2", "G3", "G4", "G5", "LURIC")

# Dataset proc dirs
dir.create(g1_dir <- file.path(proc_dir, "G1"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(g2_dir <- file.path(proc_dir, "G2"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(g3_dir <- file.path(proc_dir, "G3"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(g4_dir <- file.path(proc_dir, "G4"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(g5_dir <- file.path(proc_dir, "G5"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(luric_dir <- file.path(proc_dir, "LURIC"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(german_proc_dir <- file.path(proc_dir, "GERMAN"),
           showWarnings = FALSE, recursive = TRUE)
dir.create(wtccc_proc_dir <- file.path(proc_dir, "WTCCC"),
           showWarnings = FALSE, recursive = TRUE)
dir.create(cg_proc_dir <- file.path(proc_dir, "CG"),
           showWarnings = FALSE, recursive = TRUE)

# German
german_file_prefix <- file.path(german_proc_dir, "german")
german_merge_list_file <- file.path(german_proc_dir, "german.merge")

# German QC
german_qc1_file_prefix <- file.path(german_proc_dir, "german_QC1")
german_qc2_file_prefix <- file.path(german_proc_dir, "german_QC2")
german_qc3_file_prefix <- file.path(german_proc_dir, "german_QC3")
german_ibd_prune_in_file <- sprintf("%s.prune.in", german_qc2_file_prefix)

# German GWA
german_gwa_results_file <- file.path(gwa_dir, "german.assoc.logistic")

# German LD clumping
german_clump_in_file <- sprintf("%s.clumped", german_qc3_file_prefix)

# German dosage data
german_clumped_dosage_rds_file <- sprintf("%s_clumped_dosage.rds", german_qc3_file_prefix)
german_pruned_dosage_rds_file <- sprintf("%s_pruned_dosage.rds", german_qc3_file_prefix)

# German GWAA data
german_clumped_gwaa_rds_file <- sprintf("%s_clumped_gwaa.rds", german_qc3_file_prefix)
german_pruned_gwaa_rds_file <- sprintf("%s_pruned_gwaa.rds", german_qc3_file_prefix)

# WTCCC
wtccc_file_prefix <- file.path(wtccc_proc_dir, "wtccc")
wtccc_merge_list_file <- file.path(wtccc_proc_dir, "wtccc.merge")

# WTCCC QC
wtccc_qc1_file_prefix <- file.path(wtccc_proc_dir, "wtccc_QC1")
wtccc_qc2_file_prefix <- file.path(wtccc_proc_dir, "wtccc_QC2")
wtccc_qc3_file_prefix <- file.path(wtccc_proc_dir, "wtccc_QC3")
wtccc_ibd_prune_in_file <- sprintf("%s.prune.in", wtccc_qc2_file_prefix)

# WTCCC GWA
wtccc_gwa_results_file <- file.path(gwa_dir, "wtccc.assoc.logistic")

# German dosage data
wtccc_clumped_dosage_rds_file <- sprintf("%s_clumped_dosage.rds", wtccc_qc3_file_prefix)
wtccc_pruned_dosage_rds_file <- sprintf("%s_pruned_dosage.rds", wtccc_qc3_file_prefix)

# CG
cg_file_prefix <- file.path(cg_proc_dir, "cg")
cg_merge_list_file <- file.path(cg_proc_dir, "cg.merge")

# cg QC
cg_qc1_file_prefix <- file.path(cg_proc_dir, "cg_QC1")
cg_qc2_file_prefix <- file.path(cg_proc_dir, "cg_QC2")
cg_qc3_file_prefix <- file.path(cg_proc_dir, "cg_QC3")
cg_ibd_prune_in_file <- sprintf("%s.prune.in", cg_qc2_file_prefix)

# CG GWA
cg_gwa_results_file <- file.path(gwa_dir, "cg.assoc.logistic")

# German dosage data
cg_clumped_dosage_rds_file <- sprintf("%s_clumped_dosage.rds", cg_qc3_file_prefix)
cg_pruned_dosage_rds_file <- sprintf("%s_pruned_dosage.rds", cg_qc3_file_prefix)

# ID subsets
training_ids_rds_file <- file.path(proc_dir, "training_ids.rds")
validation_ids_rds_file <- file.path(proc_dir, "validation_ids.rds")

# Subset fam files
training_fam_rds_file <- file.path(proc_dir, "training_fam.rds")
validation_fam_rds_file <- file.path(proc_dir, "validation_fam.rds")

# ranger importance file
importance_rds_file <- file.path(proc_dir, "importance.rds")
pruned_importance_rds_file <- file.path(proc_dir, "pruned_importance.rds")
clumped_importance_rds_file <- file.path(proc_dir, "clumped_importance.rds")

# Tasks
training_task_rds_file <- file.path(proc_dir, "training_task.rds")
validation_task_rds_file <- file.path(proc_dir, "validation_task.rds")
wtccc_task_rds_file <- file.path(proc_dir, "wtccc_task.rds")
cg_task_rds_file <- file.path(proc_dir, "cg_task.rds")
task_G6_rds_file <- file.path(proc_dir, "task_G6.rds")
task_G7_rds_file <- file.path(proc_dir, "task_G7.rds")

# Benchmark
bmr_rds_file <- file.path(proc_dir, "bmr.rds")
