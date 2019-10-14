
source("init.R", chdir = TRUE)

# German data ----
# Load data ----
german_dosagedata <- copy(readRDS(german_pruned_dosage_rds_file))
gwas_results <- readRDS(gwas_results_rds)
pruned_importance <- readRDS(pruned_importance_rds_file)

pruned_importance <- pruned_importance[importance > 0]

# Define GWAS importance ----
gwas_importance <- gwas_results[DATASET == "GERMAN"][pruned_importance, on = "SNP"][, setNames(log(OR), make.names(SNP))]

saveRDS(gwas_importance, gwas_importance_rds_file)

# Define task ----
german_dosagedata[, `:=`(ID = NULL, SEX = NULL)]
data.table::setnames(german_dosagedata, make.names(names(german_dosagedata)))
task <- mlr::makeClassifTask(id = "CAD", data = german_dosagedata, target = "STATUS")

validation_ids <- readRDS(validation_ids_rds_file)
training_ids <- readRDS(training_ids_rds_file)

validation_task <- mlr::subsetTask(task,
                                   subset = validation_ids,
                                   features = pruned_importance[make.names(SNP)])
training_task <- mlr::subsetTask(task,
                                 subset = training_ids,
                                 features = pruned_importance[make.names(SNP)])

saveRDS(validation_task, validation_task_rds_file)
saveRDS(training_task, training_task_rds_file)

rm(german_dosagedata)
rm(gwas_results)
rm(pruned_importance)
rm(task)
rm(training_task)
rm(validation_task)
gc()

# WTCCC data ----
# Load data ----
wtccc_dosagedata <- copy(readRDS(wtccc_pruned_dosage_rds_file))

# Define task ----
wtccc_dosagedata[, `:=`(ID = NULL, SEX = NULL)]
data.table::setnames(wtccc_dosagedata, make.names(names(wtccc_dosagedata)))

wtccc_task <- mlr::makeClassifTask(id = "CAD_WTCCC",
                                   data = wtccc_dosagedata,
                                   target = "STATUS")

saveRDS(wtccc_task, wtccc_task_rds_file)

rm(wtccc_dosagedata)
rm(wtccc_task)
gc()

# CG data ----
# Load data ----
cg_dosagedata <- copy(readRDS(cg_pruned_dosage_rds_file))

# Define task ----
cg_dosagedata[, `:=`(ID = NULL, SEX = NULL)]
data.table::setnames(cg_dosagedata, make.names(names(cg_dosagedata)))

cg_task <- mlr::makeClassifTask(id = "CAD_WTCCC",
                                   data = cg_dosagedata,
                                   target = "STATUS")

saveRDS(cg_task, cg_task_rds_file)

rm(cg_dosagedata)
rm(cg_task)
gc()
