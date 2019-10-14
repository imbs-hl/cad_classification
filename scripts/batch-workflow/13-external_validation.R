
source("init.R", chdir = TRUE)

# Load task ----
wtccc_task <- readRDS(wtccc_task_rds_file)
cg_task <- readRDS(cg_task_rds_file)

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

# Predict ----
# WTCCC
rs_wtccc_prediction <- predict(rs_model, task = wtccc_task)
nb_wtccc_prediction <- predict(nb_model, task = wtccc_task)
rf_wtccc_prediction <- predict(rf_model, task = wtccc_task)
gb_wtccc_prediction <- predict(gb_model, task = wtccc_task)
svm_wtccc_prediction <- predict(svm_model, task = wtccc_task)

saveRDS(rs_wtccc_prediction, rs_wtccc_prediction_rds_file)
saveRDS(nb_wtccc_prediction, nb_wtccc_prediction_rds_file)
saveRDS(rf_wtccc_prediction, rf_wtccc_prediction_rds_file)
saveRDS(gb_wtccc_prediction, gb_wtccc_prediction_rds_file)
saveRDS(svm_wtccc_prediction, svm_wtccc_prediction_rds_file)

# rs_wtccc_prediction <- readRDS(rs_wtccc_prediction_rds_file)
# nb_wtccc_prediction <- readRDS(nb_wtccc_prediction_rds_file)
# rf_wtccc_prediction <- readRDS(rf_wtccc_prediction_rds_file)
# gb_wtccc_prediction <- readRDS(gb_wtccc_prediction_rds_file)
# svm_wtccc_prediction <- readRDS(svm_wtccc_prediction_rds_file)

# CG
rs_cg_prediction <- predict(rs_model, task = cg_task)
nb_cg_prediction <- predict(nb_model, task = cg_task)
rf_cg_prediction <- predict(rf_model, task = cg_task)
gb_cg_prediction <- predict(gb_model, task = cg_task)
svm_cg_prediction <- predict(svm_model, task = cg_task)

saveRDS(rs_cg_prediction, rs_cg_prediction_rds_file)
saveRDS(nb_cg_prediction, nb_cg_prediction_rds_file)
saveRDS(rf_cg_prediction, rf_cg_prediction_rds_file)
saveRDS(gb_cg_prediction, gb_cg_prediction_rds_file)
saveRDS(svm_cg_prediction, svm_cg_prediction_rds_file)

# rs_cg_prediction <- readRDS(rs_cg_prediction_rds_file)
# nb_cg_prediction <- readRDS(nb_cg_prediction_rds_file)
# rf_cg_prediction <- readRDS(rf_cg_prediction_rds_file)
# gb_cg_prediction <- readRDS(gb_cg_prediction_rds_file)
# svm_cg_prediction <- readRDS(svm_cg_prediction_rds_file)

# Performance ----
# WTCCC
rs_wtccc_performance <- performance(rs_wtccc_prediction,
                                    measures = measures)
nb_wtccc_performance <- performance(nb_wtccc_prediction,
                                    measures = measures)
rf_internal_wtccc_performance <- performance(rf_wtccc_prediction,
                                                  measures = measures)
gb_internal_wtccc_performance <- performance(gb_wtccc_prediction,
                                                  measures = measures)
svm_internal_wtccc_performance <- performance(svm_wtccc_prediction,
                                                   measures = measures)
# CG
rs_cg_performance <- performance(rs_cg_prediction,
                                 measures = measures)
nb_cg_performance <- performance(nb_cg_prediction,
                                 measures = measures)
rf_internal_cg_performance <- performance(rf_cg_prediction,
                                                  measures = measures)
gb_internal_cg_performance <- performance(gb_cg_prediction,
                                                  measures = measures)
svm_internal_cg_performance <- performance(svm_cg_prediction,
                                                   measures = measures)
