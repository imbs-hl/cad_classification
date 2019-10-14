
source("init.R", chdir = TRUE)

# Load task ----
validation_task <- readRDS(validation_task_rds_file)

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
rs_validation_prediction <- predict(rs_model, task = validation_task)
nb_validation_prediction <- predict(nb_model, task = validation_task)
rf_validation_prediction <- predict(rf_model, task = validation_task)
gb_validation_prediction <- predict(gb_model, task = validation_task)
svm_validation_prediction <- predict(svm_model, task = validation_task)

saveRDS(rs_validation_prediction, rs_validation_prediction_rds_file)
saveRDS(nb_validation_prediction, nb_validation_prediction_rds_file)
saveRDS(rf_validation_prediction, rf_validation_prediction_rds_file)
saveRDS(gb_validation_prediction, gb_validation_prediction_rds_file)
saveRDS(svm_validation_prediction, svm_validation_prediction_rds_file)

# rs_validation_prediction <- readRDS(rs_validation_prediction_rds_file)
# nb_validation_prediction <- readRDS(nb_validation_prediction_rds_file)
# rf_validation_prediction <- readRDS(rf_validation_prediction_rds_file)
# gb_validation_prediction <- readRDS(gb_validation_prediction_rds_file)
# svm_validation_prediction <- readRDS(svm_validation_prediction_rds_file)

# Performance ----
rs_internal_validation_performance <- performance(rs_validation_prediction,
                                                  measures = measures)
nb_internal_validation_performance <- performance(nb_validation_prediction,
                                                  measures = measures)
rf_internal_validation_performance <- performance(rf_validation_prediction,
                                                  measures = measures)
gb_internal_validation_performance <- performance(gb_validation_prediction,
                                                  measures = measures)
svm_internal_validation_performance <- performance(svm_validation_prediction,
                                                  measures = measures)

# Plots ----
# Plot false negative and positive rates as well as the error rate versus the threshold
d = generateThreshVsPerfData(list(riskScorer = rs_validation_prediction,
                                  naiveBayes = nb_validation_prediction,
                                  randomForest = rf_validation_prediction,
                                  gradientBoosting = gb_validation_prediction,
                                  svm = svm_validation_prediction),
                             measures = measures,
                             gridsize = 1000)
plotThreshVsPerf(d)
plotROCCurves(d, measures = list(mlr::fpr, mlr::tpr))

ggplot(d$data) +
  geom_line(aes(x = fpr, y = tpr, color = learner)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.5) +
  xlab("Falsch Positiv Rate (1-Spezifität)") + ylab("Richtig Positiv Rate (Sensitivität)")
