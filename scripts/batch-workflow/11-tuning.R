
source("init.R", chdir = TRUE)

# Load data ----
gwas_importance <- readRDS(gwas_importance_rds_file)
training_task <- readRDS(training_task_rds_file)

# Define measures ----
measures <- list(mlr::auc, # use AUC as performance measure during tuning
                 mlr::mmce, mlr::npv, mlr::fpr, mlr::f1, mlr::fnr, mlr::ssr,
                 mlr::tp, mlr::tn, mlr::gpr, mlr::lsr, mlr::acc,
                 mlr::wkappa, mlr::ppv, mlr::logloss, mlr::ber, mlr::tpr,
                 mlr::brier, mlr::gmean, mlr::fdr, mlr::tnr, mlr::qsr, mlr::bac,
                 mlr::brier.scaled, mlr::fp, mlr::fn, mlr::kappa)

# Define filter parameter set ----
filter_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeIntegerParam("fw.abs",
                                 lower = 1,
                                 upper = mlr::getTaskNFeats(training_task))
)

# Define MBO control ----
mbo_ctrl <- mlrMBO::makeMBOControl(
  save.on.disk.at = 0:(tuning_iters+1),
  impute.y.fun = function(x, y, opt.path, ...) -0.5 # This is the worst AUC
)
mbo_ctrl <- mlrMBO::setMBOControlTermination(
  mbo_ctrl,
  iters = tuning_iters,
  time.budget = tuning_time_budget
)
surrogate_lrn <- mlr::makeImputeWrapper(
  learner = mlr::makeLearner("regr.ranger", predict.type = "se", replace = FALSE),
  classes = list(numeric = mlr::imputeMax(2),
                 factor = mlr::imputeConstant("__miss__"))
)
ctrl <- mlr:::makeTuneControlMBO(learner = surrogate_lrn,
                                 mbo.control = mbo_ctrl)

# Define resampling strategies ----
resampling_description <- mlr::makeResampleDesc("CV",
                                                stratify = tuning_stratify,
                                                iters = tuning_cv_iters)

# Define risk scorer learner ----
rs_ctrl <- ctrl
rs_ctrl$mbo.control$save.file.path <- file.path(proc_dir, "rs_mbo_run.RDdata")
rs_lrn <- mlr::makeLearner("classif.riskScorer",
                           importance = gwas_importance, predict.type = "prob",
                           config = list(on.learner.error = "warn"))
rs_lrn <- mlr::makeFilterWrapper(learner = rs_lrn, fw.method = "ranger.filter")

rs_lrn <- mlr::makeTuneWrapper(rs_lrn,
                               resampling = resampling_description,
                               control = rs_ctrl,
                               measures = measures,
                               par.set = c(filter_ps, rs_ps),
                               show.info = TRUE)

# Define SVM learner ----
svm_ctrl <- ctrl
svm_ctrl$mbo.control$save.file.path <- file.path(proc_dir, "svm_mbo_run.RDdata")
svm_lrn <- mlr::makeLearner("classif.ksvm", predict.type = "prob", config = list(on.learner.error = "warn"))
svm_lrn <- mlr::makeFilterWrapper(svm_lrn, fw.method = "ranger.filter")

svm_lrn <- mlr::makeTuneWrapper(svm_lrn,
                                resampling = resampling_description,
                                control = svm_ctrl,
                                measures = measures,
                                par.set = c(filter_ps, svm_ps),
                                show.info = TRUE)

# Define RandomForest learner ----
rf_ctrl <- ctrl
rf_ctrl$mbo.control$save.file.path <- file.path(proc_dir, "rf_mbo_run.RData")
rf_lrn <- mlr::makeLearner("classif.ranger", predict.type = "prob", config = list(on.learner.error = "warn"))
rf_lrn <- mlr::makeFilterWrapper(rf_lrn, fw.method = "ranger.filter")

rf_lrn <- mlr::makeTuneWrapper(rf_lrn,
                               resampling = resampling_description,
                               control = rf_ctrl,
                               measures = measures,
                               par.set = c(filter_ps, rf_ps),
                               show.info = TRUE)

# Define Naive Bayes learner ----
nb_ctrl <- ctrl
nb_ctrl$mbo.control$save.file.path <- file.path(proc_dir, "nb_mbo_run.RData")
nb_lrn <- mlr::makeLearner("classif.naiveBayes", predict.type = "prob",
                           config = list(on.learner.error = "warn"))
nb_lrn <- mlr::makeFilterWrapper(nb_lrn, fw.method = "ranger.filter")

nb_lrn <- mlr::makeTuneWrapper(nb_lrn,
                               resampling = resampling_description,
                               control = nb_ctrl,
                               measures = measures,
                               par.set = c(filter_ps, nb_ps),
                               show.info = TRUE)

# Define XGBoost learner ----
gb_ctrl <- ctrl
gb_ctrl$mbo.control$save.file.path <- file.path(proc_dir, "gb_mbo_run.RData")
gb_lrn <- mlr::makeLearner("classif.xgboost", predict.type = "prob",
                                           print_every_n = 1000L,
                                           nthread = 1L,
                                           save_name = "/dev/null",
                                           config = list(on.learner.error = "warn"))
gb_lrn <- mlr::makeFilterWrapper(gb_lrn, fw.method = "ranger.filter")

gb_lrn <- mlr::makeTuneWrapper(gb_lrn,
                               resampling = resampling_description,
                               control = gb_ctrl,
                               measures = measures,
                               par.set = c(filter_ps, gb_ps),
                               show.info = TRUE)

# Start batchtools backend ----
parallelMap::parallelStartBatchtools(
  logging = TRUE,
  storagedir = reg_dir,
  bt.resources = list(name = "tuning",
                      ntasks = 1,
                      ncpus = 1,
                      partition = "prio,batch",
                      walltime = 0L,
                      memory = "30G",
                      account = "bayer",
                      chunks.as.arrayjobs = TRUE),
  n.chunks = 1,
  work.dir = proc_dir,
  level = "mlr.resample"
)

# Tune risk scorer ----
rs_tuning <- train(
  learner = rs_lrn,
  task = training_task
)
saveRDS(rs_tuning, rs_model_rds_file)

# Tune SVM ----
svm_tuning <- train(
  learner = svm_lrn,
  task = training_task
)
saveRDS(svm_tuning, svm_model_rds_file)

# Tune Naive Bayes ----
nb_tuning <- train(
  learner = nb_lrn,
  task = training_task
)
saveRDS(nb_tuning, nb_model_rds_file)

# Tune RandomForest ----
rf_tuning <- train(
  learner = rf_lrn,
  task = training_task
)
saveRDS(rf_tuning, rf_model_rds_file)

# Tune XGBoost ----
gb_tuning <- train(
  learner = gb_lrn,
  task = training_task
)
saveRDS(gb_tuning, gb_model_rds_file)

# Clean up ----
parallelMap::parallelStop()
