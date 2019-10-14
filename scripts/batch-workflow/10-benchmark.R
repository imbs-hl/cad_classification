
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
mbo_ctrl <- mlrMBO::makeMBOControl(impute.y.fun = function(x, y, opt.path, ...) -0.5) # This is the worst AUC
mbo_ctrl <- mlrMBO::setMBOControlTermination(mbo_ctrl, iters = bm_tuning_iters, time.budget = bm_time_budget)
surrogate_lrn <- mlr::makeImputeWrapper(mlr::makeLearner("regr.ranger", predict.type = "se", replace = FALSE),
                                        classes = list(numeric = mlr::imputeMax(2),
                                                       factor = mlr::imputeConstant("__miss__")))
ctrl <- mlr:::makeTuneControlMBO(learner = surrogate_lrn,
                                 mbo.control = mbo_ctrl)

# Define resampling strategies ----
inner <- mlr::makeResampleDesc("CV", stratify = bm_stratify, iters = bm_inner_cv_iters)
outer <- mlr::makeResampleDesc("CV", stratify = bm_stratify, iters = bm_outer_cv_iters)

# Define random forest learner ----
rf_lrn <- mlr::makeLearner("classif.ranger", predict.type = "prob", config = list(on.learner.error = "warn"))
rf_lrn <- mlr::makeFilterWrapper(learner = rf_lrn, fw.method = "ranger.filter", reduction = "pruned")
rf_lrn <- mlr::makeDownsampleWrapper(rf_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

rf_lrn <- mlr::makeTuneWrapper(rf_lrn,
                               resampling = inner,
                               control = ctrl,
                               measures = measures,
                               par.set = c(filter_ps, rf_ps),
                               show.info = TRUE)

# Define MB-MDR learner ----
mdr_lrn <- mlr::makeLearner("classif.mbmdrc", predict.type = "prob",
                            order = 2L, order.range = TRUE,
                            cv.top.results = 5000L, max.results = 5000L,
                            folds = 3L, cv.loss = "auc",
                            config = list(on.learner.error = "warn"))
mdr_lrn <- mlr::makePreprocWrapper(mdr_lrn,
                                   train = round_numerics_training,
                                   predict = round_numerics_prediction,
                                   par.vals = list())
mdr_lrn <- mlr::makeFilterWrapper(learner = mdr_lrn, fw.method = "ranger.filter")
mdr_lrn <- mlr::makeDownsampleWrapper(mdr_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

mdr_lrn <- mlr::makeTuneWrapper(mdr_lrn,
                                resampling = inner,
                                control = ctrl,
                                measures = measures,
                                par.set = c(filter_ps, mdr_ps),
                                show.info = TRUE)

# Define risk scorer learner ----
rs_lrn <- mlr::makeLearner("classif.riskScorer", importance = gwas_importance, predict.type = "prob",
                           config = list(on.learner.error = "warn"))
rs_lrn <- mlr::makeFilterWrapper(learner = rs_lrn, fw.method = "ranger.filter")
rs_lrn <- mlr::makeDownsampleWrapper(rs_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

rs_lrn <- mlr::makeTuneWrapper(rs_lrn,
                               resampling = inner,
                               control = ctrl,
                               measures = measures,
                               par.set = c(filter_ps, rs_ps),
                               show.info = TRUE)

# Define SVM learner ----
svm_lrn <- mlr::makeLearner("classif.ksvm", predict.type = "prob", config = list(on.learner.error = "warn"))
svm_lrn <- mlr::makeFilterWrapper(svm_lrn, fw.method = "ranger.filter")
svm_lrn <- mlr::makeDownsampleWrapper(svm_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

svm_lrn <- mlr::makeTuneWrapper(svm_lrn,
                                resampling = inner,
                                control = ctrl,
                                measures = measures,
                                par.set = c(filter_ps, svm_ps),
                                show.info = TRUE)

# Define GLM learner ----
glm_lrn <- mlr::makeLearner("classif.binomial", predict.type = "prob",
                            config = list(on.learner.error = "warn"))
glm_lrn <- mlr::makeFilterWrapper(glm_lrn, fw.method = "ranger.filter")
glm_lrn <- mlr::makeDownsampleWrapper(glm_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

# Define filter parameter set for GLM
filter_ps_glm <- ParamHelpers::makeParamSet(
  ParamHelpers::makeIntegerParam("fw.abs",
                                 lower = 1,
                                 upper = ceiling(mlr::getTaskSize(training_task)*bm_dw))
)

glm_lrn <- mlr::makeTuneWrapper(glm_lrn,
                                resampling = inner,
                                control = ctrl,
                                measures = measures,
                                par.set = c(filter_ps_glm, glm_ps),
                                show.info = TRUE)

# Define GLMNET learner ----
enet_lrn <- mlr::makeLearner("classif.cvglmnet", predict.type = "prob",
                             type.measure = "auc",
                             nlambda = 1000L,
                             lambda.min.ratio = 1e-5,
                             nfolds = 5,
                             config = list(on.learner.error = "warn"))
enet_lrn <- mlr::makeFilterWrapper(enet_lrn, fw.method = "ranger.filter")
enet_lrn <- mlr::makeDownsampleWrapper(enet_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

enet_lrn <- mlr::makeTuneWrapper(enet_lrn,
                                 resampling = inner,
                                 control = ctrl,
                                 measures = measures,
                                 par.set = c(filter_ps, enet_ps),
                                 show.info = TRUE)

# Define Naive Bayes learner ----
nb_lrn <- mlr::makeLearner("classif.naiveBayes", predict.type = "prob",
                           config = list(on.learner.error = "warn"))
nb_lrn <- mlr::makeFilterWrapper(nb_lrn, fw.method = "ranger.filter")
nb_lrn <- mlr::makeDownsampleWrapper(nb_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

nb_lrn <- mlr::makeTuneWrapper(nb_lrn,
                               resampling = inner,
                               control = ctrl,
                               measures = measures,
                               par.set = c(filter_ps, nb_ps),
                               show.info = TRUE)

# Define gradient boosting learner ----
gb_lrn <- mlr::makeLearner("classif.xgboost", predict.type = "prob",
                           print_every_n = 1000L,
                           nthread = 1L,
                           save_name = "/dev/null",
                           config = list(on.learner.error = "warn"))
gb_lrn <- mlr::makeFilterWrapper(gb_lrn, fw.method = "ranger.filter")
gb_lrn <- mlr::makeDownsampleWrapper(gb_lrn, dw.perc = bm_dw, dw.stratify = bm_stratify)

gb_lrn <- mlr::makeTuneWrapper(gb_lrn,
                               resampling = inner,
                               control = ctrl,
                               measures = measures,
                               par.set = c(filter_ps, gb_ps),
                               show.info = TRUE)

# Submit to cluster ----
reg <- batchtools::makeExperimentRegistry(
  conf.file = slurm_btconf_file,
  packages = c("ranger", "riskScoreR",
               "MBMDRClassifieR",
               "kernlab", "e1071",
               "glmnet", "xgboost",
               "mlr", "mlrMBO",
               "gwasFilter", "importanceFilter"),
  file.dir = file.path(reg_dir, "batchmark"),
  work.dir = file.path(proc_dir)
)

mlr::batchmark(
  tasks = training_task,
  learners = list(mdr_lrn,
                  svm_lrn,
                  rf_lrn,
                  gb_lrn,
                  enet_lrn,
                  glm_lrn,
                  nb_lrn,
                  rs_lrn),
  resamplings = outer,
  measures = measures,
  models = FALSE
)

batch_ids <- findNotDone()
batch_ids <- ajoin(batch_ids, findRunning())
batch_ids[, chunk := 1]
batchtools::submitJobs(
  ids = batch_ids,
  resources = list(name = "batchmark",
                   ntasks = 1,
                   ncpus = 1,
                   partition = "prio,batch",
                   walltime = 0L,
                   memory = "60G",
                   account = "bayer",
                   chunks.as.arrayjobs = TRUE)
)

batchtools::waitForJobs()

# Reduce results ----
na.rm <- TRUE
bmr <- reduceBatchmarkResults()

saveRDS(bmr, bmr_rds_file)
