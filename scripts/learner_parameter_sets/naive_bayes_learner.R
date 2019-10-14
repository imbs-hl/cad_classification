
# Define Naive Bayes learner parameter set ----
nb_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeNumericParam("laplace",
                                 lower = 0,
                                 upper = 10)
)
