
# Define GLMNET learner parameter set ----
enet_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeLogicalParam("standardize"),
  ParamHelpers::makeDiscreteParam("s",
                                  values = c("lambda.1se", "lambda.min")),
  ParamHelpers::makeNumericParam("alpha",
                                 lower = 0,
                                 upper = 1)
)
