
# Define GLM learner parameter set ----
glm_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeDiscreteParam("link",
                                  values = c("logit",
                                             "probit",
                                             "cloglog"))
)
