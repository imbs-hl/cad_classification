
# Define random forest learner parameter set ----
rf_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeIntegerParam("num.trees",
                                 lower = 100L,
                                 upper = 5000L,
                                 trafo = function(x) plyr::round_any(x, 100)),
  ParamHelpers::makeNumericParam("mtry.perc",
                                 lower = 0.001,
                                 upper = 0.1),
  ParamHelpers::makeIntegerParam("min.node.size",
                                 lower = 10L,
                                 upper = 100L),
  ParamHelpers::makeLogicalParam("replace")
)
