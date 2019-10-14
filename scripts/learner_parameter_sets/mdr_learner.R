
# Define MB-MDR learner parameter set ----
mdr_ps <- ParamHelpers::makeParamSet(
  ParamHelpers::makeLogicalParam("o.as.na"),
  ParamHelpers::makeIntegerParam("min.cell.size",
                                 lower = 1L,
                                 upper = 50L),
  ParamHelpers::makeNumericParam("alpha",
                                 lower = 0,
                                 upper = 1),
  ParamHelpers::makeDiscreteParam("adjustment",
                                  values = c("NONE", "CODOMINANT"))
)
