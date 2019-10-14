#' @import mlr data.table
.onLoad <- function(libname, pkgname) {
  mlr::makeFilter(
    name = "ranger.filter",
    desc = "Returns precomputed importance values from ranger",
    pkg = "importanceFilter",
    supported.tasks = c("classif", "regr", "surv"),
    supported.features = c("numerics", "factors", "ordered"),
    fun = function(task, nselect, ...) {
      SNP = NULL
      P = NULL
      dots <- list(...)
      importance <- if (!is.null(dots$reduction)) {
        switch(dots$reduction,
               clumped = clumped_importance,
               pruned_importance)
      } else {
        pruned_importance
      }
      feats = mlr::getTaskFeatureNames(task)
      importance[make.names(SNP) %in% feats,
                 stats::setNames(importance, make.names(SNP))]
    }
  )
}
