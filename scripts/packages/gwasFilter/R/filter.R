#' @import mlr data.table
.onLoad <- function(libname, pkgname) {
  mlr::makeFilter(
    name = "gwas.filter",
    desc = "Returns -log10 P values of CAD GWAS",
    pkg = "gwasFilter",
    supported.tasks = c("classif", "regr", "surv"),
    supported.features = c("numerics", "factors", "ordered"),
    fun = function(task, nselect, ...) {
      SNP = NULL
      P = NULL
      feats = mlr::getTaskFeatureNames(task)
      gwas_results[make.names(SNP) %in% feats,
                      stats::setNames(-log10(P_GWAS), make.names(SNP))]
    }
  )
}
