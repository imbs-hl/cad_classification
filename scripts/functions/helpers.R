
#' Load or create a batchtools registry
#'
#' Tries to load the batchtools registry in \code{proc_dir} with the given name.
#' If failing, tries to create a new batchtools registry.
#'
#' @param reg.name [\code{string}]\cr
#'                 The name of the batchtools registry.
#' @param packages [\code{character}]\cr
#'                 Vector of package names to load for each job.
#'
#'
#' @return A batchtools \code{\link{batchtools::Registry}} object.
#' @export
#'
load_or_create_registry <- function(reg.name, packages = character(0L),
                                    conf.file) {
  tryCatch(batchtools::loadRegistry(file.dir = file.path(reg_dir, reg.name),
                                    work.dir = proc_dir,
                                    conf.file = conf.file),
           error = function(e) {
             warning(sprintf("Error loading registry: %s.\nCreating a new registry.", as.character(e)))
             flush.console()
             tryCatch(batchtools::makeRegistry(file.dir = file.path(reg_dir, reg.name),
                                               work.dir = proc_dir,
                                               conf.file = conf.file,
                                               seed = seed,
                                               packages = c("checkmate", packages)),
                      error = function(e) stop(e))
           })
}

round_numerics_training <- function(data, target, args = list()) {
  ## Save args for transformation during prediction
  control = args
  ## Round numerical features from the data set
  data = cbind(data[target], lapply(data[, -which(names(data) == target)], FUN = function(col) {
    if (is.numeric(col)) {
      return(round(col))
    } else {
      return(col)
    }
  }))
  return(list(data = data, control = control))
}

round_numerics_prediction <- function(data, target, args = list(), control = list()) {
  ## Round numerical features from the data set
  data = as.data.frame(lapply(data, FUN = function(col) {
    if (is.numeric(col)) {
      return(round(col))
    } else {
      return(col)
    }
  }))
  return(data)
}
