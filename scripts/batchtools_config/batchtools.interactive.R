external <- FALSE
cluster.functions = batchtools::makeClusterFunctions(name = "InteractiveHooks",
                                                     submitJob = function(reg, jc) {
                                                       assertRegistry(reg, writeable = TRUE)
                                                       assertClass(jc, "JobCollection")
                                                       if (external) {
                                                         runOSCommand(Rscript(), sprintf("-e \"batchtools::doJobCollection('%s', output = '%s')\"",
                                                                                         jc$uri, jc$log.file))
                                                         makeSubmitJobResult(status = 0L, batch.id = "cfInteractive")
                                                       }
                                                       else {
                                                         doJobCollection(jc, output = jc$log.file)
                                                         makeSubmitJobResult(status = 0L, batch.id = "cfInteractive")
                                                       }
                                                     },
                                                     store.job = FALSE, fs.latency = NA_real_, hooks = list(pre.do.collection = function(reg, cache, ...) {
                                                       parallelMap::parallelLibrary(packages = reg$packages)
                                                       parallelMap::parallelSource(files = reg$source)}
                                                     ))
