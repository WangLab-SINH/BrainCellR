#' @title Run consensus1
#' @description  Using consensus1 method to clust single cell data.
#' @param norm.dat Single cell expression data input.
#' @param select.cells Numeric, use all or part of the cells for clustering, cells can be sampled to save running time in multi-round clustering. Using all cells by default.
#' @param k.nn K-nearst neighbors used for clustering. Default 15.
#' @param max.dim The number of PCA dimensions retained. Default 20.
#' @param output_dir File directory to store clustering results.
#' @param mc.cores Number of cores used for parallel processing.
#' @param de.param The differential gene expression threshold.
#' @param merge.type Determine if the DE gene score threshold should be applied to combined de.score, or de.score for up and down directions separately.
#' @param cut.method Clustering method. It can be "auto", "louvain", "hclust".
#' @param override Whether the file should be overrided.
#' @param init.result If set, the function will only find finer splits of the current clusters.
#' @param ... Other parameters.
#' @return A list with two elements: cl: cluster membership for each cell markers: top markers that seperate clusters
#' @export
#' @import scrattch.hicat
RunConsensus1 <- function (norm.dat, select.cells = colnames(norm.dat), k.nn=15, max.dim=20,
                           output_dir = "subsample_result",
                           mc.cores = 1, de.param = de_param(), merge.type = c("undirectional",
                                                                               "directional"), override = FALSE, init.result = NULL,
                           cut.method = "auto", ...) {
  niter = 1
  sampleSize = length(select.cells)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  all.cells = select.cells
  if (!is.null(init.result)) {
    all.cells = intersect(all.cells, names(init.result$cl))
  }
  run <- function(i, ...) {
    prefix = paste("iter", i, sep = ".")
    print(prefix)
    outfile = file.path(output_dir, paste0("result.", i,
                                           ".rda"))
    if (file.exists(outfile) & !override) {
      return(NULL)
    }
    select.cells = all.cells
    save(select.cells, file = file.path(output_dir, paste0("cells.",
                                                           i, ".rda")))
    result <- IterCluster(norm.dat = norm.dat,
                          select.cells = select.cells, prefix = prefix, de.param = de.param, k.nn=k.nn, max.dim=max.dim,
                          merge.type = merge.type, result = init.result)
    save(result, file = outfile)
  }

  if (mc.cores == 1) {
    sapply(1:niter, function(i) {
      run(i, ...)
    })
  }
  else {
    require(foreach)
    require(doParallel)
    cl <- makeForkCluster(mc.cores)
    registerDoParallel(cl)
    foreach(i = 1:niter, .combine = "c") %dopar% run(i)
    stopCluster(cl)
  }
  result.files = file.path(output_dir, dir(output_dir,
                                           "result.*.rda"))
  load(result.files)
  return(result)
}
