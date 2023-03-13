#' @title Iteractive clustering
#' @description  Iterative clustering algorithm for single cell RNAseq dataset.
#' @param norm.dat data input
#' @param select.cells data ref
#' @param prefix meta input
#' @param k.nn meta input
#' @param max.dim meta input
#' @param split.size meta input
#' @param result meta input
#' @param method meta input
#' @param de.param The differential gene expression threshold.
#' @param merge.type Determine if the DE gene score threshold should be applied to combined de.score, or de.score for up and down directions separately.
#' @return A numeric vector of log density.
#' @import scrattch.hicat
#' @export
IterCluster <- function (norm.dat, select.cells = colnames(norm.dat), prefix = NULL, k.nn, max.dim, de.param = de_param(), merge.type = c("undirectional",
                                                                                                                                          "directional"),
                         split.size = 10, result = NULL, method = "auto")
{
  sampleSize = length(select.cells)
  if (!is.null(prefix)) {
    print(prefix)
  }
  if (method == "auto") {
    if (length(select.cells) > 3000) {
      select.method = "louvain"
    }
    else {
      select.method = "ward.D"
    }
  }
  else {
    select.method = method
  }
  if (length(select.cells) <= 3000) {
    if (!is.matrix(norm.dat)) {
      norm.dat = as.matrix(norm.dat[, select.cells])
    }
  }
  if (is.null(result)) {
    result = OneStepCluster(norm.dat, select.cells = select.cells, k.nn=k.nn, max.dim=max.dim,de.param = de.param, merge.type = c("undirectional",
                                                                                                                                    "directional"),
                            prefix = prefix, method = select.method)
    gc()
  }
  if (!is.null(result)) {
    select.cells = intersect(select.cells, names(result$cl))
    cl = result$cl[select.cells]
    gene.mod = result$gene.mod
    markers = result$markers
    cl = setNames(as.integer(cl), names(cl))
    new.cl = cl
    cl.size = table(cl)
    to.split = names(cl.size)[cl.size >= split.size]
    if (length(to.split) > 0) {
      n.cl = 1
      for (x in sort(unique(cl))) {
        tmp.cells = names(cl)[cl == x]
        if (!x %in% to.split) {
          new.cl[tmp.cells] = n.cl
        }
        else {
          tmp.prefix = paste(prefix, x, sep = ".")
          tmp.result = IterCluster(norm.dat = norm.dat,
                                   select.cells = tmp.cells, prefix = tmp.prefix,de.param = de.param, merge.type = c("undirectional",
                                                                                                                       "directional"),
                                   max.dim=max.dim, k.nn=k.nn,
                                   split.size = split.size, method = method)
          gc()
          if (is.null(tmp.result)) {
            new.cl[tmp.cells] = n.cl
          }
          else {
            tmp.cl = tmp.result$cl
            if (length(unique(tmp.cl) > 1)) {
              new.cl[names(tmp.cl)] = n.cl + as.integer(tmp.cl)
              markers = union(markers, tmp.result$markers)
            }
          }
        }
        n.cl = max(new.cl) + 1
      }
      cl = new.cl
    }
    result = list(cl = cl, markers = markers)
    return(result)
  }
  return(NULL)
}
