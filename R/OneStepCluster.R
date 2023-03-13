#' @title One step cluster
#' @description  Creating poisson vector.
#' @details Input an integer and return the log density of a poisson distribution with lambda equals the input integer
#' @param norm.dat data input
#' @param select.cells data ref
#' @param counts meta input
#' @param method meta input
#' @param vg.padj.th meta input
#' @param dim.method meta input
#' @param max.dim meta input
#' @param rm.eigen meta input
#' @param rm.th meta input
#' @param de.param meta input
#' @param merge.type meta input
#' @param maxGenes meta input
#' @param maxGenes meta input
#' @param sampleSize meta input
#' @param k.nn meta input
#' @param prefix meta input
#' @param verbose meta input
#' @param regress.x meta input
#' @param max.cl.size meta input
#' @return A numeric vector of log density.
#' @export
#' @import WGCNA
#' @import scrattch.hicat
OneStepCluster <- function (norm.dat, select.cells = colnames(norm.dat), counts = NULL,
                            method = c("louvain", "leiden", "ward.D", "kmeans"), vg.padj.th = 0.5,
                            dim.method = c("pca", "WGCNA"), max.dim = 20, rm.eigen = NULL,
                            rm.th = 0.7, de.param = de_param(), merge.type = c("undirectional",
                                                                               "directional"), maxGenes = 3000, sampleSize = 40000,
                            max.cl.size = 500, k.nn = 50, prefix = NULL, verbose = FALSE,
                            regress.x = NULL)
{
  library(matrixStats)
  method <- match.arg(method)
  dim.method <- match.arg(dim.method)
  merge.type <- match.arg(merge.type)
  if (!is.null(regress.x)) {
    print("regression")
    tmp = lm_normalize(as.matrix(norm.dat[, select.cells]),
                       regress.x[select.cells], R_2.th = 0.1)
    norm.dat = tmp[[1]]
  }
  if (length(select.cells) > sampleSize) {
    sampled.cells = sample(select.cells, pmin(length(select.cells),
                                              sampleSize))
  }
  else {
    sampled.cells = select.cells
  }
  select.genes = row.names(norm.dat)[which(Matrix::rowSums(norm.dat[,
                                                                    select.cells] > de.param$low.th) >= de.param$min.cells)]
  if (is.null(counts)) {
    if (is.matrix(norm.dat)) {
      counts = 2^(norm.dat[, sampled.cells]) - 1
    }
    else {
      counts = norm.dat[, sampled.cells]
      counts@x = 2^(counts@x) - 1
    }
  }
  plot_file = NULL
  if (verbose & !is.null(prefix)) {
    plot_file = paste0(prefix, ".vg.pdf")
  }
  vg = find_vg(as.matrix(counts[select.genes, sampled.cells]),
               plot_file = plot_file)
  if (dim.method == "auto") {
    if (length(select.cells) > 1000) {
      dim.method = "pca"
    }
    else {
      dim.method = "WGCNA"
    }
  }
  if (dim.method == "WGCNA") {
    select.genes = as.character(vg[which(vg$loess.padj <
                                           1), "gene"])
    select.genes = head(select.genes[order(vg[select.genes,
                                              "loess.padj"], -vg[select.genes, "z"])], maxGenes)
    rd.dat = rd_WGCNA(norm.dat, select.genes = select.genes,
                      select.cells = select.cells, sampled.cells = sampled.cells,
                      de.param = de.param, max.mod = max.dim, max.cl.size = max.cl.size)$rd.dat
  }
  else {
    select.genes = as.character(vg[which(vg$loess.padj <
                                           vg.padj.th | vg$dispersion > 3), "gene"])
    select.genes = head(select.genes[order(vg[select.genes,
                                              "loess.padj"], -vg[select.genes, "z"])], maxGenes)
    if (verbose) {
      cat("Num high variance genes:", length(select.genes),
          "\n")
    }
    if (length(select.genes) < de.param$min.genes) {
      return(NULL)
    }
    rd.dat = rd_PCA(norm.dat, select.genes, select.cells,
                    sampled.cells = sampled.cells, max.pca = max.dim)$rd.dat
  }
  if (is.null(rd.dat) || ncol(rd.dat) == 0) {
    return(NULL)
  }
  if (!is.null(rm.eigen)) {
    rd.dat <- filter_RD(rd.dat, rm.eigen, rm.th, verbose = verbose)
  }
  if (is.null(rd.dat) || ncol(rd.dat) == 0) {
    return(NULL)
  }
  if (verbose) {
    print(method)
  }
  max.cl = ncol(rd.dat) * 2 + 1
  if (method == "louvain") {
    k = pmin(k.nn, round(nrow(rd.dat)/2))
    tmp = jaccard_louvain(rd.dat, k)
    if (is.null(tmp)) {
      return(NULL)
    }
    cl = tmp$cl
    if (length(unique(cl)) > max.cl) {
      tmp.means = do.call("cbind", tapply(names(cl), cl,
                                          function(x) {
                                            colMeans(rd.dat[x, , drop = F])
                                          }, simplify = F))
      tmp.hc = hclust(dist(t(tmp.means)), method = "average")
      tmp.cl = cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
      cl = setNames(tmp.cl[as.character(cl)], names(cl))
    }
  }
  else if (method == "leiden") {
    k = pmin(k.nn, round(nrow(rd.dat)/2))
    tmp = jaccard_leiden(rd.dat, k)
    if (is.null(tmp)) {
      return(NULL)
    }
    cl = tmp$cl
    if (length(unique(cl)) > max.cl) {
      tmp.means = do.call("cbind", tapply(names(cl), cl,
                                          function(x) {
                                            colMeans(rd.dat[x, , drop = F])
                                          }, simplify = F))
      tmp.hc = hclust(dist(t(tmp.means)), method = "average")
      tmp.cl = cutree(tmp.hc, pmin(max.cl, length(unique(cl))))
      cl = setNames(tmp.cl[as.character(cl)], names(cl))
    }
  }
  else if (method == "ward.D") {
    hc = hclust(dist(rd.dat), method = "ward.D")
    cl = cutree(hc, max.cl)
  }
  else if (method == "kmeans") {
    cl = kmeans(rd.dat, max.cl)$cluster
  }
  else {
    stop(paste("Unknown clustering method", method))
  }
  rd.dat.t = t(rd.dat)
  merge.result = merge_cl(norm.dat, cl = cl, rd.dat.t = rd.dat.t,
                          merge.type = merge.type, de.param = de.param, max.cl.size = max.cl.size,
                          verbose = verbose)
  gc()
  if (is.null(merge.result))
    return(NULL)
  sc = merge.result$sc
  cl = merge.result$cl
  if (length(unique(cl)) > 1) {
    if (verbose) {
      cat("Expand", prefix, "\n")
      cl.size = table(cl)
      print(cl.size)
      save(cl, file = paste0(prefix, ".cl.rda"))
    }
    de.genes = merge.result$de.genes
    markers = merge.result$markers
    cl.dat = get_cl_means(norm.dat[markers, ], cl[sample_cells(cl,
                                                               max.cl.size)])
    cl.hc = hclust(dist(t(cl.dat)), method = "average")
    cl = setNames(factor(as.character(cl), levels = colnames(cl.dat)[cl.hc$order]),
                  names(cl))
    if (verbose & !is.null(prefix)) {
      tmp = display_cl(cl, norm.dat, prefix = prefix,
                       markers = markers, max.cl.size = max.cl.size)
    }
    levels(cl) = 1:length(levels(cl))
    result = list(cl = cl, markers = markers)
    return(result)
  }
  return(NULL)
}
