#' @title Find differential expressed gene
#' @description  Find differential expressed gene.
#' @param new_data Data input, single cell expression matrix.
#' @param data_meta Cell type classification result from RunSubclassClassify.
#' @param cell_number Cell type classification result from RunSubclassClassify.
#' @return A list with three files contains DE gene reuslt.
#' @export
#' @import Seurat
RunFindDEGene <- function(new_data, data_meta, cell_number=10)
{

  mouse_data <- CreateSeuratObject(counts = new_data, min.cells = 0, min.features = 0, project = "example")
  mouse_data <- AddMetaData(mouse_data, data_meta)
  mouse_data <- NormalizeData(mouse_data, normalization.method = "LogNormalize", scale.factor = 10000)
  Idents(mouse_data) <- mouse_data$subclass_label
  all_marker_list = list()

  data1_anno = data_meta

  current_sample_data <- t(new_data)
  current_sample_meta <- data1_anno

  set.seed(1)
  type_list <- unique(current_sample_meta$cluster_label)
  new_sample_data <- current_sample_data[1,]
  new_meta_data <- c()
  for(i in 1:length(type_list)){
    temp <- current_sample_data[current_sample_meta$cluster_label==type_list[i],]
    index_id <- 1:nrow(temp)
    index_id <- sample(index_id, length(index_id))
    for(j in 1:ceiling(length(index_id)/cell_number)){
      new_meta_data <- c(new_meta_data, type_list[i])
      if(j*cell_number > length(index_id)){
        if(((j-1)*cell_number+1)==length(index_id)){
          new_sample_data <- rbind(new_sample_data, as.numeric(temp[index_id[(1+(j-1)*cell_number):length(index_id)],]))
        }else{
          new_sample_data <- rbind(new_sample_data, colMeans(temp[index_id[(1+(j-1)*cell_number):length(index_id)],]))
        }

      }else{
        new_sample_data <- rbind(new_sample_data, colMeans(temp[index_id[(1+(j-1)*cell_number):(j*cell_number)],]))
      }
    }
  }
  new_sample_data <- new_sample_data[-1,]

  rownames(data_meta) = data_meta$cell
  data1_anno = data_meta
  rownames(new_sample_data) <- 1:nrow(new_sample_data)
  new_meta = data.frame(1:nrow(new_sample_data), new_meta_data)
  colnames(new_meta) = c("cell","cluster")
  temp_subclass <- c()
  for(i in 1:nrow(new_meta)){
    temp_subclass <- c(temp_subclass, unique(data1_anno[data1_anno$cluster_label==new_meta$cluster[i],]$subclass_label))
  }
  new_meta <- cbind(new_meta, temp_subclass)
  colnames(new_meta) <- c("cell","cluster_label","subclass_label")
  rownames(new_meta) <- new_meta$cell
  new_sample_data = t(new_sample_data)
  mouse_data <- CreateSeuratObject(counts = new_sample_data, min.cells = 0, min.features = 0, project = "example")
  mouse_data <- AddMetaData(mouse_data, new_meta)
  mouse_data <- NormalizeData(mouse_data, normalization.method = "LogNormalize", scale.factor = 10000)
  Idents(mouse_data) <- mouse_data$subclass_label
  all_marker_list = list()
  num = 1
  for(i in unique(mouse_data$subclass_label)){
    temp_mouse_data = subset(mouse_data, idents = i)
    Idents(temp_mouse_data) <- temp_mouse_data$cluster_label
    #plan(workers = 6)
    mouse_cells_markers <- FindAllMarkers(temp_mouse_data, test.use = "roc",densify=T)
    mouse_cells_markers = mouse_cells_markers[mouse_cells_markers$avg_log2FC>0,]
    all_marker_list[[num]] = mouse_cells_markers
    num = num + 1
  }
  names(all_marker_list) = unique(mouse_data$subclass_label)

  rownames(new_meta) = new_meta$cell
  new_meta$class = new_meta$subclass_label

  new_pse_data <- new_sample_data


  result_list <- list(all_marker_list, new_meta, new_pse_data)
  return(result_list)




}
