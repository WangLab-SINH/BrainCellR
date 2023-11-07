#' @title Classify class and subclass
#' @description  Classify class and subclass.
#' @param new_data Data input, single cell expression data.
#' @param group_meta meta input
#' @param ref_data Reference expression data for classication.
#' @param ref_meta Reference metadata, must include column name with class_label and subclass_label.
#' @param method Method to run classification, can be choosen from Seurat or SingleR. Default is Seurat.
#' @return A dataframe with annotation of cell type major class and subclass.
#' @export
#' @import Seurat
#' @import SingleR
RunSubclassClassify <- function(new_data, group_meta, ref_data, ref_meta,method="Seurat")
{

  if(method=="Seurat"){
    common_gene = intersect(rownames(ref_data), rownames(new_data))
    if(length(common_gene)<1){
      message("There is no common gene between reference data and query data.")
      return(NULL)
    }
    new_data = new_data[common_gene,]
    reference_count <- ref_data[common_gene,]
    reference_meta <- ref_meta
    reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "example")
    reference_seurat <- AddMetaData(reference_seurat, reference_meta)

    query_count <- new_data
    query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "example")

    reference_seurat <- NormalizeData(object = reference_seurat)
    reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = 'vst', nfeatures = 2000)
    reference_seurat <- ScaleData(reference_seurat)

    query_seurat <- NormalizeData(object = query_seurat)
    query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = 'vst', nfeatures = 2000)
    query_seurat <- ScaleData(query_seurat)

    k.score = min(30, floor(ncol(query_seurat)/2))
    k.anchor = min(5, floor(ncol(query_seurat)/2))
    k.weight = min(10, floor(ncol(query_seurat)/2))
    error_flag = F
    result = tryCatch({
      sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = query_seurat,
                                         dims = 1:30,k.score = k.score, k.anchor = k.anchor)
    }, error = function(e){
      error_flag = T
    })
    if(!exists("sim.anchors")){
      sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = query_seurat,
                                         dims = 1:30,k.score = k.score, k.anchor = k.anchor, project.query = T)
    }

    ##replace Group with the actual column name from meta
    k.weight = min(floor(nrow(sim.anchors@anchors)/2),k.weight)
    k.weight = max(k.weight,3)
    predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$class_label,
                                dims = 1:30,k.weight=k.weight)
    query_seurat <- AddMetaData(object = query_seurat, metadata = predictions)

    predict_meta = data.frame(rownames(query_seurat@meta.data), query_seurat@meta.data$predicted.id)

    rownames(predict_meta) = predict_meta$rownames.query_seurat.meta.data.
    colnames(predict_meta) = c("group","class_label")
    for(i in unique(group_meta$group)){
      temp = predict_meta[rownames(group_meta[group_meta$group==i,]),]
      type = names(table(temp$class_label)[table(temp$class_label)==max(table(temp$class_label))])
      if(max(table(temp$class_label)) <= (nrow(temp)/2)){
        predict_meta[rownames(group_meta[group_meta$group==i,]),]$class_label = NA
      }else{
        predict_meta[rownames(group_meta[group_meta$group==i,]),]$class_label = type
      }

    }
    predict_meta = predict_meta[!is.na(predict_meta$class_label),]
    new_data = new_data[,rownames(predict_meta)]

    all_label = unique(predict_meta$class_label)
    all_meta = predict_meta[1,]
    for(i in all_label){
      reference_count <- ref_data[common_gene,][,ref_meta$class_label==i]
      reference_meta <- ref_meta[ref_meta$class_label==i,]
      reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "example")
      reference_seurat <- AddMetaData(reference_seurat, reference_meta)

      #load query
      query_count <- new_data[,predict_meta$class_label==i]

      query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "example")
      query_seurat <- NormalizeData(object = query_seurat)

      #standard pipeline
      reference_seurat <- NormalizeData(object = reference_seurat)
      reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = 'vst', nfeatures = 2000)
      reference_seurat <- ScaleData(reference_seurat)

      query_seurat <- NormalizeData(object = query_seurat)
      query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = 'vst', nfeatures = 2000)
      query_seurat <- ScaleData(query_seurat)

      ##prediction###
      k.score = min(30, floor(ncol(query_seurat)/2))
      k.anchor = min(5, floor(ncol(query_seurat)/2))
      k.weight = min(10, floor(ncol(query_seurat)/2))
      error_flag = F
      result = tryCatch({
        sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = query_seurat,
                                           dims = 1:30,k.score = k.score, k.anchor = k.anchor)
      }, error = function(e){
        error_flag = T
      })
      if(!exists("sim.anchors")){
        sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = query_seurat,
                                           dims = 1:30,k.score = k.score, k.anchor = k.anchor, project.query = T)
      }
      ##replace Group with the actual column name from meta
      k.weight = min(floor(nrow(sim.anchors@anchors)/2),k.weight)
      k.weight = max(k.weight,3)
      predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$subclass_label,
                                  dims = 1:30,k.weight=k.weight)
      query_seurat <- AddMetaData(object = query_seurat, metadata = predictions)

      predict_meta1 = data.frame(rownames(query_seurat@meta.data), query_seurat@meta.data$predicted.id)

      rownames(predict_meta1) = predict_meta1$rownames.query_seurat.meta.data.
      colnames(predict_meta1) = c("group","class_label")
      all_meta <- rbind(all_meta, predict_meta1)
    }

    all_meta = all_meta[-1,]
    rownames(all_meta) = all_meta$group
    colnames(all_meta) = c("cell","subclass")
    all_meta_temp = all_meta
    for(i in unique(group_meta$group)){
      temp = all_meta[rownames(group_meta[group_meta$group==i,]),]
      type = names(table(temp$subclass)[table(temp$subclass)==max(table(temp$subclass))])
      if(max(table(temp$subclass)) <= (nrow(temp)*0.8)){
        all_meta[rownames(group_meta[group_meta$group==i,]),]$subclass = NA
      }else{
        all_meta[rownames(group_meta[group_meta$group==i,]),]$subclass = type
      }

    }
    all_meta = all_meta[!is.na(all_meta$subclass),]
  }else if(method=="SingleR"){
    common_gene = intersect(rownames(ref_data), rownames(new_data))
    if(length(common_gene)<1){
      message("There is no common gene between reference data and query data.")
      return(NULL)
    }
    new_data = new_data[common_gene,]
    reference_count <- ref_data[common_gene,]
    reference_meta <- ref_meta
    query_count <- new_data

    predictions <- SingleR(test=query_count, assay.type.test=1,
                           ref=reference_count, labels=reference_meta$subclass_label)




    predict_meta = data.frame(rownames(predictions), predictions$labels)

    rownames(predict_meta) = predict_meta[,1]
    colnames(predict_meta) = c("group","class_label")
    for(i in unique(group_meta$group)){
      temp = predict_meta[rownames(group_meta[group_meta$group==i,]),]
      type = names(table(temp$class_label)[table(temp$class_label)==max(table(temp$class_label))])
      if(max(table(temp$class_label)) <= (nrow(temp)/2)){
        predict_meta[rownames(group_meta[group_meta$group==i,]),]$class_label = NA
      }else{
        predict_meta[rownames(group_meta[group_meta$group==i,]),]$class_label = type
      }

    }
    predict_meta = predict_meta[!is.na(predict_meta$class_label),]
    new_data = new_data[,rownames(predict_meta)]

    all_label = unique(predict_meta$class_label)
    all_meta = predict_meta[1,]
    for(i in all_label){
      reference_count <- ref_data[common_gene,][,ref_meta$class_label==i]
      reference_meta <- ref_meta[ref_meta$class_label==i,]

      #load query
      query_count <- new_data[,predict_meta$class_label==i]
      predictions <- SingleR(test=query_count, assay.type.test=1,
                             ref=reference_count, labels=reference_meta$subclass_label)


      predict_meta1 = data.frame(rownames(predictions), predictions$labels)

      rownames(predict_meta1) = predict_meta1[,1]
      colnames(predict_meta1) = c("group","class_label")
      all_meta <- rbind(all_meta, predict_meta1)
    }

    all_meta = all_meta[-1,]
    rownames(all_meta) = all_meta$group
    colnames(all_meta) = c("cell","subclass")
    all_meta_temp = all_meta
    for(i in unique(group_meta$group)){
      temp = all_meta[rownames(group_meta[group_meta$group==i,]),]
      type = names(table(temp$subclass)[table(temp$subclass)==max(table(temp$subclass))])
      if(max(table(temp$subclass)) <= (nrow(temp)*0.8)){
        all_meta[rownames(group_meta[group_meta$group==i,]),]$subclass = NA
      }else{
        all_meta[rownames(group_meta[group_meta$group==i,]),]$subclass = type
      }

    }
    all_meta = all_meta[!is.na(all_meta$subclass),]
  }


  data_meta = data.frame(all_meta$cell, predict_meta[rownames(all_meta),]$class_label, all_meta$subclass,group_meta[rownames(all_meta),]$group)
  rownames(data_meta) = data_meta$all_meta.cell
  colnames(data_meta) = c("cell","class_label","subclass_label","cluster_label")
  # write.csv(data_meta, paste0(file_path_root, "new_meta1.csv"))
  return(data_meta)
}
