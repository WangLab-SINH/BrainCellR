#' @title Get cell type name
#' @description  Creating poisson vector.
#' @details Input an integer and return the log density of a poisson distribution with lambda equals the input integer
#' @param input_data data input
#' @param cell_type_file meta input
#' @param all_marker_list de gene input
#' @return A numeric vector of log density.
#' @export
GetCellTypeName <- function(input_data, cell_type_file, all_marker_list){
  # cell_type_file = read.csv(cell_type_file,row.names=1)
  # all_marker_list = readRDS(all_marker_list)
  total_all_marker_list = all_marker_list[[1]][1,]
  for(i in 1:length(all_marker_list)){
    total_all_marker_list = rbind(total_all_marker_list, all_marker_list[[i]])
  }
  total_all_marker_list = total_all_marker_list[-1,]
  print(length(unique(total_all_marker_list$cluster)))
  # result_file_name = "G:/PFC/mouse_M1/de_gene_consensus_1_2_t_subclass.csv"
  # result_file_name1 = "G:/PFC/mouse_M1/de_gene_consensus_1_2_t_subclass_old.csv"
  # result_file_name2 = "G:/PFC/mouse_M1/de_gene_consensus_1_2_t_subclass_raw.csv"
  # result_score_file_name = "G:/PFC/mouse_M1/scoreconsensus_1_2_t_subclass.csv"
  data_meta = cell_type_file

  cl = data_meta$cluster_label
  names(cl) = rownames(data_meta)
  rm_cluster = names(table(cl))[table(cl)<2]
  data_meta = data_meta[!data_meta$cluster_label %in% rm_cluster,]
  input_data = input_data[,rownames(data_meta)]


  mouse_data = mouse_data <- CreateSeuratObject(counts = input_data, min.cells = 0, min.features = 0, project = "example")
  # mouse_data = mouse_data[,rownames(cell_type_file)]
  mouse_data <- AddMetaData(mouse_data, data_meta)


  cl = data_meta$cluster_label
  names(cl) = rownames(mouse_data@meta.data)
  #temp_exp = as.matrix(mouse_data@assays$RNA@data)
  temp_exp = Matrix(mouse_data@assays$RNA@data, sparse=TRUE)
  temp_exp = as_matrix(temp_exp)
  medianExpr = do.call("cbind", tapply(names(cl), cl, function(x) rowMedians(temp_exp[,x], na.rm=T)))
  rm(temp_exp)
  gc()
  rownames(medianExpr) = rownames(mouse_data@assays$RNA@data)

  #first
  Idents(mouse_data) <- mouse_data$subclass_label
  num = 1
  binary_data_frame = data.frame(0,0,0,0,0,0,0,0,0)
  binary_data_frame_raw = data.frame(0,0,0,0,0,0,0,0,0)
  binary_data_frame_pvalue_fc = data.frame(0,0,0,0,0,0,0,0,0)
  for(i in unique(names(all_marker_list))){
    print(num)
    temp_mouse_data = subset(mouse_data, idents = i)
    Idents(temp_mouse_data) <- temp_mouse_data$cluster_label
    temp_mouse_data_exp = as.matrix(temp_mouse_data@assays$RNA@data)
    temp_cluster = as.character(unique(temp_mouse_data@meta.data$cluster_label))

    for(j in temp_cluster){
      current_temp_cluster = data.frame(0,0,0,0,0,0,0,0,0)
      current_temp_cluster_raw = data.frame(0,0,0,0,0,0,0,0,0)
      current_temp_cluster_pvalue_fc = data.frame(0,0,0,0,0,0,0,0,0)
      temp_score = 0
      if(nrow(all_marker_list[[num]]) > 0){
        if(nrow(all_marker_list[[num]][all_marker_list[[num]]$cluster==j,])>0){
          temp_current_marker_list_temp = all_marker_list[[num]][all_marker_list[[num]]$cluster==j,]
          temp_current_marker_list_temp = temp_current_marker_list_temp[temp_current_marker_list_temp$avg_log2FC>1,]
          temp_current_marker_list_temp = temp_current_marker_list_temp[temp_current_marker_list_temp$pct.1>0.5,]
          if(nrow(temp_current_marker_list_temp) > 0){
            for(l in 1:nrow(temp_current_marker_list_temp)){
              g = temp_current_marker_list_temp[l,"gene"]
              temp_medianExpr = medianExpr[,temp_cluster]
              # median_sec_median_score = max(temp_medianExpr[g,colnames(temp_medianExpr)!=j]) / temp_medianExpr[g,j]
              # if(median_sec_median_score == Inf){
              #   median_sec_median_score = 1
              # }
              # median_sec_median_score = 1 - median_sec_median_score

              temp_score = 0
              other_temp_cluster = temp_cluster[temp_cluster!=j]
              temp_median = median(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==j,])])
              temp_percentage = sum(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==j,])]!=0) / length(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==j,])])
              all_temp_other_percentage = c()
              for(k in other_temp_cluster){
                temp_other_median = median(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==k,])])
                temp_other_percentage = sum(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==k,])]==0)
                all_temp_other_percentage = c(all_temp_other_percentage, temp_other_percentage)
                temp222 = (1 - temp_other_median / temp_median) * temp_other_percentage/length(temp_mouse_data_exp[g,rownames(temp_mouse_data@meta.data[temp_mouse_data@meta.data$cluster_label==k,])])
                if(!is.na(temp222)){
                  if(temp222 < 0){
                    temp222 = 0
                  }
                }
                if(!is.nan(temp222)){
                  if(temp222 < 0){
                    temp222 = 0
                  }
                }

                temp_score = temp_score + temp222
              }
              temp_score = temp_score / (length(temp_cluster) - 1)
              current_temp_cluster = rbind(current_temp_cluster, c(j, temp_score,g,as.character(temp_current_marker_list_temp[l,1:6])))
              current_temp_cluster_raw = rbind(current_temp_cluster_raw, c(j, temp_score,g,as.character(temp_current_marker_list_temp[l,1:6])))
              current_temp_cluster_pvalue_fc = rbind(current_temp_cluster_pvalue_fc, c(j, temp_score,g,as.character(temp_current_marker_list_temp[l,1:6])))
            }
          }

        }
      }

      current_temp_cluster = current_temp_cluster[-1,]
      current_temp_cluster$X0.1 = as.numeric(current_temp_cluster$X0.1)
      current_temp_cluster = current_temp_cluster[!is.nan(current_temp_cluster$X0.1),]
      current_temp_cluster = current_temp_cluster[!is.infinite(current_temp_cluster$X0.1),]
      current_temp_cluster = current_temp_cluster[order(current_temp_cluster$X0.1, decreasing = T),]
      if(nrow(current_temp_cluster)>0){
        current_temp_cluster = current_temp_cluster[1:min(20, nrow(current_temp_cluster)),]
      }else{
        current_temp_cluster1 = data.frame(j,NA,NA,NA,NA,NA,NA,NA,NA)
        colnames(current_temp_cluster1) = colnames(current_temp_cluster)
        current_temp_cluster = current_temp_cluster1
      }

      current_temp_cluster_raw = current_temp_cluster_raw[-1,]
      current_temp_cluster_raw$X0.1 = as.numeric(current_temp_cluster_raw$X0.1)
      current_temp_cluster_raw = current_temp_cluster_raw[!is.nan(current_temp_cluster_raw$X0.1),]
      current_temp_cluster_raw = current_temp_cluster_raw[!is.infinite(current_temp_cluster_raw$X0.1),]
      if(nrow(current_temp_cluster_raw)>0){
        current_temp_cluster_raw = current_temp_cluster_raw[1:min(20, nrow(current_temp_cluster_raw)),]
      }else{
        current_temp_cluster1 = data.frame(j,NA,NA,NA,NA,NA,NA,NA,NA)
        colnames(current_temp_cluster1) = colnames(current_temp_cluster_raw)
        current_temp_cluster_raw = current_temp_cluster1
      }



      current_temp_cluster_pvalue_fc = current_temp_cluster_pvalue_fc[-1,]
      current_temp_cluster_pvalue_fc$X0.1 = as.numeric(current_temp_cluster_pvalue_fc$X0.1)
      current_temp_cluster_pvalue_fc = current_temp_cluster_pvalue_fc[!is.nan(current_temp_cluster_pvalue_fc$X0.1),]
      current_temp_cluster_pvalue_fc = current_temp_cluster_pvalue_fc[!is.infinite(current_temp_cluster_pvalue_fc$X0.1),]
      current_temp_cluster_pvalue_fc$pvalue_index = order(current_temp_cluster_pvalue_fc$X0.3)
      current_temp_cluster_pvalue_fc$fc_index = order(current_temp_cluster_pvalue_fc$X0.4, decreasing = T)
      current_temp_cluster_pvalue_fc$new_index = current_temp_cluster_pvalue_fc$pvalue_index + current_temp_cluster_pvalue_fc$fc_index
      current_temp_cluster_pvalue_fc = current_temp_cluster_pvalue_fc[order(current_temp_cluster_pvalue_fc$new_index),]
      current_temp_cluster_pvalue_fc = current_temp_cluster_pvalue_fc[,-c(10,11,12)]
      if(nrow(current_temp_cluster_pvalue_fc)>0){
        current_temp_cluster_pvalue_fc = current_temp_cluster_pvalue_fc[1:min(20, nrow(current_temp_cluster_pvalue_fc)),]
      }else{
        current_temp_cluster1 = data.frame(j,NA,NA,NA,NA,NA,NA,NA,NA)
        colnames(current_temp_cluster1) = colnames(current_temp_cluster_pvalue_fc)
        current_temp_cluster_pvalue_fc = current_temp_cluster1
      }



      binary_data_frame = rbind(binary_data_frame, current_temp_cluster)
      binary_data_frame_raw = rbind(binary_data_frame_raw, current_temp_cluster_raw)
      binary_data_frame = rbind(binary_data_frame, current_temp_cluster)
      binary_data_frame_pvalue_fc = rbind(binary_data_frame_pvalue_fc, current_temp_cluster_pvalue_fc)
    }

    num = num + 1
  }

  binary_data_frame$X0.1 = as.numeric(binary_data_frame$X0.1)
  binary_data_frame = binary_data_frame[!is.nan(binary_data_frame$X0.1),]
  binary_data_frame = binary_data_frame[!is.infinite(binary_data_frame$X0.1),]

  binary_data_frame = binary_data_frame[-1,]

  binary_data_frame_old = binary_data_frame

  binary_data_frame_new = binary_data_frame[1,]
  for(k in 1:nrow(binary_data_frame)){
    current_gene = binary_data_frame[k,3]
    current_cell_type = binary_data_frame[k,1]
    if(!is.na(current_gene)){
      if(medianExpr[current_gene, current_cell_type] >= (sort(as.numeric(medianExpr[current_gene, ]), decreasing = T)[ceiling(length(medianExpr[current_gene, ])/ 10)])){
        binary_data_frame_new <- rbind(binary_data_frame_new, binary_data_frame[k,])
      }
    }

  }

  binary_data_frame_new = binary_data_frame_new[-1,]
  binary_data_frame = binary_data_frame_new

  binary_gene = data.frame(0,0,0,0,0)
  score_list = c()
  for(i in unique(binary_data_frame$X0)){
    temp = binary_data_frame[binary_data_frame$X0==i,]
    binary_gene = rbind(binary_gene, temp$X0.2[1:5])
    score_list = c(score_list, temp$X0.1[1])
  }
  binary_gene = binary_gene[-1,]
  rownames(binary_gene) = unique(binary_data_frame$X0)

  data_meta_cluster_level = unique(data.frame(data_meta$class,data_meta$subclass_label,data_meta$cluster_label))
  rownames(data_meta_cluster_level) = data_meta_cluster_level$data_meta.cluster_label
  colnames(data_meta_cluster_level) = c("class","subclass","group")
  data_meta_final = data.frame(0,0,0,0,0,0,0,0)
  for(i in 1:nrow(data_meta_cluster_level)){
    if(rownames(data_meta_cluster_level)[i] %in% rownames(binary_gene)){
      temp = c(as.character(data_meta_cluster_level[i,]), as.character(binary_gene[rownames(data_meta_cluster_level)[i],]))
    }else{
      temp = c(as.character(data_meta_cluster_level[i,]), NA,NA,NA,NA,NA)
    }
    data_meta_final = rbind(data_meta_final, temp)
  }
  data_meta_final = data_meta_final[-1,]

  colnames(data_meta_final) = c("class","subclass","cluster","gene1","gene2","gene3","gene4","gene5")
  data_meta_final$cluster_new = data_meta_final$cluster
  for(k in 1:nrow(data_meta_final)){
    temp_name = paste0(paste0(paste0(paste0(data_meta_final[k,"class"],"_"),data_meta_final[k,"subclass"]),"_"),data_meta_final[k,"gene1"])
    data_meta_final[k,"cluster_new"] = temp_name
  }
  # write.csv(data_meta_final, result_file_name)
  return(data_meta_final)





}

#' @title Change to matrix
#' @description  Creating poisson vector.
#' @details Input an integer and return the log density of a poisson distribution with lambda equals the input integer
#' @param mat data input
#' @return A numeric vector of log density.
#' @export
as_matrix <- function(mat){

  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])

  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x

  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }

  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}
