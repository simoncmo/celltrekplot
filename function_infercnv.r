### Helper FUNCTION for inferCNV
GetInferCNVSubcluster = function(infercnv_obj){
    imap(infercnv_obj@tumor_subclusters$subclusters, function(celltype_result, celltype){
        imap(celltype_result, function(obj, group){
            data.frame(Cell = obj %>% names, infercnv_subcluster = group, infercnv_celltype = celltype)
        }) %>% bind_rows()
    }) %>% bind_rows() %>% 
    column_to_rownames('Cell')
}

setRownames = function(df, rownames){
    rownames(df) = rownames
    return(df)
}

# Bind col when there are different length
# Might need to check if col exisit already 
bind_cols_dfs = function(df1, df2){
  # different length
  if(nrow(df1) != nrow(df2)){
      df_long  = if(nrow(df1) > nrow(df2)) df1 else df2
      df_short = if(nrow(df1) > nrow(df2)) df2 else df1
      missing_rows = setdiff(rownames(df_long), rownames(df_short))
      #df_short_cols = names(df_short)
      df_short_final = bind_rows(df_short, data.frame(row.names = missing_rows))
  }else{
      df_long = df1
      df_short_final = df2
  }
  if(intersect(colnames(df_long), colnames(df_short_final)) %>% length !=0) stop('Duplicated col found! Not merging. Please double check inputs')
  bind_cols(df_long, df_short_final[rownames(df_long), ,drop=F])
}
