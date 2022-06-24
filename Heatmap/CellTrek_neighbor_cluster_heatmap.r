## HEATMAP
# 6/23/2022
#### This function take plot obj and sc, then plot heatmap of Neighborhood most variable genes
#### DEMO #####
# obj_use: Plot object from "PrepareSpatialPlotObject"
# sn_use: the sc/snRNA object
# MakeClusterHeatmap(obj_use, sn_use, assay = 'SCT', n_gene = 100)

## [Internal] 1. expression subgroup 
GetExpressionData = function(plot_obj, sc_obj, assay = 'SCT'){
    # Use MakeNames on scRNA to match that in cellTrek
    sc_obj = RenameCells(sc_obj, new.names = colnames(sc_obj) %>% make.names())
    
    # Shared ids 
    shared_ids = intersect(colnames(sc_obj), plot_obj$Neighbor_cluster_table$Origin)
    
    ## Subset scRNA
    sc_subset  = subset(sc_obj, cells = shared_ids)

    ## Get Variable genes
    sc_subset = sc_subset %>% FindVariableFeatures()

    ## Get expression data
    DefaultAssay(sc_subset) = assay
    exp_df = FetchData(sc_subset, vars = VariableFeatures(sc_subset))

    ## remove zero genes
    exp_df = exp_df[exp_df %>% colSums() != 0]
    
    # Subset/reorder both table
    exp_df_plot = exp_df[shared_ids, ]
    #sn_celltrek_plot_meta = plot_obj$Neighbor_cluster_table %>% column_to_rownames('Origin') %>% .[shared_ids, ]
    
    exp_df_plot
}

##  [Internal]  2 Select genes to plot - Top n most variable genes
GetTopGeneMatrix = function(expr_matrix, n_gene = 100){
    top_var_genes = apply(expr_matrix, 2, var, na.rm = T) %>% 
        .[.>0] %>% # filter out var 0 genes
        sort(., decreasing = T) %>% 
        .[1:n_gene] %>% names # Top n genes
    expr_matrix[,top_var_genes] # column is ordered by variance of gene
}

## 3. Plot
MakeClusterHeatmap = function(plot_obj, sc_obj, assay = 'SCT', n_gene = 100, gene_to_label){
    ## Get Expression mtx
    expression_matrix = GetExpressionData(plot_obj, sc_obj, assay)
    expression_matrix = GetTopGeneMatrix(expression_matrix, n_gene) 
    
    ### PLOT
    # 1. Gene to mark
    if(missing(gene_to_label)){
        message('By default label top most variable genes')
        gene_to_label = colnames(expression_matrix)[1:20] # Matrix ordered. Top 20 most variable genes
    }
    
    # 2. Column split
    plot_meta = plot_obj$Neighbor_cluster_table %>% 
        column_to_rownames('Origin') %>% .[rownames(expression_matrix), ]
    
    # 3. Make Gene label marks
    gene_mark  = intersect(names(expression_matrix), gene_to_label)
    gene_mark_pos = match(gene_mark, names(expression_matrix))
    right_anno = rowAnnotation(foo = anno_mark(at = gene_mark_pos, labels = gene_mark))

    ## 4. PLOT
    Heatmap(expression_matrix %>% scale %>% t,  right_annotation = right_anno,
            name = 'Scaled\nExpression',
            show_column_names = F,
            show_row_names = F,
            show_row_dend = F,
            column_split = plot_meta$Cell_cluster,
                cluster_column_slices = F
            )
}

