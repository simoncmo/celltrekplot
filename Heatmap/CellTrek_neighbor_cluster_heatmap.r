## HEATMAP
# 6/23/2022
#### This function take plot obj and sc, then plot heatmap of Neighborhood most variable genes
#### DEMO #####
# obj_use: Plot object from "PrepareSpatialPlotObject"
# sn_use: the sc/snRNA object
# cell_type_plt = 'immune'
# GetExpressionData(obj_use, sn_use, cell_type = cell_type_plt) %>% 
#     GetTopGeneMatrix %>% 
#     MakeClusterHeatmap


# 1. expression subgroup
GetExpressionData = function(plot_obj, sc_obj, assay = 'SCT', cell_type){
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

## 2 Select genes to plot
## Idea A: Top n most variable genes
GetTopGeneMatrix = function(expr_matrix, n_gene = 100){
    top_var_genes = apply(expr_matrix, 2, var, na.rm = T) %>% 
        .[.>0] %>% # filter out var 0 genes
        sort(., decreasing = T) %>% 
        .[1:n_gene] %>% names # Top n genes
    expr_matrix = expr_matrix[,top_var_genes] # column is ordered by variance of gene
}

## 3. Plot
MakeClusterHeatmap = function(expression_matrix, gene_to_label){
    ## Gene to mark
    if(missing(gene_to_label)){
        message('By default label top most variable genes')
        gene_to_label = colnames(expression_matrix)[1:20] # Matrix ordered. Top 20 most variable genes
    }
    
    # Make Gene label marks
    gene_mark  = intersect(names(expression_matrix), gene_to_label)
    gene_mark_pos = match(gene_mark, names(exp_top_df))
    right_anno = rowAnnotation(foo = anno_mark(at = gene_mark_pos, labels = gene_mark))

    ## PLOT
    Heatmap(exp_top_df %>% scale %>% t,  right_annotation = right_anno,
            name = 'Scaled\nExpression',
            show_column_names = F,
            show_row_names = F,
                show_row_dend = F,
            column_split = sn_celltrek_plot_meta$Cell_cluster,
                cluster_column_slices = F
            )
}

