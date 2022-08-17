## HEATMAP
library(ComplexHeatmap)
# 6/23/2022
#### This function take plot obj and sc, then plot heatmap of Neighborhood most variable genes
#### DEMO #####
# obj_use: Plot object from "PrepareSpatialPlotObject"
# sn_use: the sc/snRNA object
# MakeClusterHeatmap(obj_use, sn_use, assay = 'SCT', n_gene = 100)

## [Internal] 1. expression subgroup 
## This ONLY EXTRACT expression of Target cells
GetExpressionData = function(plot_obj, sc_obj, assay = 'SCT', cell_mode = c('target_cells','all'), gene_mode = c('Variable','DEG'), gene_use =NULL){
    # Use MakeNames on scRNA to match that in cellTrek
    sc_obj = RenameCells(sc_obj, new.names = colnames(sc_obj) %>% make.names())
    
    # Choose Cell
    cell_mode = match.arg(cell_mode)
    if(cell_mode == 'target_cells'){
        # Target CElls only
        shared_ids = intersect(colnames(sc_obj), plot_obj$Neighbor_cluster_table$Origin)
    }else{
        # All Cells
        shared_ids = intersect(colnames(sc_obj), plot_obj$meta_selected$id_new)
    }

    # Mode
    # If supply gene_use : ignore gene_mode by pass this whole part
    if(is.null(gene_use)){
        # Choose Genes
        gene_mode = match.arg(gene_mode)
        if(gene_mode =='Variable'){
            message('Variable MODE')
            ## Subset scRNA
            sc_obj  = subset(sc_obj, cells = shared_ids)

            ## Get Variable genes
            sc_obj = sc_obj %>% FindVariableFeatures()
            gene_use = VariableFeatures(sc_obj)
        }else{
            ## DEG
            if(is.null(plot_obj[['Markers']])){
                stop('Please run GetClusterMarkers on the plot object first!')
            }
            message('DEG MODE')
            ## Get DEGs
            gene_use = plot_obj$Markers$Cluster_DEG_summary$Gene
        }
    }
    
    ## Get expression data
    DefaultAssay(sc_obj) = assay
    exp_df = FetchData(sc_obj, vars = gene_use, cells = shared_ids)

    ## remove zero genes
    exp_df = exp_df[,exp_df %>% colSums() != 0]
    
    # Subset/reorder both table
    exp_df_plot = exp_df[shared_ids, ]
    #sn_celltrek_plot_meta = plot_obj$Neighbor_cluster_table %>% column_to_rownames('Origin') %>% .[shared_ids, ]
    
    exp_df_plot
}

##  [Internal]  2 Select genes to plot - Top n most variable genes
GetTopGeneMatrix = function(expr_matrix, n_gene = 100){
    n_gene = min(ncol(expr_matrix), n_gene)
    top_var_genes = apply(expr_matrix, 2, var, na.rm = T) %>% 
        .[.>0] %>% # filter out var 0 genes
        sort(., decreasing = T) %>% 
        .[1:n_gene] %>% names # Top n genes
    expr_matrix[,top_var_genes] # column is ordered by variance of gene
}

## 3. Plot
MakeClusterHeatmap = function(plot_obj, sc_obj, assay = 'SCT', limit_gene_show = T, n_gene_show = 100, gene_to_label,
                              # for GetExpressionData
                             cell_mode = c('target_cells','all'), gene_mode = c('Variable','DEG'), gene_list = NULL,
                              # Mark
                              n_mark_marker = 5 # DEG to mark per group
                             ){
    ## parameter 
    cell_mode = match.arg(cell_mode)
    gene_mode = match.arg(gene_mode)
    
    ## Get Expression mtx
    gene_to_plot = unique(unlist(gene_list)) 
    expression_matrix = GetExpressionData(plot_obj, sc_obj, assay, cell_mode = cell_mode, gene_mode = gene_mode, gene_use =gene_to_plot)
    
    ## limit # gene to plot
    expression_matrix = if(!limit_gene_show) expression_matrix else GetTopGeneMatrix(expression_matrix, n_gene_show) 
    
    ### PLOT
    # 1. Gene to mark
    if(missing(gene_to_label)){
        if(gene_mode == 'Variable'){
            message('Varialbe mode: By default label top most variable genes')
            gene_to_label = colnames(expression_matrix)[1:20] # Matrix ordered. Top 20 most variable genes
        }else{
            message('DEG mode: By default label top n gene in each group')
            gene_to_label = GetTopDegPerGroup(plot_obj, n_mark_marker)
        }
        
    }
    
    # 2. Column split (cell)
    cell_cluster_df = plot_obj$Neighbor_cluster_table %>% select(Origin, Cell_cluster) %>% distinct %>% 
        column_to_rownames('Origin') %>% .[rownames(expression_matrix), ,drop = F]
    
    # 3. Make Gene label marks
    gene_mark  = intersect(names(expression_matrix), gene_to_label)
    gene_mark_pos = match(gene_mark, names(expression_matrix))
    right_anno = rowAnnotation(foo = anno_mark(at = gene_mark_pos, labels = gene_mark))
    
    # 4. Row split (gene)
    if(!is.null(gene_list)){ # Custom gene list- group gene that way
        message('User supplied Gene list. Use item name to split row')
        gene_groups = gene_list %>% imap(~rep(.y, length(.x)) %>% setNames(.x)) %>% reduce(c) # Gene as key, group as value
        # Note this way.. gene in multiple list, only 'First' occurance will match
        row_groups = gene_groups[colnames(expression_matrix)] # keep only gene in the matrix
    }else if(gene_mode == 'DEG'){
        # Get Gene Groups
        row_groups = plot_obj$Markers$Cluster_DEG_summary$Groups
    }else if(gene_mode == 'Variable'){ # For Variable 
        # Use Kmeans , k = # of 'groups'
        message('Variable Mode use kmeans for each row')
        tmp_k = expression_matrix %>% scale %>% t %>% kmeans(x = ., centers = length(unique(cell_cluster_df$Cell_cluster)))
        row_groups = tmp_k$cluster
    }
 
    ## 4. PLOT
    message('Note: Mark position will be wrong in interactive interphse. Need output to pdf')
    Heatmap(expression_matrix %>% scale %>% t, 
            name = 'Scaled\nExpression',
            # Col
            show_column_dend = F,
            show_column_names = F,
            # Row
            show_row_names = F,
            show_row_dend = F,
            right_annotation = right_anno,
            # Col Split
            column_split = cell_cluster_df$Cell_cluster,
            cluster_column_slices = F,
            # Row Split
            row_split = row_groups,
            row_title_rot = 0,
            # Top
            column_title_gp = gpar(fill = plot_obj$Palettes$palette_cluster, lwd=0, font = 2)
            )
}


## MARKERs
## [API] This function Get Markers from Neighborhood based target cell clusters
GetClusterMarkers = function(scRNA, plot_obj, assay = 'SCT', p_cutoff = 0.05, fc_cutoff = 0.2){
    # sn obj with matching name 
    scRNA = RenameSeuratCells(scRNA)
    
    # Get clutser
    Neighbor_cluster_df = plot_obj$Neighbor_cluster_table %>% 
            filter(Origin %in% colnames(scRNA)) # All Target cells in snRNA only
    all_clusters        = Neighbor_cluster_df$Cell_cluster %>% unique
    all_target_ids = Neighbor_cluster_df$Origin 
    
    # DEG
    deg_df = map(all_clusters, function(cluster){
        message(str_glue('running {cluster} ..'))
        cell1_ids = Neighbor_cluster_df %>% subset(Cell_cluster %in% cluster) %>% pull(Origin) # Cell ids
        cell2_ids = setdiff(all_target_ids, cell1_ids) # All other cells 
        FindMarkers(scRNA@assays[[assay]], 
                    cells.1 = cell1_ids, cells.2 = cell2_ids) %>% 
            mutate(Cluster = cluster,
                   cell1_count = length(cell1_ids),
                   cell2_count = length(cell2_ids)
                  ) %>% 
            rownames_to_column('Gene')
    }) %>% bind_rows()
    
    ## Filter
    deg_sig = deg_df %>% filter(p_val_adj < p_cutoff, abs(avg_log2FC) > fc_cutoff) 
    deg_sig_genes = deg_df %>% arrange(desc(avg_log2FC)) %>% pull(Gene)
    
    ## Summary table
    # Get Row label
    deg_group_df = deg_sig %>% mutate(Group = ifelse(avg_log2FC>0, str_glue('High in {Cluster}'), str_glue('Low in {Cluster}'))) %>% 
        group_by(Gene) %>% summarise(Groups = toString(Group), avg_log2FC = mean(avg_log2FC))
    
    ## Save to plot object
    plot_obj$Markers = list(Cluster_DEG_table = deg_df, Cluster_DEG_table_filtered = deg_sig,
                            Cluster_DEG = deg_sig_genes, Parameters = list(p_cutoff = p_cutoff, fc_cutoff=fc_cutoff),
                            Cluster_DEG_summary = deg_group_df
                           )
    return(plot_obj)
}

## [Internal]
# Get top DEGs per group 
# This is used in select gene for anno_mark for Heatamp
GetTopDegPerGroup = function(plot_obj, top_n = 5){
    plot_obj$Markers$Cluster_DEG_summary %>% 
        group_by(Groups) %>% 
        slice_max(order_by = abs(avg_log2FC), n = top_n) %>% 
        pull(Gene)

}