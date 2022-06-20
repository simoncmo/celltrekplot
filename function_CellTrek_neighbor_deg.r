## 5/23/2022
## DEG function 
## This function Select Neighborhood cells from CellTrek obj
## and find cell type specific DEG vs non-neighbor cells
## Requied :  
# - Celltrek obj
# - Graph_obj
# - group.by: Column name of Celltype in meta.data
# - ident: Celltype of Interests
# Optional
# - n_order : hop/order of to consider as neighbor. Default 2. No recommand over 3
# - cell_count_cutoff : minumun cell count in neighbor cells to run DEG analysis. default 50
# - Assay : The assay used to run DEG. Default is 'R A'
# - return_as : default retun both DEG And META Table for downstream analysis. Set to 'markers_only' if only need the deg 
FindNeighborMarkers = function(celltrek_obj, graph_obj, group.by, ident, n_order = 2, cell_count_cutoff = 50, p_cutoff = 0.05, Assay = 'RNA', return_as= c('both','markers_only')){
    neighbor_df = GetNeighborByIdent(obj = celltrek_obj, 
            obj_graph = graph_obj, 
            group.by =  group.by,
            ident = ident,
            n_order = n_order,
            exclude_self=TRUE) # exclude_self : remove query point themself

    # Define Neighbor cells
    meta_df = celltrek_obj@meta.data %>% mutate(Group = ifelse(.data[[group.by]] %in% ident, 'Target_cell', 'Not_neighbor')) 
    meta_df[neighbor_df, 'Group'] = 'Neighbor'
    #celltrek_obj@meta.data = meta_df
    
    # Select Cell Types to do DEG
    cell_count_cutoff = 50
    cell_test_df = meta_df %>% 
        count(Group, .data[[group.by]]) %>% filter(n > cell_count_cutoff) %>% # filter by minimun cell count to do DEG
        count(.data[[group.by]]) %>% filter(n >1) %>% # Keep only cell type that BOTH neighbor + not_neighbor > cell_count_cutoff %>% 
        pull(.data[[group.by]])
    message(str_glue('Identified {length(cell_test_df)} Cell Types passed cell count cutoff: {cell_count_cutoff}'))
    message(str_glue('Cell Types to be tested: {toString(cell_test_df)}'))

    # Select cells to keep
    cell_keep = meta_df %>% filter(.data[[group.by]] %in% cell_test_df) %>% rownames
    
    # DEGs 
    deg_df = map(cell_test_df, function(celltype){
        cell_neighbor     = meta_df %>% filter( .data[[group.by]] %in% celltype, Group == 'Neighbor') %>% rownames()
        cell_not_neighbor = meta_df %>% filter( .data[[group.by]] %in% celltype, Group == 'Not_neighbor') %>% rownames()
        FindMarkers(celltrek_obj@assays[[Assay]], cells.1 = cell_neighbor, cells.2 = cell_not_neighbor) %>% ## Currently doen't have SCT!!!!
            rownames_to_column('Gene') %>% 
            mutate(Target_celltype = ident, 
                   Neighbor_celltype = celltype, 
                   n_cell_neighbor = length(cell_neighbor),
                   n_cell_not_neighbor = length(cell_not_neighbor)
                  )
    }) %>% bind_rows()
    # Filter
    deg_df = deg_df %>% filter(p_val_adj < p_cutoff)
    message(str_glue('Done! Found {length(unique(deg_df$Gene))} unique Markers at p < {p_cutoff}'))
    
    # return
    return_as = match.arg(return_as)
    result = if(return_as=='both') list(markers = deg_df, meta_table = meta_df) else deg_df
    return(result)
}
    

## 5/23/2022
## DEG function 
## This function Select Neighborhood cells from CellTrek obj
## and find cell type specific DEG vs non-neighbor cells
## Requied :  
# - Celltrek obj
# - Graph_obj
# - group.by: Column name of Celltype in meta.data
# - ident: Celltype of Interests
# Optional
# - n_order : hop/order of to consider as neighbor. Default 2. No recommand over 3
# - cell_count_cutoff : minumun cell count in neighbor cells to run DEG analysis. default 50
# - Assay : The assay used to run DEG. Default is 'R A'
# - return_as : default retun both DEG And META Table for downstream analysis. Set to 'markers_only' if only need the deg 
FindNeighborMarkers = function(celltrek_obj, graph_obj, group.by, ident, n_order = 2, cell_count_cutoff = 50, p_cutoff = 0.05, Assay = 'RNA', return_as= c('all','markers_only')){
    neighbor_df = GetNeighborByIdent(obj = celltrek_obj, 
            obj_graph = graph_obj, 
            group.by =  group.by,
            ident = ident,
            n_order = n_order,
            exclude_self=TRUE) # exclude_self : remove query point themself

    # Define Neighbor cells
    meta_df = celltrek_obj@meta.data %>% mutate(Group = ifelse(.data[[group.by]] %in% ident, 'Target_cell', 'Not_neighbor')) 
    meta_df[neighbor_df, 'Group'] = 'Neighbor'
    #celltrek_obj@meta.data = meta_df
    
    # Select Cell Types to do DEG
    cell_count_cutoff = 50
    cell_test_df = meta_df %>% 
        count(Group, .data[[group.by]]) %>% filter(n > cell_count_cutoff) %>% # filter by minimun cell count to do DEG
        count(.data[[group.by]]) %>% filter(n >1) %>% # Keep only cell type that BOTH neighbor + not_neighbor > cell_count_cutoff %>% 
        pull(.data[[group.by]])
    message(str_glue('Identified {length(cell_test_df)} Cell Types passed cell count cutoff: {cell_count_cutoff}'))
    message(str_glue('Cell Types to be tested: {toString(cell_test_df)}'))

    # Select cells to keep
    cell_keep = meta_df %>% filter(.data[[group.by]] %in% cell_test_df) %>% rownames
    
    # DEGs 
    deg_df = map(cell_test_df, function(celltype){
        cell_neighbor     = meta_df %>% filter( .data[[group.by]] %in% celltype, Group == 'Neighbor') %>% rownames()
        cell_not_neighbor = meta_df %>% filter( .data[[group.by]] %in% celltype, Group == 'Not_neighbor') %>% rownames()
        FindMarkers(celltrek_obj@assays[[Assay]], cells.1 = cell_neighbor, cells.2 = cell_not_neighbor) %>% ## Currently doen't have SCT!!!!
            rownames_to_column('Gene') %>% 
            mutate(Target_celltype = ident, 
                   Neighbor_celltype = celltype, 
                   n_cell_neighbor = length(cell_neighbor),
                   n_cell_not_neighbor = length(cell_not_neighbor)
                  )
    }) %>% bind_rows()
    
    # Filter
    deg_filtered_df = deg_df %>% filter(p_val_adj < p_cutoff)
    message(str_glue('Done! Found {length(unique(deg_filtered_df$Gene))} unique Markers at p < {p_cutoff}'))
    
    # return
    return_as = match.arg(return_as)
    result = if(return_as=='all') list(raw_marker_table = deg_df, marker_table = deg_filtered_df, meta_table = meta_df) else deg_df
    return(result)
}
    
################################################
### DEG Neighbor Expression plot
################################################
## Required:
# - cell column: column with 'cell type info'
## Optional:
# - meta_table: meta_table extra column use to spit plot the not in celltrek_celltrek_obj. 
# Must have same row with celltrek_celltrek_obj and has split.by column

### DEMO
# 1. Plot Split by Cell Type
# SplitSpatialCellTrekPlot(obj, feature = 'Slco1a1', split.by = 'cell_type_after_transfer', pt_size = 0.1)
# 2. Plot DEG result 
# Use MakeNeighborMarkerVlnplot wrapper
## More general version
SplitSpatialCellTrekPlot = function(celltrek_obj, feature, split.by, meta_table, plot_title, ...){
    p = SpatialCellTrekPlot(celltrek_obj, group.by = feature, label = F, ...)

    # Get split by info and add to plot
    if(missing(meta_table)) meta_table = celltrek_obj@meta.data
    p$data[,split.by] = meta_table[rownames(celltrek_obj@meta.data), split.by, drop=F]

    # Split plot
    if(missing(plot_title)) plot_title = str_glue('{feature}: Spatial Expression Split by {split.by}')
    p + facet_wrap(.~.data[[split.by]]) + 
        plot_annotation(title = plot_title, theme = theme(plot.title = element_text(face = 'bold', size = 15, hjust =0.5))) +
        scale_color_gradientn(colours = c('gray85', RColorBrewer::brewer.pal(9, 'YlOrRd')[2:9]), name = 'Expression') + 
        theme(axis.ticks = element_blank(), axis.text  = element_blank())
}

# - marker_result_list : Result list from FindNeighborMarker Function. Must keep return.as = 'both'
# This is wrapper function for Plotting DEG result
MarkerSplitSpatialCellTrekPlot = function(celltrek_obj, marker_result_list, feature, split_cell_col = NULL, ...){
    # Get meta table
    meta_table = marker_result_list$meta_table
    # Change and Order group
    target_cell = marker_result_list$marker_table$Target_celltype %>% unique 
    group_order = c('Target_cell',"Neighbor",  "Not_neighbor")
    meta_table = meta_table %>% mutate(Group = factor(Group, levels = group_order), 
                                       Group = forcats::fct_recode(Group, {{target_cell}} := 'Target_cell')) # Change 'Target cell' to its actually name
    
    # Plot title
    plot_title = str_glue('Neighbor Marker {feature}')
    
    p_st = SplitSpatialCellTrekPlot(celltrek_obj, feature, split.by = 'Group', meta_table, plot_title, ...)
    # Split by cel
    if(!is.null(split_cell_col)){
        # Add Cell type
        p_st$data = bind_cols(p_st$data, FetchData(obj, vars = split_cell_col)[rownames(p_st$data), ,drop=F])
        p_st = p_st + facet_grid(cell_type_Abbr~Group)
    }
    return(p_st)
}



#################################
# Plot function
#################################
# 5/25/2022
# Wrapper function to visualize DEG result
MarkerBubblePlot = function(marker_result_list, palette='', fc_cutoff = 0.5, p_cutoff = 0.05, max.overlaps = 10, expand_x_multi = 0.2){
    palette = ExtendColorVector(palette, target_class = marker_result_list$raw_marker_table$Neighbor_celltype)
    rank_marker_df = marker_result_list$raw_marker_table %>% 
         mutate(Label = ifelse(abs(avg_log2FC) > fc_cutoff & p_val_adj < p_cutoff, Gene, ''))    # Label Gene pass cutoffs

    jitter_pos = position_jitter(seed = 42)    # Set position for jitter points
    rank_marker_df %>% 
        ggplot(aes(x = Neighbor_celltype, y = avg_log2FC, color = Neighbor_celltype)) + 
            geom_jitter(aes(size = -log10(p_val_adj)),alpha = 0.7, position = jitter_pos) + 
            geom_hline(yintercept = c(-fc_cutoff,fc_cutoff), linetype ='dashed') + 
            ggrepel::geom_text_repel(aes(label = Label), max.overlaps = max.overlaps, force = 30, position = jitter_pos) +
            cowplot::theme_cowplot() + 
            scale_color_manual(values = palette)+
            scale_size(range = c(0.1,2)) + 
            scale_x_discrete(expand = expansion(mult = expand_x_multi, add = 0)) +
            labs(x ='', title = 'Marker Plot')+
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}