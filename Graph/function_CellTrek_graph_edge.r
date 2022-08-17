# Edge related function
# update 6/19/2022

#########################
## Example
## plot_info_obj = GetEdgeCellTypeTable(plot_info_obj) %>% AddEdgeGroup
#########################



### API ###################################
# 1. Get Edge Type 
GetEdgeCellTypeTable = function(plot_info_obj, cell_group = c('Tumor'), others_name = 'TME'){
    # Parameter
    celltype_column = plot_info_obj$Parameters$cell_column
    # Celltype
    celltypes = plot_info_obj$meta_selected[[cell_column]]
    # Generate Edge table
    edge_celltype_table = plot_info_obj$Edge_table %>% 
        mutate(cell1_celltype = celltypes[from_idx], 
               cell2_celltype = celltypes[to_idx]) %>% 
        mutate(Connect_type = ifelse(cell1_celltype == cell2_celltype, 'Cis','Trans')) %>%
        rowwise() %>% 
        mutate(Edge_celltype = str_c(sort(c(cell1_celltype, cell2_celltype)), collapse="-")) 

    ## Edge summary table
    edge_summary_list = list(CellType = edge_celltype_table %>% count(Connect_type, Edge_celltype) %>% arrange(desc(n))
                            )
    
    ## Add back
    plot_info_obj$Edge_table   = edge_celltype_table
    plot_info_obj$Edge_summary = edge_summary_list
    message("Save to $Edge_table and Edge_summary$CellType tables in the object")
    return(plot_info_obj)
}


## 2. Get Edge Group
AddEdgeGroup = function(plot_info_obj, cell_group = c('Tumor'), others_name = 'TME'){
    # Parameter
    edge_celltype_table = plot_info_obj$Edge_table
    ## Get Cell Type Group
    edge_celltype_table = GetEdgeInteractionType(edge_celltype_table, cell_group, others_name)

    ## Edge summary table?
    edge_summary_list = list(CellType = edge_celltype_table %>% count(Connect_type, Edge_group, Edge_celltype) %>% arrange(desc(n)),
                             Group   = edge_celltype_table %>% count(Connect_type, Edge_group) %>% arrange(desc(n))
                            )
    
    ## Add back
    plot_info_obj$Edge_table   = edge_celltype_table
    plot_info_obj$Edge_summary = edge_summary_list
    message("Save to $Edge_table and Edge_summary$CellType and Edge_summary$Group tables in the object")
    return(plot_info_obj)
}

########## INternal

# Generate Cell Type Group 
GetEdgeInteractionType = function(edge_celltype_table, cell_group = c('Tumor'), others_name = 'TME'){
    message(str_glue('Generate cell type group using {toString(cell_group)}. Not selected Cell Type as: {others_name}'))
    edge_celltype_table %>% mutate(cell1_celltype_group = ifelse(cell1_celltype %in% cell_group, cell1_celltype, others_name), 
                                   cell2_celltype_group = ifelse(cell2_celltype %in% cell_group, cell2_celltype, others_name)) %>% 
    rowwise() %>% 
    mutate(Edge_group = str_c(sort(c(cell1_celltype_group, cell2_celltype_group)), collapse="-"))
}


########### PLOT Function
#### Edge based plot
### [API] This function make Delaunay graph and color each edge by either node ('cell1' or 'cell2')
EdgePlotByCell = function(plot_info_obj, color_edge_by = c('cell1','cell2'), segment_size= 0.5, title_size=20, plt_title = '',
                                 highlight_neighbor = F, background_edge_color = 'gray70'){
    # Parameters
    cell_column = plot_info_obj$Parameters$cell_column
    palette_celltype = plot_info_obj$Palettes$palette_celltype
    # Celltype
    celltypes = plot_info_obj$meta_selected[[cell_column]]
    # Edge
    Color_edge_table = plot_info_obj$Edge_table %>% 
        mutate(Celltype1 = celltypes[from_idx], Celltype2 = celltypes[to_idx])
    # Choose Cell1 or Cell2 to color edges
    color_edge_by = match.arg(color_edge_by)
    color_column = if(color_edge_by == 'cell1') 'Celltype1' else 'Celltype2'
        
    # Highlight Neighbor (Turn everything else gray)
    background_cells = if(highlight_neighbor) plot_info_obj$Background_table %>% rownames else ''
    Color_edge_table = Color_edge_table %>% mutate({{ color_column }} := ifelse(.data[[color_edge_by]] %in% background_cells, 
                                                                               'background',
                                                                                .data[[color_column]]
                                                                               ))
    
    palette_use = c(palette_celltype, 'background' = background_edge_color)
    # Plot
    p = ggplot()

    p = p + geom_segment(data = Color_edge_table, 
                         aes(x = x1, y = y1, xend = x2, yend = y2, 
                             color = .data[[color_column]] ), 
                         size = segment_size)
    p + scale_fill_manual(values = palette_celltype) + 
        scale_color_manual(values = palette_use) + 
        labs(title = plt_title) + 
    theme_void()+

        theme(aspect.ratio = 1, plot.title = element_text(size = title_size, face = "bold", hjust = 0.5))
}