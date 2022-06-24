################################
## Spatial Cluster Plot
################################
## DEMO
# sample = '20210129-AKI3M'
# cell_type_plt = 'Havcr1+Krt20+'
# km_use = 3
# palette_use_cluster = MakeClusterPalette(k = km_use, celltype = cell_type_plt)

# # Step 1"
# tmp_plt_info_obj = PrepareSpatialPlotObject(celltrek_list[[sample]], 
#                                             graph_list[[sample]],
#                                             km_use = km_use,
#                                             palette_cluster = palette_use_cluster,
#                                             cell_column = cell_column,
#                                             palette = col_cell_full, celltype_highlight = cell_type_plt, n_order = 1)

### Step 2 - Make 3 different kind of plot
# options(repr.plot.width = 20, repr.plot.height = 10)
# p0 = MakeDelaunayGraphPlot(tmp_plt_info_obj, cell_column = cell_column) 
# p1 = p0 %>% AddClusterHighlight(tmp_plt_info_obj)
# p2 = p0 %>% AddDistanceToPoint(tmp_plt_info_obj, center_pt = c(13000,-10000)) 

############ API ########################
# Step 1 Generate Plot obj contain - table and palette for plotting
# THe Spatial plot function uses Table and Palette here
# Tables: Neighbor_table''Target_table''Background_table''Edge_table''Neighbor_cluster_table''meta_selected', Neighbor_cluster_table_full
# Parameters: Target_celltype, n_order, cell_column
# Palettes: palette_celltype, palette_cluster, palette_w_cluster

PrepareSpatialPlotObject = function(obj, plt_graph, cell_column, 
                                    cell_plot, flip_x = F, flip_y = T, celltype_highlight=NA, n_order=3, km_use=3,# For MakeSpatialTableList
                                   palette = '', palette_cluster = '' # For GetSpatialPaletteList
                                   ){
    # Make table list
    table_list = MakeSpatialTableList(obj, plt_graph, cell_column, cell_plot, flip_x, flip_y, celltype_highlight, n_order, km_use) 
    # Make palette list
    palette_list = GetSpatialPaletteList(table_list, palette, palette_cluster)
    #return(table_list)
    #return(palette_list)
    result_obj = c(table_list, list('Palettes' = palette_list))
    return(result_obj)
}


########### Internal ####################
## STEP 1a - Get Table For the Plot
MakeSpatialTableList = function(obj, plt_graph, cell_column, cell_plot, flip_x = F, flip_y = T, 
                                      celltype_highlight=NA, n_order=2, km_use=3)  # setdefaul km = 3
{
    #####################
    ## Select cell to plot
    if (missing(cell_plot)) 
        cell_plot = obj@meta.data[[cell_column]] %>% unique
    meta_selected = obj@meta.data %>% filter(.data[[cell_column]] %in% cell_plot)
    
    #####################
    ## Coordinate,
    if (flip_x) 
        meta_selected$coord_y = -meta_selected$coord_y
    if (flip_y) 
        meta_selected$coord_x = -meta_selected$coord_x
    
    #####################
    # Check if cell exist 
    if(!celltype_highlight %in% meta_selected[[cell_column]]) stop('Cell type not found')
    
    #####################
    # Get Target cell cluster (Based on Neighbor Profile)
    neighbor_cluster_df = GetNeighborCellTypeTable(obj, plt_graph, cell_column) %>% 
        filter(Origin_celltype %in% c(celltype_highlight)) %>% ClusterAndReorderTable(n_km = km_use)
    # Get ALL Tables
    table_list = MakeMetadataClusterTable(obj, plt_graph, meta_selected, n_order, cell_column, celltype_highlight, neighbor_cluster_df)
    # Add meta_select. keep full table
    table_list = c(table_list, list(meta_selected = meta_selected,
                                    Neighbor_cluster_table_full = neighbor_cluster_df
                                   ))
    # Add Parameters used to generate table 'Target celltype'
    parameter_list = list(Target_celltype = celltype_highlight, 
                          cell_column = cell_column,
                            n_order = n_order,
                          km_use=km_use
                         )
    table_list = c(table_list, list('Parameters' = parameter_list))
    
    return(table_list)
}

## Step 1b Get Palette list from Table list
GetSpatialPaletteList = function(table_list, palette='', palette_cluster=''){
    cell_column = table_list$Parameters$cell_column
    #####################
    palette_list = list()
    # Palette 
    palette_list$palette_celltype = ExtendColorVector(color_vector = palette, target_class = table_list$meta_selected[[cell_column]])
    # Cluster Palette 
    palette_list$palette_cluster = ExtendColorVector(palette_list$palette_cluster, target_class =  table_list$Neighbor_cluster_table[['Cell_cluster']])
    ## Color 
    palette_list$palette_w_cluster = c(palette_list$palette_celltype, palette_list$palette_cluster)
    return(palette_list)
}

############ API ########################
# STEP 2 - Base Plot
## Normal Delaunay Highlight cell and neighbor
## 6/17/2022
MakeDelaunayGraphPlot = function(plot_info_obj, sample_name,
     # Toggle parameters
    show_all_points = F, # Show all point of CellTrek output. Will suppress all parameter other than show_edge below                       
    show_edge = T, show_background_pts=T, show_neighbor_pts=T, show_target_pts=T, 
    title_size = 15, segment_size = 0.2, pt_size = 0.7, 
    vertex.size = 1.2, edge.color = "gray80", ...) 
{
    meta_selected = plot_info_obj$meta_selected
    #####################
    # Plot Parameter
    cell_column = plot_info_obj$Parameters$cell_column
    cell_highlight = plot_info_obj$Parameters$Target_celltype
    # Format title
    plt_title = if(show_all_points) 'CellTrek spots' else str_glue('Highlight {cell_highlight} + neighbor of {plot_info_obj$Parameters$n_order} orders')
    if(!missing(sample_name)) plt_title = str_c(sample_name, plt_title, sep=': ')
    
   
    #####################
    # Palettes
    palette_celltype   = plot_info_obj$Palettes$palette_celltype
    palette_cluster    = plot_info_obj$Palettes$palette_cluster
    palette_w_cluster  = plot_info_obj$Palettes$palette_w_cluster
    
    # Plot
    p = ggplot()
    ## 1. EDGE
    p = if(!show_edge) p else p + geom_segment(data = plot_info_obj$Edge_table, aes(x = x1, y= y1, xend = x2, yend = y2), size = segment_size, color = 'gray70')  # Edge
    if(show_all_points){
        p = p + geom_point(data = plot_info_obj$meta_selected, aes(x = coord_y, y = coord_x, color = .data[[cell_column]]), size = pt_size*2, shape = 20,)  # All CellTrekPoints
        p = p + scale_color_manual(values = palette_celltype)
    }else{
        ## 2. Background pts 
        p = if(!show_background_pts) p else p + geom_point(  data = plot_info_obj$Background_table, aes(x = coord_y, y = coord_x), size = pt_size, color = 'gray70', alpha = 0.2)  # none highlight pts
        ## 3. Neighbor pt
        p = if(!show_neighbor_pts) p else p + geom_point(  data = plot_info_obj$Neighbor_table  , aes(x = coord_y, y = coord_x, color = .data[[cell_column]]), size = pt_size*4, shape = 21, stroke = 0.8, fill ='white')  # Neighbor pts
        ## 4. Target cell pt
        p = if(!show_target_pts) p else p + geom_point(  data = plot_info_obj$Target_table , aes(x = coord_y, y = coord_x, fill = .data[[cell_column]]), size = pt_size*5,  shape = 21, stroke = 0.2, color = 'gray10')  # Neighbor pts
        ## Palettes
        p = p + scale_fill_manual(values = palette_w_cluster) + 
            scale_color_manual(values = palette_w_cluster) 
    }
    ## 6. Others 
    p + labs(title = plt_title) + 
        theme_void() + 
        theme(aspect.ratio = 1, plot.title  = element_text(size = title_size, face = 'bold', hjust = 0.5)) 
}

############ API ########################
# Step 3 - Optional, highlight Cluster of Taret Cell type
AddClusterHighlight = function(p_st, plot_info_obj, pt_size = 2.5, stroke = 0.2, remove_neighbor_cells = T){
    ## 5. Cluster
    cell_highlight = plot_info_obj$Parameters$Target_celltype
    target_cluster_stroke_color =  plot_info_obj$Palettes$palette_celltype[[cell_highlight]]
    p_st = p_st + geom_point(  data = plot_info_obj$Neighbor_cluster_table,    # Cell Cluster pts
                      aes(x = x, y = y, fill = Cell_cluster), size = pt_size,  
                      shape = 21, stroke = stroke, color = target_cluster_stroke_color) 
    # Remove neighbor cells to make it cleaner
    # Neighbor cell is layer 3
    if(remove_neighbor_cells){
        p_st$layers = p_st$layers[setdiff(1:length(p_st$layers), 3)] # remove layer 3
    }
    return(p_st)
}


############ API ########################
## Step 4: [Optional] - Add distance to a designated Point 
AddDistanceToPoint = function(p_st, plot_info_obj ,center_pt = c(0,0),
                              pt_size = 1, show_neighbor_pts = F
                             ){
    ### Target cell's cluster table 
    target_cluster_df = plot_info_obj$Neighbor_cluster_table
    target_celltype   = plot_info_obj$Parameters$Target_celltype
    
    # Palettes 
    #palette_celltype   = plot_info_obj$Palettes$palette_celltype
    palette_cluster    = plot_info_obj$Palettes$palette_cluster
    palette_w_cluster  = plot_info_obj$Palettes$palette_w_cluster
    
    # Get coord + target point distance table
    coord_target_df = GetCoordToPointDistance(target_cluster_df, target_point = center_pt, x_column = 'x', y_column = 'y') 
    
    # Box Plot
    p_box = DistanceBoxPlot(coord_target_df, group.by = 'Cell_cluster') + 
        scale_fill_manual(values = palette_cluster) + 
       scale_color_manual(values = palette_cluster) #+ RotatedAxis()
    
    # Format title
    plt_title = str_glue('Target cell {target_celltype} cluster to point {toString(center_pt)}')
    
    # Plot
    p_st = p_st +
        # Line : Target to Point
        geom_curve(data = coord_target_df, aes(x = x, y = y, xend = target_x, yend = target_y, colour = Cell_cluster, alpha = (1/log10(distance))/2), size = 0.5, curvature = 0) + 
            scale_alpha_continuous(range = c(0.1,0.8)) +
        # Target Cluster
        geom_point(  data = target_cluster_df  , aes(x = x, y = y, fill = Cell_cluster), size = pt_size*3,  shape = 21, stroke = 0.2, color = 'gray10') + # Cell Cluster pts
        # Point
        geom_point(  data = data.frame(x = center_pt[[1]], y = center_pt[[2]]), aes(x=x,y=y), size = pt_size*5, shape = 23, stroke = 1, color = 'red', fill = 'blue') +  # Center point
        #scale_fill_manual(values = palette_w_cluster) + 
        #scale_color_manual(values = palette_w_cluster) + 
        labs(title = plt_title) + 
        cowplot::theme_cowplot() + 
        theme(aspect.ratio = 1, plot.title  = element_text(size = 10, face = 'bold', hjust = 0.5)) 
    
    # Remove neighborhood layer - layer 3
    if(!show_neighbor_pts) p_st$layers = p_st$layers[setdiff(1:length(p_st$layers), 3)]
    
    p_st|p_box
    
}



############ INTERNAL #########################
## Interal for making table list
MakeMetadataClusterTable = function (obj, obj_graph, meta_selected, n_order, cell_column, 
    celltype_highlight, neighbor_cluster_df) 
{
    table_ids_list = list(Neighbor = GetNeighborByIdent(obj, 
        obj_graph, n_order = n_order, group.by = cell_column, 
        idents = celltype_highlight, exclude_self = TRUE), Target = meta_selected %>% 
        filter(.data[[cell_column]] %in% celltype_highlight) %>% 
        rownames())
    table_list = SplitMetaByIdList(meta_selected, table_ids_list)
    plt_layout = meta_selected[, c("coord_y", "coord_x")] %>% 
        setNames(c("x", "y"))
    # Edge
    edge_df = GetDelaunayEdgeDistanceTable(obj_graph, plt_layout)
    # Edge - Add cell type
    
    table_list = c(table_list, list(Edge_table = edge_df))
    if (is.null(neighbor_cluster_df)) 
        return(table_list)
    neighbor_cluster_df = neighbor_cluster_df[, c("Origin", "Origin_celltype", 
        "kmeans")] %>% distinct %>% mutate(Cell_cluster = paste0(Origin_celltype, 
        "- clst", kmeans))
    neighbor_cluster_df = left_join(neighbor_cluster_df, plt_layout %>% 
        rownames_to_column("Origin"), by = "Origin")
    table_list = c(table_list, list(Neighbor_cluster_table = neighbor_cluster_df))
    return(table_list)
}

## Get Meta split by Target Id and Neighbor Ids
SplitMetaByIdList = function(meta, list_of_ids){
    background_ids = setdiff(rownames(meta), unlist(list_of_ids))
    list_of_ids = c(list_of_ids, list(Background = background_ids))
    table_list = map(list_of_ids, function(idx) meta[idx, ]) %>% setNames(str_c(names(.),'_table'))
    return(table_list)
}
                     
# Get Distance of each point in df to target point
# Return in Table format
GetCoordToPointDistance = function(coord_df, target_point = c(0,0), x_column = 'coord_x', y_column = 'coord_y'){
    coord_df %>% mutate(target_x = target_point[[1]], target_y = target_point[[2]],
                        distance = sqrt( (.data[[x_column]] - target_x)^2 + (.data[[y_column]] - target_y)^2 ))
}
                     

# Make Distance Boxplot from Coord_df with target distance 
DistanceBoxPlot = function(coord_target_df, group.by ){
    coord_target_df[[group.by]] = as.character(coord_target_df[[group.by]])
    ggplot(data = coord_target_df, aes(x = .data[[group.by]], y = distance, fill = .data[[group.by]])) + 
    geom_boxplot(outlier.shape = NA, color = 'gray30', alpha = 0.4, width = 0.5) + 
    geom_jitter(aes(color = .data[[group.by]])) + 
    cowplot::theme_cowplot() + 
    labs(title = 'Distance distribution') + 
    theme(aspect.ratio = 1)
}
