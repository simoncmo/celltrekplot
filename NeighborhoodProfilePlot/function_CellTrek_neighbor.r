library(igraph)
library(tripack)
library(viridis)


################################################
## Delauney and Neighborhood Percentage
## Added 5/12/2022
################################################
## DEMO
# cell_column = 'cell_type_Abbr' # Column for the cell type
# celltype_plt = c('PT(dual_identity)','PT(S1+S2)','PT(S3)') # remove this value if want to plot All cell types 
# obj_graph = MakeCellTrekDelaunayGraph(celltrek_obj)
# MakeCellNeighborhoodRatioPlot(obj = celltrek_obj, 
#                          obj_graph = obj_graph,
#                          cell_column  = cell_column, 
#                          celltype_plt = celltype_plt) 
################################################
#' Make Delaunay Graph object from Celltrek result or Visium ST obj
#'
#' This function create Delaunay Graph object from Celltrek result or Visium ST obj
#' @param obj The CellTrek of ST object
#' @param trim_graph Whether trim the graph based on maximum distance allow [default = T]
#' @param max_distance Maximum distance to make a connection/edge [defulat = T]
#' @param obj_type Specify whether its a CellTrek or Visium ST object, takes either 'CellTrek','ST' [default = 'CellTrek']
#' @return A igraph obj
#' @export
#' @examples
#' MakeCellTrekDelaunayGraph(celltrek_obj, trim_graph = T, max_distance = 800, obj_type = 'CellTrek')
MakeCellTrekDelaunayGraph = function (obj, trim_graph = T, max_distance = 800, obj_type = c('CellTrek','ST')) 
{
    # Choose type
    obj_type = match.arg(obj_type)
    message(paste('Using', obj_type, 'mode.'))
    # Get coord
    coord_df <- if(obj_type == 'CellTrek') obj@meta.data[, c("coord_x", "coord_y")] else GetTissueCoordinates(obj)
    coord_df = coord_df %>% setNames(c("x", "y"))
    obj_trimesh = tri.mesh(coord_df)
    obj_neighbor = neighbours(obj_trimesh)
    obj_adjmtx_list = obj_neighbor %>% setNames(rownames(obj@meta.data))
    obj_graph = graph_from_adj_list(adjlist = obj_adjmtx_list, 
        mode = "all")
    if (trim_graph) 
        obj_graph = TrimDelaunayGraph(obj_graph, coord_df, max_distance)
    return(obj_graph)
}

## Get n-hop, n-order neighbor
GetAllNeighbors = function(obj, obj_graph, quary_pt, n_order = 1, exclude_quary = T){
    all_ids = rownames(obj@meta.data)
    pt_idx = match(quary_pt, rownames(obj@meta.data)) # Get Index of spots
    # Get neighbor based on size
    neighbor_idx = ego(obj_graph, nodes = pt_idx, order = n_order) %>% unlist
    # get neighbor. If exclude_quary, will exclude quary point themselves 
    neighbor_pt  = if(!exclude_quary) all_ids[neighbor_idx] else all_ids[setdiff(neighbor_idx, pt_idx)] #
    return(neighbor_pt)
}
## Get n-hop, n-order neighbor by ident
GetNeighborByIdent = function(obj, obj_graph, group.by, idents, n_order = 1, exclude_self = T){
    # Extract quary spots 
    quary_pt = obj@meta.data %>% filter(.data[[group.by]] %in% idents) %>% rownames()
    GetAllNeighbors(obj, obj_graph, quary_pt, n_order, exclude_quary = exclude_self)
}

########################################
## 5/17/2022
## Neighborhood Ratio BarPlot
########################################
MakeCellNeighborhoodRatioPlot = function(obj, obj_graph, n_order = 1, cell_column, celltype_plt, neighbor_celltype, exclude_self=T, 
                                        plot_type = c('bar','matrix'), palette, sort_by_percent = T){
    if(missing(celltype_plt)) celltype_plt = obj@meta.data[,cell_column]  %>% unique # Plot All Cell Types
    # Step 1 Get All neighbor points 
    celltype_neighbor_list = map(celltype_plt, function(celltype){
        GetNeighborByIdent(obj, obj_graph, n_order = n_order, group.by = cell_column, idents = celltype, exclude_self=exclude_self) # exclude_self : remove query point themself
    }) %>% setNames(celltype_plt)
    
    # Step 2 Get Neighborhood Cell Type Percentage
    celltype_percent_df = imap(celltype_neighbor_list, function(ids, celltype){
        obj@meta.data[ids, ] %>% 
           mutate(Source_cell = celltype) %>% 
           count(Source_cell, .data[[cell_column]]) %>%
           mutate(Percent = n/sum(n))
    }) %>% bind_rows()
    # Optional
    if(!missing(neighbor_celltype)) celltype_percent_df = celltype_percent_df %>% filter(.data[[cell_column]] %in% neighbor_celltype)
    # Reorder cells 
    if(sort_by_percent){
        # # Method 1 - Deprecating
        #  celltype_order = celltype_percent_df %>% group_by(.data[[cell_column]]) %>% summarize(Total = sum(Percent)) %>% 
        #                 arrange(Total) %>% pull(.data[[cell_column]]) %>% unique 
        # # Method 2 - Deprecating
        # celltype_order = celltype_percent_df %>% arrange(Percent) %>% pull(Source_cell) %>% unique
        # Method 3 - hclust
        percent_mtx   = celltype_percent_df %>% pivot_wider(id_cols = .data[[cell_column]], names_from = Source_cell, values_from = Percent, values_fill = 0) %>% 
            column_to_rownames(cell_column)
        
        hclust_result = hclust(percent_mtx %>% as.matrix %>% dist)
        celltype_order = hclust_result$labels[hclust_result$order]
        
        # Change order
        celltype_percent_df[['Source_cell']] = factor(celltype_percent_df[['Source_cell']], levels = celltype_order)
        celltype_percent_df[[cell_column]] = factor(celltype_percent_df[[cell_column]], levels = celltype_order)
    }
    # Step3 Plot
    plot_type = match.arg(plot_type)
    p = celltype_percent_df %>% 
        ggplot() + 
            cowplot::theme_cowplot() + 
            theme(aspect.ratio = 1) + 
            RotatedAxis()
    if(plot_type == 'matrix'){
        p = p + geom_tile(aes(x = Source_cell, y = .data[[cell_column]], fill = Percent)) + 
            scale_fill_viridis() 
    }else if(plot_type == 'bar'){
        if(missing(palette)) palette = ''
        target_class = if(sort_by_percent) celltype_order else celltype_percent_df[[cell_column]]
        palette = ExtendColorVector(palette, target_class)
        p = p + geom_bar(aes(x = Source_cell, y = Percent, fill = .data[[cell_column]]), stat = 'identity') +
            scale_fill_manual(values = palette)
    }
    p
}

##############################
## Make Delaunay Plot from CellTrek obj
########## DEMO ##############
# MakeDelaunayPlot(obj = result_list_celltrek$cell_type_Abbr$W12_M, cell_column = cell_column,
#                  cell_plot = c('PT(dual_identity)','PT(S3)','PT(S1+S2)'),
#                 palette = col_cell_type, vertex.size = 1.5)
##############################
## This is deprecating --------
# MakeDelaunayPlot = function(obj, cell_column, palette, cell_plot, flip_x=F, flip_y=T,
#                             vertex.size = 1.2, edge.color = 'gray80',...){
#     if(missing(cell_plot)) cell_plot = obj@meta.data[[cell_column]] %>% unique # plot all 
#     # Select data
#     meta_selected = obj@meta.data %>% filter(.data[[cell_column]] %in% cell_plot)
#     cell_ids      = meta_selected %>% rownames()
#     # Adjust coordinate. Since xy switch by default. flip_y -> -coord_y
#     if(flip_x) meta_selected$coord_y = -meta_selected$coord_y
#     if(flip_y) meta_selected$coord_x = -meta_selected$coord_x
#     # Make Graphobj
#     plt_graph  = MakeCellTrekDelaunayGraph(obj[,cell_ids])
#     # Color
#     if(missing(palette)){ 
#         palette = ExtendColorVector(color_vector = '', 
#                                     target_class = obj@meta.data[[cell_column]])}
#     color_cell_type = palette[ meta_selected[[cell_column]] ]
#     # Set color
#     V(plt_graph)$color = color_cell_type
#     # Set layout
#     plt_layout = meta_selected[cell_ids,c('coord_y','coord_x')] %>% as.matrix()
#     # Plot
#     plot(plt_graph, vertex.label =NA, vertex.size = vertex.size, 
#          vertex.frame.color=NA, edge.color = 'gray80', layout = plt_layout, ...)
# }

##############################
## Simpler function to plot Delaunay obj - but uses PLOT instead of GGPLOT
## This is deprecating --------..
# PlotDelaunayGraph = function(obj, plt_graph, cell_column, palette, cell_plot, flip_x = F, flip_y = T, 
#     vertex.size = 1.2, edge.color = "gray80", ...) 
# {
#     if (missing(cell_plot)) 
#         cell_plot = obj@meta.data[[cell_column]] %>% unique
#     meta_selected = obj@meta.data %>% filter(.data[[cell_column]] %in% 
#         cell_plot)
#     cell_ids = meta_selected %>% rownames()
#     if (flip_x) 
#         meta_selected$coord_y = -meta_selected$coord_y
#     if (flip_y) 
#         meta_selected$coord_x = -meta_selected$coord_x
    
#     if (missing(palette)) {
#         palette = ExtendColorVector(color_vector = "", target_class = obj@meta.data[[cell_column]])
#     }
#     color_cell_type = palette[meta_selected[[cell_column]]]
#     V(plt_graph)$color = color_cell_type
    
#     # Set up layout for plt
#     plt_layout = meta_selected[cell_ids, c("coord_y", "coord_x")] %>% 
#         as.matrix()
    
#     # Plot
#     plot(plt_graph, vertex.label = NA, vertex.size = vertex.size, 
#         vertex.frame.color = NA, edge.color = "gray80", layout = plt_layout, 
#         ...)
# }
##############################
## 5/17/2022
## Simpler function to plot Delaunay obj use GGPLOT version
## Use This instead!
##############################
## DEMO
# cell_column = 'cell_type_Abbr' # Column for the cell type
# obj_graph = MakeCellTrekDelaunayGraph(celltrek_obj)
# PlotDelaunayGraphGG(obj, obj_graph, cell_column = cell_column, palette = col_cell_type)
##############################

# PlotDelaunayGraphGG = function(obj, plt_graph, cell_column, palette, cell_plot, flip_x = F, flip_y = T, 
#     title_size = 15, segment_size = 0.2, pt_size = 0.7, 
#     vertex.size = 1.2, edge.color = "gray80", ...) 
# {
#     if (missing(cell_plot)) 
#         cell_plot = obj@meta.data[[cell_column]] %>% unique
#     meta_selected = obj@meta.data %>% filter(.data[[cell_column]] %in% cell_plot)
#     cell_ids = meta_selected %>% rownames()
#     if (flip_x) 
#         meta_selected$coord_y = -meta_selected$coord_y
#     if (flip_y) 
#         meta_selected$coord_x = -meta_selected$coord_x
    
#     if (missing(palette)) {
#         palette = ExtendColorVector(color_vector = "", target_class = obj@meta.data[[cell_column]])
#     }
#     color_cell_type = palette[meta_selected[[cell_column]]]
#     V(plt_graph)$color = color_cell_type
    
#     # Set up layout for plt
#     plt_layout = meta_selected[cell_ids, c("coord_y", "coord_x")]  %>% setNames(c('x','y'))
    
#     # Get Edge/segment df
#     edge_df = GetDelaunayEdgeDistanceTable(plt_graph, plt_layout)

#     # Plot
    
#     ggplot() + 
#         geom_segment(data = edge_df, aes(x = x1, y= y1, xend = x2, yend = y2), size = segment_size, color = 'gray60') + 
#         geom_point(  data = meta_selected, aes(x = coord_y, y = coord_x, color = .data[[cell_column]]), size = pt_size) + 
#         scale_color_manual(values = color_cell_type) + 
#         labs(title = 'Delaunay Triangulation Graph') + 
#         theme_void() + 
#         theme(aspect.ratio = 1, plot.title  = element_text(size = title_size, face = 'bold', hjust = 0.5)) 
# }



############################################################
## Internal Helper Function for Getting Edge and Trim Delaunay Graph
## Subroutine - Get Edge distacne df
############################################################
GetDelaunayEdgeDistanceTable = function(graph, layout){
    # Set up X Y Coordinate
    coord_match_df = layout %>% as.data.frame
    coord_x  = coord_match_df$x %>% setNames(seq(1, nrow(coord_match_df)))
    coord_y  = coord_match_df$y %>% setNames(seq(1, nrow(coord_match_df)))
    cell_ids = rownames(coord_match_df) %>% setNames(seq(1, nrow(coord_match_df)))
    
    # Get Edge table and x, y 
    edge_df = get.edgelist(graph) %>% 
        as.data.frame %>% 
        setNames(c('from_idx', 'to_idx')) %>% 
        mutate(edge_name = str_c(from_idx, to_idx ,sep='|')) %>% 
        mutate(x1 = coord_x[from_idx], y1 = coord_y[from_idx]) %>% 
        mutate(x2 = coord_x[to_idx], y2 = coord_y[to_idx]) %>% 
        mutate(distance = sqrt( ( x1 - x2 ) ^ 2 + ( y1 - y2 ) ^ 2 )) %>% 
        mutate(cell1 = cell_ids[from_idx], cell2 = cell_ids[to_idx])     # Add Cell ID for To and From as Well
    return(edge_df)
}

## Subroutine to trim the graph by edge size
TrimDelaunayGraph = function(graph, layout, max_distance = 500){
    # Get distance
    graph_dist_df = GetDelaunayEdgeDistanceTable(graph, layout)
    # Select edges to remove
    edge_remove = graph_dist_df %>% filter(distance > max_distance) %>% pull(edge_name)
    # Get new graph
    graph %>% delete_edges(edge_remove)    
}
