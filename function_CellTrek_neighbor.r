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
# Delany triangulation 
MakeCellTrekDelaunayGraph = function(obj){
    coord_df <- obj@meta.data[,c('coord_x','coord_y')] %>% setNames(c('x','y'))
    obj_trimesh = tri.mesh(coord_df)
    obj_neighbor = neighbours(obj_trimesh) # Get neighbor list
    
    # Turn into ids and adjacency matrix
    obj_adjmtx_list = obj_neighbor %>% setNames(rownames(obj@meta.data)) 
    # Make Graph
    obj_graph = graph_from_adj_list(
        adjlist =  obj_adjmtx_list, mode = 'all'
    )
    return(obj_graph)
}

## Get n-hop, n-order neighbor
GetAllNeighbors = function(obj, obj_graph, quary_pt, n_order = 1){
    all_ids = rownames(obj@meta.data)
    pt_idx = match(quary_pt, rownames(obj@meta.data)) # Get Index of spots
    # Get neighbor based on size
    neighbor_idx = ego(obj_graph, nodes = pt_idx, order = n_order) %>% unlist
    neighbor_pt  = all_ids[neighbor_idx]
    return(neighbor_pt)
}
## Get n-hop, n-order neighbor by ident
GetNeighborByIdent = function(obj, obj_graph, group.by, idents, n_order = 1){
    # Extract quary spots 
    quary_pt = obj@meta.data %>% filter(.data[[group.by]] %in% idents) %>% rownames()
    GetAllNeighbors(obj, obj_graph, quary_pt, n_order)
}

MakeCellNeighborhoodRatioPlot = function(obj, obj_graph, n_order = 1, cell_column, celltype_plt, neighbor_celltype,
                                        plot_type = c('bar','matrix'), palette){
    if(missing(celltype_plt)) celltype_plt = obj@meta.data[,cell_column]  %>% unique # Plot All Cell Types
    # Step 1 Get All neighbor points 
    celltype_neighbor_list = map(celltype_plt, function(celltype){
        GetNeighborByIdent(obj, obj_graph, n_order = n_order, group.by = cell_column, idents = celltype)
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
        palette = ExtendColorVector(palette, celltype_percent_df[[cell_column]])
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
MakeDelaunayPlot = function(obj, cell_column, palette, cell_plot, flip_x=F, flip_y=T,
                            vertex.size = 1.2, edge.color = 'gray80',...){
    if(missing(cell_plot)) cell_plot = obj@meta.data[[cell_column]] %>% unique # plot all 
    # Select data
    meta_selected = obj@meta.data %>% filter(.data[[cell_column]] %in% cell_plot)
    cell_ids      = meta_selected %>% rownames()
    # Adjust coordinate. Since xy switch by default. flip_y -> -coord_y
    if(flip_x) meta_selected$coord_y = -meta_selected$coord_y
    if(flip_y) meta_selected$coord_x = -meta_selected$coord_x
    # Make Graphobj
    plt_graph  = MakeCellTrekDelaunayGraph(obj[,cell_ids])
    # Color
    if(missing(palette)){ 
        palette = ExtendColorVector(color_vector = '', 
                                    target_class = obj@meta.data[[cell_column]])}
    color_cell_type = palette[ meta_selected[[cell_column]] ]
    # Set color
    V(plt_graph)$color = color_cell_type
    # Set layout
    plt_layout = meta_selected[cell_ids,c('coord_y','coord_x')] %>% as.matrix()
    # Plot
    plot(plt_graph, vertex.label =NA, vertex.size = vertex.size, 
         vertex.frame.color=NA, edge.color = 'gray80', layout = plt_layout, ...)
}