# PLot expression 
CellTrekFeaturePlot = function(celltrek_obj, feature, pt_size = 0.7, color_max = 0.85, # set max to 85%
                               plot_title ='',
                               # Edges
                               show_edge =T, plot_info_obj, segment_size = 0.2 ){
    # Empty plot
    p = ggplot()
    # Edge
    if(show_edge){
          p = p + geom_segment(data = plot_info_obj$Edge_table, aes(x = x1, y= y1, xend = x2, yend = y2), size = segment_size, color = 'gray70')  # Edge
      }
    # Get coord and expression
    coord_table = plot_info_obj$meta_selected[,c('coord_x','coord_y')]
    expr_table  = FetchData(celltrek_obj, vars = feature)
    expr_max    = expr_table[, feature] %>% max
    expr_coord_table = bind_cols(coord_table[rownames(expr_table), ], expr_table)
    
    # Plot
    p = p + geom_point(data= expr_coord_table, aes(y = -coord_x, x = coord_y, color = .data[[feature]])) + 
        geom_point(size = pt_size) + 
        scale_color_gradientn(colours = c('gray30',Seurat:::SpatialColors(100)), limits = c(0, expr_max*color_max)) + 
        cowplot::theme_cowplot() + 
        labs(title = plot_title) + 
        theme(aspect.ratio = 1)
    p
}
