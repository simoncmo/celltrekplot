# PLot Identity
# Simple function only takes Celltrek object and identity
CellTrekSpatialDimPlot = function(celltrek_obj, group.by, pt_size = 0.7, 
                               plot_title ='',
                               # Label
                               label = F, label_size = 5
                               # Edges
                               #show_edge =T, plot_info_obj, segment_size = 0.2 
                                 ){
    # Empty plot
    p = ggplot()
    
    # Get coord and expression
    coord_table = celltrek_obj@meta.data[,c('coord_x','coord_y')]
    ident_table = FetchData(celltrek_obj, vars = group.by)
    ident_coord_table = bind_cols(coord_table[rownames(ident_table), ], ident_table)
    
    # Plot
    p = p + geom_point(data= ident_coord_table, aes(y = -coord_x, x = coord_y, color = .data[[group.by]]),
                      size = pt_size) + 
        #scale_color_gradientn(colours = c('gray30',Seurat:::SpatialColors(100)), limits = c(0, expr_max*color_max)) + 
        cowplot::theme_cowplot() + 
        labs(title = plot_title) + 
        theme(aspect.ratio = 1)
    # Label coordinates
    if(label){
        label_coord =  ident_coord_table %>% group_by(.data[[group.by]]) %>% 
            summarize(coord_x = mean(coord_x), coord_y = mean(coord_y))
    p = p + geom_text(data = label_coord, aes(y = -coord_x, x = coord_y, label = .data[[group.by]]),
                     size = label_size)
    }
    
    p
}
