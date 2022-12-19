### Plot ST and CellTrek together 
CoPlotSTCellTrek = function(celltrek_obj, st_obj, ptsize = 0.1, stptsize = 1){
    # Make table
    ST_table = st_obj@images[[1]]@coordinates[,c('imagerow','imagecol')] %>% 
        setNames(c('coord_x','coord_y')) %>% 
        mutate(assay='ST')
    
    celltrek_table = celltrek_obj@meta.data[,c('coord_x','coord_y')] %>% mutate(assay = 'celltrek')
    # Flip y to reorient
    ST_table$coord_x = -ST_table$coord_x
    celltrek_table$coord_x = -celltrek_table$coord_x
    
    # Plot
    ggplot(data.frame(), aes(x = coord_y, y = coord_x, color = assay)) + 
        geom_point(data = ST_table, size = stptsize, alpha = 0.3) + 
        geom_point(data = celltrek_table, size = ptsize) + 
        cowplot::theme_cowplot() + 
        coord_fixed()
}

PlotBothStCellTrek = function(celltrek_obj, st_obj){
    
}