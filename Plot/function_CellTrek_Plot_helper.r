## PLOT FUNCTIONs for CellTrek
# INTERNAL 
CellTrekPlot = function(data_df, group.by, coord_cols, pt_size = 1, label=T, label_size = 5, label_col = 'black'){
    # Get coord name
    coord_x = coord_cols[[1]]
    coord_y = coord_cols[[2]]
    # Text label 
    text_df = data_df %>%  group_by(.data[[group.by]]) %>% summarize({{coord_x}} := mean(.data[[coord_x]]), 
                                                                     {{coord_y}} := mean(.data[[coord_y]]))
    # Plot
    p = data_df %>% 
        ggplot(aes(x = .data[[coord_x]], y = .data[[coord_y]], color = .data[[group.by]])) + 
        geom_point(size = pt_size) + 
        cowplot::theme_cowplot()+
        theme(aspect.ratio = 1)
    # enlarge legend dotsize. Current doesn't support changing size here
    p = p + guides(color = guide_legend(override.aes = list(size = 5)))
    # label
    if(!label) p else p + geom_text(data = text_df, aes(x = .data[[coord_x]], y =.data[[coord_y]], 
                                      label = .data[[group.by]]), color = label_col, size = label_size)

}

#### API
UmapCellTrekPlot = function(celltrek_obj, group.by, pt_size = 1, label=T, label_size = 5, label_col = 'black', filter_na = F){
    data_df = FetchData(celltrek_obj, c('umap_1','umap_2',group.by)) # Get DATA
    data_df = if(filter_na) data_df %>% filter(!is.na(.data[[group.by]])) else data_df
    CellTrekPlot(data_df, group.by, coord_cols = c('umap_1','umap_2'), pt_size=pt_size, label=label, label_size = label_size, label_col = label_col) # Plot
}

SpatialCellTrekPlot = function(celltrek_obj, group.by, pt_size = 1, label=T, label_size = 5, label_col = 'black', filter_na = F){
    data_df = FetchData(celltrek_obj, c('coord_y','coord_x',group.by)) # Get DATA
    # Flip x and y and add "-"" to new y
    tmp_new_x = data_df$coord_y
    data_df$coord_y = -data_df$coord_x
    data_df$coord_x = tmp_new_x
    # filter
    data_df = if(filter_na) data_df %>% filter(!is.na(.data[[group.by]])) else data_df
    #Plot
    CellTrekPlot(data_df, group.by, coord_cols = c('coord_x','coord_y'), pt_size=pt_size, label=label, label_size = label_size, label_col = label_col)
}

## Added 5/12/2022
## Countour
CellTrekContour = function(obj, cell_column, cell_highlight, show_countour = T, 
                           countour_col = 'blue', fill_countour =F, 
                           flip_x = F, flip_y = T, switch_xy = T, # Coordinates
                           pt_size = 0.5){
    # Adjust x,y 
    cell_df = obj@meta.data
    if(switch_xy){
        tmp = cell_df$coord_x
        cell_df$coord_x = cell_df$coord_y
        cell_df$coord_y = tmp
    }
    if(flip_x) cell_df$coord_x = -cell_df$coord_x
    if(flip_y) cell_df$coord_y = -cell_df$coord_y
    # Select Cell
    if(missing(cell_highlight)){
        highlight_df = cell_df
        plt_title = 'Density of All Spots'
    }else{
        highlight_df = cell_df %>% filter(.data[[cell_column]] %in% cell_highlight)
        plt_title = str_glue('Density of {toString(cell_highlight)}')
    }
    
    # Plot
    p = ggplot(data = cell_df, aes(x=coord_x, y=coord_y, color = .data[[cell_column]])) 
    # Fill
    if(fill_countour){
        p = p + stat_density_2d(data = highlight_df,aes(fill = ..density..), geom = "raster", contour = FALSE) +
            scale_fill_distiller(palette=1, direction=1) 
    }
    # Countour
    if(show_countour){
        p = p + geom_density_2d(data = highlight_df, aes(x=coord_x, y=coord_y), color =countour_col) 
    }
    # Points
    p = p +
        geom_point(size = pt_size)+
        cowplot::theme_cowplot() + 
        theme(aspect.ratio = 1) + 
        labs(title = plt_title)
    
    return(p)
}
