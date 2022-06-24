## DEVeloping
## Step 3 : 
## Make Neighborhood cell type profile plot
## Use plot_info_obj
MakeNeighborProfilePlotV2 = function(plot_info_obj, title_sample = '', reorder_profile=T){
    # Parameters
    Target_celltype = plot_info_obj$Parameters$Target_celltype
    km_use = plot_info_obj$Parameters$km_use
    # Table
    neighbor_df = plot_info_obj$Neighbor_cluster_table_full
    
    if(!all(c('hc_order', 'kmeans') %in% names(neighbor_df))) stop('Missing hclust and keamns. Run ClusterAndReorderTable first')
    # Palette
    palette_use = plot_info_obj$Palettes$palette_celltype
    palette_use = ExtendColorVector(palette_use, target_class = neighbor_df[['Neighbor_celltype']])
    palette_cluster = plot_info_obj$Palettes$palette_cluster
    
    
    
    ###########################
    # Top plot: KM count Plt
    ## check how many cell in each km group
    cluster_count_df = neighbor_df[,c('Origin','kmeans')] %>% ungroup %>% distinct %>% count(kmeans) %>% 
        mutate(kmeans = as.character(kmeans),
            CellType = Target_celltype,
              Percent = n /sum(n),
              Cell_cluster = paste0(all_of(Target_celltype),'- clst',kmeans),
              Label = paste0(Cell_cluster, '\n', n, ' (', round(Percent,3)*100, '%)')) %>% 
        mutate(Cell_cluster = forcats::fct_reorder(Cell_cluster, desc(kmeans)))  # reorder label
    ## plot
    p_kmpercent = cluster_count_df %>% ggplot(aes(x = n, y = CellType, fill = Cell_cluster)) + 
        geom_bar(stat ='identity') + 
        geom_text(aes(label = Label), position = position_stack(vjust = 0.5)) +
        labs(x = 'Number of Target Cell', title= 'Target Cell Cluster') +
        scale_x_continuous(expand = c(0,0)) + 
        scale_fill_manual(values = palette_cluster) +
        cowplot::theme_cowplot()
    
    ########################### 
    # Middle: Individual Profile matrix plt
    p_profilemtx = neighbor_df  %>% 
        ggplot(aes(x = Origin, y = percent, fill = Neighbor_celltype)) + 
            facet_grid(.~kmeans, scales='free', space = 'free')+
            geom_bar(stat = 'identity') + 
            scale_fill_manual(values = palette_use) + 
            labs(x = 'Individaul Target Cell', y = 'Percentage', title = 'Individual Cell Neighbor Cell Type Profile') + 
            cowplot::theme_cowplot()+
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),)
    if(reorder_profile) p_profilemtx = ReorderProfilePlotByCellRatio(p_profilemtx)

    ###############
    ### Bottom : Profile by cluster
    cluster_level_df = neighbor_df %>%
        group_by(Origin_celltype, kmeans, Neighbor_celltype) %>% 
        summarize(total = sum(n) , .groups ='keep') %>% 
        group_by(kmeans) %>% 
        mutate(Percent = total/sum(total)) %>% 
        mutate(Neighbor_celltype = forcats::fct_reorder(Neighbor_celltype, total, .fun = 'sum')) %>%# reorder Cell
        mutate(Label = ifelse(Percent > 0.20, str_glue('{round(Percent, 2)*100}%'), ''))

    p_profile_cluster = cluster_level_df  %>% 
        ggplot(aes(x = total, y = Origin_celltype, fill = Neighbor_celltype)) + 
            facet_grid(.~kmeans, scales='free', space = 'free')+
            geom_bar(stat ='identity') + 
            geom_text(aes(label = Label), position= position_stack(vjust = 0.5)) +
            scale_x_continuous(expand = c(0,0)) + 
            labs(x = 'Number of Neighbor Cell Types', y = '', title = 'Cluster Level Neighbor Cell Type Profile') + 
            scale_fill_manual(values = palette_use) +
            cowplot::theme_cowplot() & NoLegend()
    
    p = p_kmpercent + p_profilemtx + p_profile_cluster + plot_layout(height = c(1,2,1)) +  
        plot_annotation(title = str_glue('{title_sample} - {Target_celltype}: Neighbor Cell Profile'),
                       theme = theme(plot.title = element_text(face = 'bold', size = 20)))
    p
}




########################################################################
## Internal Functinos
########################################################################

# This function reorder Cellid and Cell type in p_profile plot 
ReorderProfilePlotByCellRatio = function(p_profile){
    # Select top cell in each group, order by Cell by it
    kmean_groups = p_profile$data$kmeans %>% unique
    # Get 'Weight' for each group
    cell_weight_table = map(kmean_groups, function(k){
        tmp_df = p_profile$data %>% 
            filter(kmeans == k) %>% 
            group_by(kmeans, Neighbor_celltype) %>% summarize(total = sum(n), .groups='keep') %>% 
            arrange(desc(total))
        tmp_df %>% ungroup %>% mutate(Weight = 20^(nrow(tmp_df):1)) # Add weight to help sorting   
    }) %>% bind_rows()

    ## Step1.1 : Cell order
    tmp_cell_order_table = cell_weight_table %>% group_by(Neighbor_celltype) %>% summarize(total = sum(total)) %>% arrange(desc(total))
    tmp_cell_order = tmp_cell_order_table$Neighbor_celltype

    ## Step2: Calculate weight
    new_order_table = left_join(p_profile$data, cell_weight_table, by=c('kmeans','Neighbor_celltype')) %>%
        mutate(Score = percent*Weight) %>% 
        group_by(kmeans, Origin) %>% 
        summarize(Per_cell_score = sum(Score), .groups = 'keep') %>% 
        arrange(kmeans, desc(Per_cell_score)) %>%
        mutate(Origin = as.character(Origin))
    new_id_order =  new_order_table %>% pull(Origin)

    ## Step3: reoder data
    p_profile$data = p_profile$data %>% mutate(Origin = factor(Origin, levels = new_id_order))
    p_profile$data = p_profile$data %>% mutate(Neighbor_celltype = factor(Neighbor_celltype, levels = tmp_cell_order))
    return(p_profile)
}

##################################
# Below moved from function_CellTrek_neighbor_plot. R 
# Need reivew what to keep
# 6/20/2022
##################################


## 5/18/2022
## Individual Cell Neighborhood analysis
## To show DMG cell location/neighboring cells

library(patchwork)
######### INTERNAL ####################
## Get Connection table. Origin = Cell1 id, Neighbor = Cell2 id, Connect column just mean these 2 point connect
GetConnectionTable = function(obj, graph_obj){
    ########################
    # Extract Neighboring cells
    ########################
    mtx = igraph::as_adjacency_matrix(graph_obj)
    ## https://stackoverflow.com/questions/12029177/sparse-matrix-to-a-data-frame-in-r
    ## reformat sparce matrix
    summ = Matrix::summary(mtx)
    cell_ids = rownames(obj@meta.data)
    connection_df = data.frame(
               Origin        = cell_ids[summ$i],
               Neighbor      = cell_ids[summ$j],
               Connect       = ifelse(summ$x==1,T,F)) # will all be true. i is cell1, j is cell2
    connection_df
}
#GetConnectionTable(ST, st_graph) 

######### API #########################
# Step 1 Get neighborhood table
GetNeighborCellTypeTable = function(obj, graph_obj, cell_column){
    ########################
    # Extract Neighboring cells
    ########################
    connection_df = GetConnectionTable(obj, graph_obj)

    ########################
    # Get Cell Type count and percentage
    ########################
    # Get cell type
    meta_cell_df = obj@meta.data[,cell_column, drop=F] %>% rownames_to_column('Neighbor')
    neighbor_celltype_df = left_join(connection_df, meta_cell_df, by ='Neighbor')

    # Convert to cell count/percent
    neighbor_celltype_count_df = neighbor_celltype_df %>% group_by(Origin) %>% 
        count(.data[[cell_column]]) %>% mutate(percent = n / sum(n))

    # Rename Cell Column
    neighbor_celltype_count_df = neighbor_celltype_count_df %>% 
        dplyr::rename(Neighbor_celltype = .data[[cell_column]])
    
    # Add Origin Cell Type
    meta_cell_origin_df = meta_cell_df %>% setNames(c('Origin', 'Origin_celltype'))
    neighbor_celltype_count_df = left_join(neighbor_celltype_count_df, meta_cell_origin_df, by = 'Origin') %>% 
        select(starts_with('Origin'), everything())

    return(neighbor_celltype_count_df)
}

### Reorder and Cluster
# Step 2 Reorder Neighborhood df
ClusterAndReorderTable = function(neighbor_count_df, n_km=3, celltype_col = 'Neighbor_celltype'){
    ########################
    ### Reorder and Cluster
    ########################
    # Covert to matrix form for hclust
    celltype_mtx = list()
    celltype_mtx$count   = neighbor_count_df %>% pivot_wider(id_cols = Origin, names_from = .data[[celltype_col]], 
                                                                       values_from = n, values_fill = 0) %>% column_to_rownames('Origin')
    celltype_mtx$percent = neighbor_count_df %>% pivot_wider(id_cols = Origin, names_from = .data[[celltype_col]], 
                                                                       values_from = percent, values_fill = 0) %>% column_to_rownames('Origin')


    # Get hclust order - takes a half minute 
    # currently use count - could use percent next 
    cell_order_df = GetHclustKmeansTable(celltype_mtx$count, n_km)

    # Add order to LONG count df
    if(any(c('hclust','kmeans') %in% names(neighbor_count_df))){
        message('Found hclust kmeans column in the neighborhood df, removing them')
        neighbor_count_df = neighbor_count_df %>% select(!c('hc_order', 'kmeans'))
    }
    meta_celltype_count_order_df = left_join(neighbor_count_df, cell_order_df, by='Origin')

    # Order origin - this take long for some reason
    # Might try just arrange next time
    message('Reordering sample..')
    meta_celltype_count_order_df = meta_celltype_count_order_df %>% 
        mutate(Origin = factor(Origin, levels =  cell_order_df$Origin[cell_order_df$hc_order]))
 
    message('Done!')
    return(meta_celltype_count_order_df)
}



## Step 3 : 
## Make Neighborhood cell type profile plot
MakeNeighborProfilePlot = function(neighbor_df, palette_use='', Origin_celltype = '', title_sample = '', palette_cluster = ''){
    if(!all(c('hc_order', 'kmeans') %in% names(neighbor_df))) stop('Missing hclust and keamns. Run ClusterAndReorderTable first')
    palette_use = ExtendColorVector(palette_use, target_class = neighbor_df[['Neighbor_celltype']])

    # Profile matrix plt
    p_profilemtx = neighbor_df  %>% 
    ggplot(aes(x = Origin, y = percent, fill = Neighbor_celltype)) + 
        facet_grid(.~kmeans, scales='free', space = 'free')+
        geom_bar(stat = 'identity') + 
        scale_fill_manual(values = palette_use) + 
        cowplot::theme_cowplot()+
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),)
    
    # KM count Plt
    ## check how many cell in each km group
    cluster_count_df = neighbor_df[,c('Origin','kmeans')] %>% ungroup %>% distinct %>% count(kmeans) %>% 
        mutate(kmeans = as.character(kmeans),
            CellType = Origin_celltype,
              Percent = n /sum(n),
              Cell_cluster = paste0(all_of(Origin_celltype),'- clst',kmeans),
              Label = paste0(Cell_cluster, '\n', n, ' (', round(Percent,3)*100, '%)'))
    
    palette_cluster = ExtendColorVector(palette_cluster, target_class = cluster_count_df[['Cell_cluster']]) # Palette 
    
    ## plot
    p_kmpercent = cluster_count_df %>% ggplot(aes(x = CellType, y = n, fill = Cell_cluster)) + 
        geom_bar(stat ='identity') + 
        geom_text(aes(label = Label), position = position_stack(vjust = 0.5)) +
        scale_fill_manual(values = palette_cluster) +
        cowplot::theme_cowplot()
    
    # Plot Both
    wrap_plots(m = p_profilemtx, k = p_kmpercent, design = 'mmmmmmk', guides = 'collect') + 
        plot_annotation(title = str_glue('{title_sample} - {Origin_celltype}: Neighbor Cell Profile'),
                       theme = theme(plot.title = element_text(face = 'bold', size = 20)))
    
}


########################
### Internal 
### Reorder and Cluster
########################
## Get hclust and keam table based on a mtx
GetHclustKmeansTable = function(mtx, n_km, cluster_method = c('kmeans','hclust')){
    # Get hclust order - takes a half minute 
    message('Calculating hclust..')
    tmp_hc = hclust(dist(mtx))
    
    # Cluster
    cluster_method = match.arg(cluster_method)
    if(cluster_method == 'kmeans'){
        # Get kmean groups
        message(str_glue('Calculating kmeans, using km = {n_km} ..'))
        tmp_km = kmeans(mtx, centers = n_km)
        cluster_lab = tmp_km$cluster
    }else{
        message(str_glue('Using cutree and hclust for clustering, using k = {n_km} ..'))
        cluster_lab = cutree(tmp_hc, k = n_km)
    }


    # Turn hc and km into table
    cell_order_df = data.frame(Origin = names(tmp_km$cluster),
                                   hc_order = tmp_hc$order,
                                   kmeans = cluster_lab  # need to change label for this column later
                                  )
    return(cell_order_df)
}


MakeClusterPalette = function(k, celltype=''){
    target_class = paste0(all_of(celltype),'- clst',seq_len(k))
    ExtendColorVector(color_vector = '', target_class = target_class)
}