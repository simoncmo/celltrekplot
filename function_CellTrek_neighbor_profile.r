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
