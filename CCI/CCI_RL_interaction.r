### Cell Cell Interactino
## R-L interaction
# 7/2/2022


# [API] 1. CCI - first get score for RL per each dot connected by edges
GetReceptorLigandDistScoreTable = function(celltrek_obj, plot_info_obj, receptor, ligand){
    # parameter
    edge_df = plot_info_obj$Edge_table
    # Get Expression
    exp_df = FetchData(celltrek_obj, vars = c(receptor, ligand)) %>% 
        setNames(c('Receptor','Ligand')) %>% 
        rownames_to_column('cell_id')
    # Add to Edge table
    edge_exp_df = merge(edge_df, 
          exp_df %>% setNames(str_c('cell1_', names(.))), 
          by.x = 'cell1', by.y = 'cell1_cell_id', all.x=T) %>% 
    merge(exp_df %>% setNames(str_c('cell2_', names(.))), 
          by.x = 'cell2', by.y = 'cell2_cell_id', all.x=T)
    # Add scores - Cell1 R exp * Cell2 L exp OR Cell2 R exp * Cell1 L exp OR
    edge_exp_df = edge_exp_df %>% mutate(C1C2RL_score = cell1_Receptor * cell2_Ligand,
                                         C2C1RL_score = cell2_Receptor * cell1_Ligand)
    edge_exp_df
}

# [API] 2. Get scores from RL table 
GetInteractionScores = function(RLscore_table, group.by = c('group','celltype')){
    group.by = match.arg(group.by)
    colname_use = str_c('Edge_', group.by)
    # Get score
    RLscore_table %>% 
        group_by(.data[[colname_use]]) %>% 
        summarize(C1C2RL_score = mean(C1C2RL_score), 
                  C2C1RL_score = mean(C2C1RL_score))
}

# [API] 3. This get the 1. DistScore table, 2. score per group, 3. score per cell types
GetRLscores = function(celltrek_obj, plot_info_obj, receptor, ligand){
    # Get tables
    RLscore_table = GetReceptorLigandDistScoreTable(celltrek_obj, plot_info_obj, receptor, ligand)
    RLscore_outs = list(RLscore_table = RLscore_table, 
        RLscore_cell_type = GetInteractionScores(RLscore_table, group.by = 'celltype'),
        RLscore_group = GetInteractionScores(RLscore_table, group.by = 'group'))

    # Make obj
    CCI_list = list(dummy_name =RLscore_outs) %>% setNames(str_glue('{receptor}_{ligand}'))
    
    # Check if already have CCI and if so integrate, else initiate
    CCI_list = if(!('CCI' %in% names(plot_info_obj) )){
        list()
    }else{
        plot_info_obj$CCI
    }
    #return(CCI_list)
    # Add 
    CCI_list[[str_glue('{receptor}_{ligand}')]] = RLscore_outs
    #return(CCI_list)
    # Add to PLOT object
    plot_info_obj$CCI = CCI_list
    plot_info_obj
}

## [API] 4. Make CellTrek RL plot
MakeCellTrekRLscore = function(plot_obj, pair_name = NA, direction = c('C1C2','C2C1')){
    if(!'CCI' %in% names(plot_obj)) stop('No CCI cound. Use GetRLscores first')
    if(!pair_name %in% names(plot_obj$CCI)){
        warning(str_glue('{pair_name} not found in CCI slot. Use first pair {names(plot_obj$CCI)[[1]]} instead '))
        pair_name = names(plot_obj$CCI)[[1]]
    }
    # Set score direction
    direction = match.arg(direction)
    direction_column = str_c(direction, 'RL_score')
    
    # expression range 
    score_table = plot_obj$CCI[[pair_name]]$RLscore_table
    exp_max = max(score_table[[direction_column]])
    # Plot
    p = ggplot()

    p = p + geom_segment(data = score_table, 
                         aes(x = x1, y = y1, xend = x2, yend = y2, 
                             color = .data[[direction_column]] ))  
                         #size = segment_size)
    p + scale_color_gradientn(colours = c('gray90','orange','purple'), limits = c(0,0.9*exp_max))+ 
            #scale_color_manual(values = palette_use) + 
            #labs(title = plt_title) + 
            theme_void()+
            theme(aspect.ratio = 1, plot.title = element_text(size = 3, face = "bold", hjust = 0.5))
} 