### Add Meta to snRNA or CellTrek
AddSTMetaToSNByDistance = function(distance_table, ST, sn_obj, 
                                   dist_max = 100, 
                                   id_col = 'id_new',
                                   ST_column_keep = 'histology'){
    # For SnOnly
    distance_table  = distance_table %>% mutate(id_new = id_new %>% str_remove('\\.[1-9]$') %>% str_replace('\\.','-'))
    
    # Run
    AddSTMetaToObjByDistance(distance_table, ST, sn_obj, dist_max, id_col = id_col, ST_column_keep, data_type = 'single cell')
}

### Add Meta to snRNA or CellTrek
AddSTMetaToCellTrekByDistance = function(distance_table, ST, sn_obj, 
                                         dist_max = 100, 
                                         id_col = 'id_new',
                                         ST_column_keep = 'histology'){
    # Run
    AddSTMetaToObjByDistance(distance_table, ST, sn_obj, dist_max, id_col = id_col, ST_column_keep, data_type = 'celltrek')
}

# General function
AddSTMetaToObjByDistance = function(distance_table, ST, obj, dist_max = 100, id_col = 'id_new',
                                    ST_column_keep = 'histology', data_type = c('single cell','celltrek')){
    # Filter by distance
    distance_table  = distance_table %>% filter(distance < dist_max)
    
    # Add meta of ST to distance table
    ST_dist_meta = merge(distance_table, 
                         ST@meta.data %>% rownames_to_column('ST_id'), 
                         by = 'ST_id') 

    # Make unique
    ST_dist_meta = ST_dist_meta %>% 
        select(c(all_of(id_col),'distance', all_of(ST_column_keep))) %>% 
        distinct %>%
        group_by(.data[[id_col]])%>% 
        summarise(
            avg_distance_to_closet_ST = mean(distance),
            {{ST_column_keep}} := toString(sort(unique(.data[[ST_column_keep]])))) %>% 
        ungroup %>%
        column_to_rownames(id_col)

    meta_add = ST_dist_meta[,c('avg_distance_to_closet_ST', ST_column_keep)]
    
    # Add To celltrek
    sn_obj = AddMetaData(obj, metadata = meta_add)
    sn_obj
}
