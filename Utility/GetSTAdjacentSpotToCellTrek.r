## Get cell close to st spot
# Currently take a couple minutes.
GetSTAdjacentCellTrek = function(celltrek_obj, ST, ST_spots, multithread = F){
    mapfun = if(multithread){ furrr::future_pmap }else{pmap}
    # Get ST Coordinates 
    st_cyst_coord = ST@images[[1]]@coordinates[,c('imagerow','imagecol')] %>% 
        .[ST_spots, ] %>%
        setNames(c('st_x','st_y')) %>% 
        rownames_to_column('ST_id')
    
    # Do calculation : 
    # Compare each CellTrek cell distance 
    mapfun(celltrek_obj@meta.data, function(coord_x, coord_y, id_new, ...){
        mapfun(st_cyst_coord, function(st_x, st_y, ST_id, ...){
            distance = sqrt((coord_x - st_x)^2 + (coord_y - st_y)^2)
            #if(distance < dist_threshold){
                # Found dist 
                data.frame(
                    # CellTrek
                    id_new = id_new,
                    celltrek_x = coord_x,
                    celltrek_y = coord_y,
                    # ST
                    st_x = st_x,
                    st_y = st_y,
                    ST_id = ST_id,
                    # distance
                    distance = distance
                )
            #    break
            #}
            
        }) %>% bind_rows() 
    }) %>% bind_rows() 
    
}