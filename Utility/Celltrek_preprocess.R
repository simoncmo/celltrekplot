# Moved from 'function_CellTrek_wrapper.R'
############################################################################################
######### Preprocess ############################################################
############################################################################################      
## Need to make sure tow metadata match (espeiclaly same data type for column with same name) since it will bind row from both meta
CheckSharedColumns = function(sn, ST){
    shared_cols = intersect(sn@meta.data %>% names, ST@meta.data %>% names)
    
    walk(shared_cols, function(col){
    is_same = sn@meta.data[,col] %>% class  ==     ST@meta.data[,col] %>% class
    if(!is_same){
        warning(str_glue('STOP!! Found column type not matching: {col}. Please fix before continue'))
        print(str_c(col, ', sn = ',sn@meta.data[,col] %>% class, ', ST = ', ST@meta.data[,col] %>% class,', Same = ', is_same))
    }else{
        print(str_glue('{col}: Same'))
    }
    })
}
FixColumnToIntegar = function(obj, columns){
    for(col in col_to_fix){
        obj@meta.data[, col] = obj@meta.data[, col] %>% as.character %>% as.integer
    }
    return(obj)
}
