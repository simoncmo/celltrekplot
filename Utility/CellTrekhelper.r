################################################
# General helper
################################################
#' Subset Celltrek to keep cell from scRNA-seq object only
#'
#' This function keep cell from scRNA-seq only. 
#' This function uses id_raw == id_new in the meta.data to filter cells
#' @param celltrek_obj CellTrek result
#' @return CellTrek result without any duplicated cells
#' @export
#' @examples
#' 



RenameSeuratCells <- function(obj){
    Seurat::RenameCells(object = obj, new.names = Cells(obj) %>% make.names)
}


################################################################################
## Remove cells
################################################################################
## OLD VERSION
RemoveCellNotInCellTrek = function(...){
    warning('Old func. Use RemoveScRNANotInCellTrek instead')
    RemoveScRNANotInCellTrek(...)
}

# For CellTrek
RemoveCelltrekDuplication <- function(celltrek_obj) { # Was "CelltrekRemoveDuplication"
  cell_keeps = celltrek_obj@meta.data %>% filter(id_raw == id_new) %>% pull(id_raw)
  subset(celltrek_obj, cells = cell_keeps)
}

# For scRNA
RemoveScRNANotInCellTrek = function(scRNA_obj, CellTrek_obj){
    scRNA_obj = RenameSeuratCells(scRNA_obj) # Make sure renamed
    celltrek_used_cells = CellTrek_obj@meta.data$id_raw %>% unique # Extract cells
    scRNA_obj_new = subset(scRNA_obj, cells = celltrek_used_cells) # Subset
    message(str_glue("Kept {ncol(scRNA_obj_new)} cells in the new scRNA object"))
    return(scRNA_obj_new)
}
################################################################################
################################################################################

## This function Add SCT from sc/snRNA to cellTrek
AddAssayToCellTrek = function (celltrek_obj, scRNA, assay = "SCT") 
{
    sc_new = RenameCells(scRNA, new.names = make.names(colnames(scRNA)))
    # create celltrek id without tailing .1 to .9
    celltrek_id_new_notail = colnames(celltrek_obj) %>% str_remove('\\.[0-9]$')
    celltrek_meta = celltrek_obj@meta.data
    celltrek_assay = sc_new@assays[[assay]]@data %>% .[, celltrek_id_new_notail] # this will allow same cell name to duplicate
    colnames(celltrek_assay) = celltrek_meta$id_new # rename to celltrek id
    celltrek_obj[[assay]] = CreateAssayObject(celltrek_assay)
    DefaultAssay(celltrek_obj) = assay
    celltrek_obj
}

## This function Add Metadata from sn to CellTrek
AddSnMetaToCellTrek = function(celltrek_obj, sn, ident){
    if(ident %in% colnames(celltrek_obj@meta.data)){
        warning(str_glue('{ident} found in celltrek! Please double check to prevent duplication'))
        warning('Return original object')
        return(celltrek_obj)
    }
    
    sn_table = sn@meta.data %>% {data.frame(id_raw = make.names(rownames(.)),
                                           col = .[[ident]]
                                          )} %>% setNames(c('id_raw', ident))
    row_celltrekmeta = rownames(celltrek_obj@meta.data)
    celltrek_obj@meta.data = left_join(celltrek_obj@meta.data, 
             sn_table, 
              by = 'id_raw'
             ) %>% `rownames<-`(row_celltrekmeta)
    return(celltrek_obj)
}