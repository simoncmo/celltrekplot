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

RemoveCelltrekDuplication <- function(celltrek_obj) {
  cell_keeps = celltrek_obj@meta.data %>% filter(id_raw == id_new) %>% pull(id_raw)
  subset(celltrek_obj, cells = cell_keeps)
}