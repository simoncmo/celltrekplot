################################
# Add CNV gene level observation to Seurat obj
################################

# [API] Add to CNV observation to 
AddCNVAssayToSeurat = function(sn, cnv_obs_table, add_filtered_table = T,
                               assay_name = 'cnv_all', filtered_assay_name = 'cnv_with_event'){
    
    message(str_glue('Total genes: {nrow(cnv_obs_table)}'))
    message(str_glue('Adding whole CNV table as {assay_name}'))
    sn[[assay_name]] = CreateCNVAssayObject(sn, cnv_obs_table) 
    
    # Add filtered table 
    if(add_filtered_table){
        # Create filter table
        gene_with_cnv_events = apply((infercnv_obs != 1), 1, any) %>% names(.)[.] 
        message(str_glue('Total filtered genes: {nrow(gene_with_cnv_events)}'))
        infercnv_obs_with_event = infercnv_obs[gene_with_cnv_events, ]
        
        message(str_glue('Adding filtered CNV table as {filtered_assay_name}'))
        sn[[filtered_assay_name]] = CreateCNVAssayObject(sn, infercnv_obs_with_event) 
    }
                
    return(sn)
    
}

# Create obj and fill missing CELLS
# Add 'reference cells' back as NA
CreateCNVAssayObject = function(sn, cnv_obs_table){
    non_tumor_cells = setdiff(colnames(sn), colnames(cnv_obs_table))
    # Create ref table
    features = rownames(cnv_obs_table)
    ref_df = data.frame(matrix(NA, nrow = length(features), 
                      ncol = length(non_tumor_cells)),
               row.names = features 
              ) %>% setNames(non_tumor_cells)
    # Combine and create Assay
    bind_cols(cnv_obs_table, ref_df) %>% 
        .[, colnames(sn)] %>%  # reorder
        Seurat::CreateAssayObject(.)
}

