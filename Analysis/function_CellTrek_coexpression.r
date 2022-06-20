library(viridis) ## Used in CoEpxressionPlot

############################################################################################
######### Preprocess ############################################################
############################################################################################      
## Need to make sure tow metadata match (espeiclaly same data type for column with same name) since it will bind row from both meta
CheckSharedColumns = function(sn, ST){
    shared_cols = intersect(sn@meta.data %>% names, ST@meta.data %>% names)
    
    walk(shared_cols, function(col){
    is_same = sn@meta.data[,col] %>% class  ==     ST@meta.data[,col] %>% class
    if(!is_same){
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

### NEED MANUAL Assign Cluster: Workflow like below
# CheckSharedColumns(sn, ST)
# col_to_fix = c('SCT_snn_res.0.5', 'seurat_clusters','nCount_SCT')
# ST = FixColumnToIntegar(ST, col_to_fix)
# print('After fixed ..')
# CheckSharedColumns(sn, ST)


############################################################################################
######### CoExpression Analysis ############################################################
############################################################################################
# Turn GS into 
#cc_result = obj_subset_scoexp_res_cc
CoExpressionPlot = function(cc_result){
    if(length(cc_result$gs) == 0 ){
    print("Nothing in gs, cut off might be too stringent. Nothing to plot")
    return()
    }
    
    ccresult_k_df = cc_result$gs  %>% imap(function(genes, group){
        data.frame(row.names = genes, G = rep(group, length(genes)))
    }) %>% bind_rows()

    # mtx to plt
    pheat_mtx = cc_result$wcor[rownames(ccresult_k_df), rownames(ccresult_k_df)]

    # pheatmap
    pheatmap::pheatmap(pheat_mtx, 
                       clustering_method='ward.D2', 
                       annotation_row=ccresult_k_df, 
                       show_rownames=F, show_colnames=F, 
                       treeheight_row=10, treeheight_col=10, annotation_legend = T, fontsize=8,
                       color=viridis(10), main='Spatial co-expression')
}


## Function To Cluster All cells
# Example: 
# cc_result = consensuscluster_multiple_cells(celltrek_obj, sample_name = sample_id, cell_type_col = sn_col, 
#                                 root_path = '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/CellTrek/AKI/')

consensuscluster_multiple_cells = function(obj_celltrek, sample_name = 'sample_1', cell_types, cell_type_col = 'cell_type_Abbr', assay = 'RNA', root_path){
    # path
    if(missing(root_path))  root_path  = getwd()
    cc_result_list = list()
    if(missing(cell_types)) cell_types = obj_celltrek@meta.data[[cell_type_col]] %>% unique 
    
    for(cell in cell_types){
        print(str_glue('Plotting {sample_name}, cell: {cell}'))
        # subset 
        id_celltype = obj_celltrek@meta.data %>% filter( .data[[cell_type_col]] == cell ) %>% rownames
        obj_subset <- subset(obj_celltrek, cells = id_celltype) 
        obj_subset@assays[[assay]]@scale.data <- matrix(NA, 1, 1)

        # We select top 2000 variable genes (exclude mitochondrial, ribosomal and high-zero genes)
        obj_subset <- FindVariableFeatures(obj_subset)
        vst_df <- obj_subset@assays[[assay]]@meta.features %>% data.frame %>% mutate(id=rownames(.))
        nz_test <- apply(as.matrix(obj_subset[[assay]]@data), 1, function(x) mean(x!=0)*100)
        hz_gene <- names(nz_test)[nz_test<20]
        mt_gene <- grep('^Mt-', rownames(obj_subset), value=T)
        rp_gene <- grep('^Rpl|^Rps', rownames(obj_subset), value=T)
        vst_df <- vst_df %>% dplyr::filter(!(id %in% c(mt_gene, rp_gene, hz_gene))) #%>% 
                         #arrange(., -vst.variance.standardized)
        feature_temp <- vst_df$id[1:2000]

        # Create and Move to Folder
        dir.create(str_glue('{root_path}/consensuscluster_result/{sample_name}/{cell}/'), recursive = T)
        setwd(str_glue('{root_path}/consensuscluster_result/{sample_name}/{cell}/'))

        # We use scoexp to do the spatial-weighted gene co-expression analysis.
        obj_subset_scoexp_res_cc <- CellTrek::scoexp(celltrek_inp=obj_subset, 
                                                    assay=assay , approach='cc', 
                                                    gene_select = feature_temp, sigm=140, 
                                                    avg_cor_min=.4, zero_cutoff=3, min_gen=40, max_gen=400)
        # SAVE cc result
        # save RDS with module score
        saveRDS(obj_subset_scoexp_res_cc, str_glue('cluster_result.rds'))
        cc_result_list[[cell]] = list()                         
        cc_result_list[[cell]]$cc_result   = obj_subset_scoexp_res_cc # Save as a list too for output
                         
        ## check if any result to plot
        if(length(obj_subset_scoexp_res_cc$gs)==0){
            print("Nothing in GS. Nothing to Plot") 
            next # Skip All steps below
        }

        # CoExpression Heatmap Plot
        pdf('pheatmap_coexpression.pdf')
        print(CoExpressionPlot(obj_subset_scoexp_res_cc))
        dev.off()
        
        # ModuleScore
        obj_subset <- AddModuleScore(obj_subset, features=obj_subset_scoexp_res_cc$gs, name='CC_', nbin=10, ctrl=50, seed=42)
        pdf('module_score_plots.pdf')
        ## First we look into the coexpression module based on the scRNA-seq embedding
        p1 = FeaturePlot(obj_subset, grep('CC_', colnames(obj_subset@meta.data), value=T)) & theme(aspect.ratio = 1)
        ## Next we investigate the module scores at the spatial level.
        p2 = SpatialFeaturePlot(obj_subset, grep('CC_', colnames(obj_subset@meta.data), value=T)) & theme(aspect.ratio = 1) 
        print(p1)
        print(p2)
        dev.off()
        
        # save RDS with module score
        cc_result_list[[cell]]$seurat_obj  = obj_subset 
                         
        # Return back to root
        setwd(root_path)
    }
    return(cc_result_list)
}
                         
