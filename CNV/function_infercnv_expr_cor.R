################################
# Find Feature that have similar cnv - expression pattern
# AKA expression affectly mainly by CNV
# Tumor only
################################
GetCNVExpressionCorrelation = function(sn, cells_include){
    if(!'cnv_with_event' %in% Assays(sn)) stop('Please Add CNV assay cnv_with_event first')
    # Get data
    cnv_mtx = sn$cnv_with_event %>% GetAssayData %>% as.matrix
    feature_mtx = sn$SCT %>% GetAssayData %>% as.matrix
    
    # Extract shared featurs
    feature_shared = intersect(rownames(feature_mtx), rownames(cnv_mtx))

    # Extract shared 'tumor cells' from cnv assay
    tumor_cells = cnv_mtx %>% is.na %>% apply(2, all) %>% .[!.] %>% names# cell with Not ALL NA values
    if(!missing(cells_include)) tumor_cells = intersect(tumor_cells, cells_include)
    
    # Subset both tables
    feature_mtx_tumor = feature_mtx[feature_shared, tumor_cells]
    cnv_mtx_tumor     = cnv_mtx[feature_shared, tumor_cells]
    
    # Correlation test
    map(feature_shared, function(feature){
        cor_result = cor.test(feature_mtx_tumor[feature, ], 
             cnv_mtx_tumor[feature, ])
        data.frame(
            gene = feature,
            cor  = cor_result$estimate,
           p_val = cor_result$p.value) %>% remove_rownames()
    }) %>% bind_rows() %>% 
            mutate(p_val_adj = p.adjust(p = p_val, method = 'fdr')) %>% arrange(desc(cor))
}


################################
## Plot function

################################
# CNV palette
## Color palette for CNV level 0 to 3 by 0.5
# Needed for "FeaturePlotExpCnv"
cnv_palette = c(RColorBrewer::brewer.pal(9, 'PuBu')[c(6,9)] %>% rev, 
           'gray85',
            RColorBrewer::brewer.pal(9, 'YlOrRd')[c(4,6,7,9)]
           ) %>% setNames(seq(0, 3, by = 0.5))


# Plot function
# Need "cnv_palette"
FeaturePlotExpCnv = function(sn, gene, cnv_key = 'cnvwithevent', expr_key = 'sct', shink_cnv_scale = T, ident_col, label_ident = T){
        p0 = FeaturePlot(sn, features = str_glue('{expr_key}_{gene}')) + 
            scale_color_gradientn(colors = c('gray80', RColorBrewer::brewer.pal(9, 'YlOrRd')[3:9])) + 
            labs(color = 'snRNA\nExpr', title = gene, subtitle = str_glue("Expression ({expr_key})"))
        # CNV
        p1 = FeaturePlot(sn, features = str_glue('{cnv_key}_{gene}')) + 
            labs(color = 'CNV', title = gene, subtitle = str_glue("CNV ({cnv_key})"))

        # Find CNV range 
        cnv_range = p1$data[str_glue('{cnv_key}_{make.names(gene)}')] %>% 
                range(., na.rm =T) %>% # Get min max of expression range
                {seq(from = .[[1]], to =.[[2]] , by= 0.5)} #%>% # so that . doesn't become first argument
        cnv_range_chr = as.character(cnv_range) # So that the palette index will be right
            
        # Shrink palette 
        p1 = if(shink_cnv_scale){
            p1 + scale_color_gradientn(colors = cnv_palette[cnv_range_chr], breaks = cnv_range)    
        }else{
            p1 + scale_color_gradientn(colors = cnv_palette, limits = c(0,3), breaks = as.numeric(names(cnv_palette)))
        }
    
        # final 
        p_final = if(!missing(ident_col)){
            p_dim = DimPlot(sn, group.by = ident_col, label = label_ident) + 
                    labs(title = gene, subtitle = ident_col)
             p0 | p1 | p_dim
        }else{
            (p0 | p1)
        }
     p_final & coord_fixed()
}


################################################################
# CellTrek Version - Consider consolidate this with to "FeaturePlotExpCnv" into same internal function
################################################################
# Plot function
# Need "cnv_palette"
CellTrekFeaturePlotExpCnv = function(celltrek_obj, gene, cnv_key = 'cnvwithevent', expr_key = 'sct', shink_cnv_scale = T, ST_obj){
        # Make plots
        p_sct  = CellTrekFeaturePlot(celltrek_obj = celltrek_obj, feature = c(str_glue('{expr_key}_{gene}'))) + 
                    scale_color_gradientn(colors = c('gray80', RColorBrewer::brewer.pal(9, 'YlOrRd')[3:9])) + 
                    labs(color = 'snRNA\nExpr', title = gene, subtitle = str_glue("Expression ({expr_key})"))
        p_cnv  = CellTrekFeaturePlot(celltrek_obj = celltrek_obj, feature = c(str_glue('{cnv_key}_{gene}'))) +
                    labs(color = 'CNV',  title = gene, subtitle = str_glue("CNV ({cnv_key})"))

        # Find CNV range 
        cnv_range = FetchData(celltrek_obj, vars = c(str_glue('{cnv_key}_{gene}')) ) %>% # Not sure why but need c() to work
                range(., na.rm =T) %>% # Get min max of expression range
                {seq(from = .[[1]], to =.[[2]] , by= 0.5)} #%>% # so that . doesn't become first argument
        cnv_range_chr = as.character(cnv_range) # So that the palette index will be right
            
        # Shrink palette 
        p_cnv = if(shink_cnv_scale){
            p_cnv + scale_color_gradientn(colors = cnv_palette[cnv_range_chr], breaks = cnv_range)    
        }else{
            p_cnv + scale_color_gradientn(colors = cnv_palette, limits = c(0,3), breaks = as.numeric(names(cnv_palette)))
        }
        # final 
        p_final = if(!missing(ST_obj)){
            p_ST = SpatialFeaturePlot(ST_obj, gene, image.alpha=0, stroke = NA) + 
                scale_fill_gradientn(colors = c('gray80', RColorBrewer::brewer.pal(9, 'YlOrRd')[3:9])) + 
                labs(title = gene, subtitle = 'Visium ST', fill = 'ST\nExpr')
             p_sct | p_cnv | p_ST
        }else{
            (p_sct | p_cnv)
        }
     p_final & coord_fixed() & cowplot::theme_cowplot() & 
            labs(x = '', y = '') &
            theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())
}