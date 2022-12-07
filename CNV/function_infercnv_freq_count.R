####################################
# Count Frequency of CNV events from Average cnv mtx
## DEMO:
## avg_cnv_freq_mtx = GetCNVFrequencyMatrix(avg_cnv_mtx)
####################################
# Get cluster level CNV frequency
GetCNVFrequencyMatrix = function(avg_cnv_mtx, remove_cnv_one = T, aggr_level = c('Sample','Cluster')){
    ## 0. Select which level to aggregate at
    aggr_level = match.arg(aggr_level)
    avg_cnv_mtx = if(aggr_level == 'Sample'){
        avg_cnv_mtx %>% mutate(Sample = str_remove(rownames(.), '_CNV_[0-9]')) # Sample level
    }else{ # Cluster level
        avg_cnv_mtx %>% mutate(Sample = rownames(.)) # Cluster level
    }
    ## Unique 'sample'
    n_sample = avg_cnv_mtx$Sample %>% unique %>% length
    
    ## 1. Count CNV event per gene
    avg_cnv_freq_mtx = avg_cnv_mtx %>%  
        group_by(Sample) %>% 
        pivot_longer(cols = -Sample, names_to = 'Gene', values_to = 'CNV') %>% 
        distinct %>% 
        # Count event 
        ungroup() %>%
        count(Gene, CNV) %>% 
        pivot_wider(id_cols = Gene, names_from = CNV, values_from = n) 
    
    ## 2. Filter out CNV = 1
    if(remove_cnv_one) avg_cnv_freq_mtx = avg_cnv_freq_mtx %>% select(-`1`)
    
    ## 3. Summarize frequency and sort
    avg_cnv_freq_mtx %>% 
        mutate(
            sample_size = n_sample, 
            amp_events = rowSums(.[,str_subset(names(.),'1.5|2|3'),drop=F], na.rm=T),
            del_events = rowSums(.[,str_subset(names(.),'0'),drop=F], na.rm=T),
            total_events = amp_events + del_events ) %>% 
        arrange(desc(total_events))
    
}
