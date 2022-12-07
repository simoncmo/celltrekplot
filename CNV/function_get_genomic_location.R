################################################
# METHOD 1: from Biomart 
## Usage:
## chr_loc_df = GetChromeTableBioMart(gene_use)
################################################
# downside: No Cytogenic band location
## https://www.notion.so/Gene-to-Genomic-location-8255ff4de4814b4ba1f9a4f01c636d73#38a8683b6f3c41ba8b73b8b8826cb1f3
GetChromeTableBioMart = function(gene_list){
    library(biomaRt)
    grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="https://grch37.ensembl.org",
    path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    ensembl = useDataset("hsapiens_gene_ensembl",mart=grch37)
    t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",
    "external_gene_name",'chromosome_name','start_position','end_position'), mart = ensembl)
    gene1<-gene_list  #c("GENE1", "GENE2") #vector of your gene names
    gene1<-as.data.frame(gene1)
    colnames(gene1)[1]<-"external_gene_name"
    my_ids.version <- merge(gene1, t2g, by= 'external_gene_name')
    return(my_ids.version)
}


################################################
# Method 2. From Organism dplyr and another package
# Downside : Need extra installation and loading of the data.base
# Usage: chr_cyto_df = GetCytogenicTableUCSC(gene_use)
################################################
GetCytogenicTableUCSC = function(gene_use){
   library(Organism.dplyr)
    
   # Check if data loaded. If not run it!
   check_var = ls(envir = .GlobalEnv, pattern = 'src_human')
    if(length(check_var) == 0){
        print('src_human not found. Creating one!')
        src_human <<- src_ucsc("Human") # Super assign to .GlobalEnv
    }
    
   # https://support.bioconductor.org/p/100890/
    df = tbl(src_human, "id") %>% 
        filter(symbol %in% gene_use) %>% 
        dplyr::select(map, symbol) %>% 
        distinct() %>% as.data.frame 
    out = data.frame(Gene = gene_use,
                     Cytogenic = column_to_rownames(df, 'symbol') %>% .[gene_use, ]
                    ) %>% mutate(Cytogenic_short = str_extract(Cytogenic, '[0-9]{1,2}[p,q]'))
    return(out)
}
