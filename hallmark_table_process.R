library(tidyverse)

# Process table
halmark_raw = read_tsv('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/CellTrek/shared/table/halmark_Eduard_Mustafa.tsv')
halmark_long = halmark_raw %>% pivot_longer(cols = everything(), names_to = 'pathway', values_to = 'gene') %>%
    filter(!is.na(gene))
halmark_list = split(x = halmark_long$gene, f = halmark_long$pathway)
write_tsv(halmark_long, '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/CellTrek/shared/table/halmark_Eduard_Mustafa_long.tsv')