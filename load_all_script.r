# Read and load all scripts
library(tidyverse)
this_script <- '/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/CellTrek/shared/script/load_all_script.r'# getSrcDirectory(function(x) {x}) - this failed .... manually set path
path_of_this_script = dirname(this_script) 
name_of_this_script = basename(this_script)

# Get all scripts 
all_scripts = list.files(path_of_this_script, recursive = T, pattern = '*.r', full.names = T) %>%
    str_subset(name_of_this_script, negate = T) # Remove this own script from the set 

# Source
walk(all_scripts, source)