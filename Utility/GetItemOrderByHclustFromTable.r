# Get column item order using hclust
GetItemOrderByHclustFromTable = function(table, sample_column, feature_column){
    table %>% 
    count(.data[[sample_column]], .data[[feature_column]]) %>% 
    pivot_wider(id_cols    = .data[[sample_column]], 
                names_from = .data[[feature_column]], 
                values_from = n, 
                values_fill = 0) %>% 
    column_to_rownames(sample_column) %>% 
    dist() %>% hclust() %>% {.$labels[.$order]}
}
