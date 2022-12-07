source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/CellTrek/shared/script/Utility/GetItemOrderByHclustFromTable.r')

MakeOrderedBarPlot = function(table, sample_column, feature_column, order_by = c('descend','ascend','hclust')){
    # count data
    count_df = table %>% count(.data[[sample_column]], .data[[feature_column]])
    # Get sample order 
    order_by = match.arg(order_by)
    # Sum table
    sum_df = count_df %>% group_by(.data[[sample_column]]) %>% summarize(total = sum(n))
    sample_order = if(order_by=='ascend'){
        sum_df %>% arrange(total) %>% pull(.data[[sample_column]])
    }else if(order_by =='descend'){
        sum_df %>% arrange(desc(total)) %>% pull(.data[[sample_column]])
    }else{
        GetItemOrderByHclustFromTable(table, sample_column, feature_column)
        }
    count_df = count_df %>% mutate({{ sample_column }} := factor(.data[[sample_column]], levels = sample_order))
    # plot
    count_df %>% 
        ggplot(aes(x = .data[[sample_column]], y = .data[['n']], fill = .data[[feature_column]])) +
        geom_bar(stat = 'identity') + 
        cowplot::theme_cowplot() + 
        RotatedAxis()
}