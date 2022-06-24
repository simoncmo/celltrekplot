#' Plot Degree distribution of a graph obj
#'
#' This function Plot Degree distribution of a graph obj.
#' @param graph_obj A graph obj
#' @return A ggplot histogram for degree of distribution
#' @export
#' @examples
#' PlotDegreeDistribution(graph_obj)
PlotDegreeDistribution = function(graph_obj){
    deg_freq = graph_obj %>% degree_distribution() # From igraph
    # Format
    graph_dist_df = data.frame(degree = 0:(length(deg_freq)-1), freq = deg_freq, count = deg_freq * vcount(graph_obj))
    # plot
    graph_dist_df %>% ggplot(aes(x = degree, y = count)) + 
    geom_bar(stat = 'identity') + 
    scale_x_continuous(breaks = 0:(length(deg_freq)-1)) +
    labs(title = 'Graph QC: Degree Distribution', x = 'degree(# of neighboring cells)', y = 'Cell Count') +
    cowplot::theme_cowplot()
}