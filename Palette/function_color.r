#' Get discrete color using Seurat Discrete color function
#'
#' This function Get n discrete color using All palettes from Seurat's package Discrete color function
#' @param n The number of colors needed.
#' @return A vector of colors.
#' @export
#' @examples
#' get_discrete_colors(32)
get_discrete_colors = function(n){
    palettes_discrete_list <- list(alphabet = c("#F0A0FF", "#0075DC", "#993F00", 
        "#4C005C", "#191919", "#005C31", "#2BCE48", "#FFCC99", 
        "#808080", "#94FFB5", "#8F7C00", "#9DCC00", "#C20088", 
        "#003380", "#FFA405", "#FFA8BB", "#426600", "#FF0010", 
        "#5EF1F2", "#00998F", "#E0FF66", "#740AFF", "#990000", 
        "#FFFF80", "#FFE100", "#FF5005"), alphabet2 = c("#AA0DFE", 
        "#3283FE", "#85660D", "#782AB6", "#565656", "#1C8356", 
        "#16FF32", "#F7E1A0", "#E2E2E2", "#1CBE4F", "#C4451C", 
        "#DEA0FD", "#FE00FA", "#325A9B", "#FEAF16", "#F8A19F", 
        "#90AD1C", "#F6222E", "#1CFFCE", "#2ED9FF", "#B10DA1", 
        "#C075A6", "#FC1CBF", "#B00068", "#FBE426", "#FA0087"), 
        glasbey = c("#0000FF", "#FF0000", "#00FF00", "#000033", 
            "#FF00B6", "#005300", "#FFD300", "#009FFF", "#9A4D42", 
            "#00FFBE", "#783FC1", "#1F9698", "#FFACFD", "#B1CC71", 
            "#F1085C", "#FE8F42", "#DD00FF", "#201A01", "#720055", 
            "#766C95", "#02AD24", "#C8FF00", "#886C00", "#FFB79F", 
            "#858567", "#A10300", "#14F9FF", "#00479E", "#DC5E93", 
            "#93D4FF", "#004CFF", "#F2F318"), polychrome = c("#5A5156", 
            "#E4E1E3", "#F6222E", "#FE00FA", "#16FF32", "#3283FE", 
            "#FEAF16", "#B00068", "#1CFFCE", "#90AD1C", "#2ED9FF", 
            "#DEA0FD", "#AA0DFE", "#F8A19F", "#325A9B", "#C4451C", 
            "#1C8356", "#85660D", "#B10DA1", "#FBE426", "#1CBE4F", 
            "#FA0087", "#FC1CBF", "#F7E1A0", "#C075A6", "#782AB6", 
            "#AAF400", "#BDCDFF", "#822E1C", "#B5EFB5", "#7ED7D1", 
            "#1C7F93", "#D85FF7", "#683B79", "#66B0FF", "#3B00FB"), 
        stepped = c("#990F26", "#B33E52", "#CC7A88", "#E6B8BF", 
            "#99600F", "#B3823E", "#CCAA7A", "#E6D2B8", "#54990F", 
            "#78B33E", "#A3CC7A", "#CFE6B8", "#0F8299", "#3E9FB3", 
            "#7ABECC", "#B8DEE6", "#3D0F99", "#653EB3", "#967ACC", 
            "#C7B8E6", "#333333", "#666666", "#999999", "#CCCCCC"))
    # Get color
    palettes_discrete_list %>% unlist %>% unique %>% .[1:n]
}