# CellTrek Plot
Extra Neighborhood analysis and plotting functions for CellTrek obj


## Delauney and Neighborhood Percentage
- `MakeCellTrekDelaunayGraph` - Delany triangulation and make graph using iGraph
- `GetAllNeighbors` - Get n-hop, n-order neighbor
- `GetNeighborByIdent` - Get n-hop, n-order neighbor using ident in the celltrek obj

`MakeCellNeighborhoodRatioPlot` - Make Neighborhood plot. support bar or matrix
```R
cell_column = 'cell_type_Abbr' # Column for the cell type
celltype_plt = c('PT(dual_identity)','PT(S1+S2)','PT(S3)') # remove this value if want to plot All cell types 
obj_graph = MakeCellTrekDelaunayGraph(celltrek_obj)
MakeCellNeighborhoodRatioPlot(obj = celltrek_obj, 
                         obj_graph = obj_graph,
                         cell_column  = cell_column, 
                         plot_type = 'bar',
                         celltype_plt = celltype_plt) 
```
> **_NOTE:_**  This is order 3 (3 hops) neighborhood cell type proportion 
<img src="https://user-images.githubusercontent.com/54045654/168443340-3a0229b9-b5bd-4857-922c-8b3bfe7807ef.png" alt="drawing" width="450"/>


## Make Delaunay Plot from CellTrek obj
- `MakeDelaunayPlot` - Make CellTrek ST plot and connnect each dot by edges of Delaunay graph
```R
MakeDelaunayPlot(obj = result_list_celltrek$cell_type_Abbr$W12_M, cell_column = cell_column,
                  cell_plot = c('PT(dual_identity)','PT(S3)','PT(S1+S2)'),
                  palette = col_cell_type, vertex.size = 1.5)
```
> **_NOTE:_**  Current version haven't consider maximum neighborhood distance
<img src="https://user-images.githubusercontent.com/54045654/168443265-37276938-98ca-4c5c-916d-7ccf92489877.png" alt="drawing" width="450"/>
