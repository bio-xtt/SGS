# SGS
A package used to load the Seurat single cell analysis object into SGS cellbrowser

## Installation
the package can install from github like this:

```
# install from github
devtools::install_github("bio-xtt/SGS") 

# install from source
install.packages("/home/bio-xtt/Desktop/SGSR_0.1.0.tar.gz", type=source)

```


## Usage

```

# load the package
library(SGSR)

# load the example datasets: pbmc_small
data("pbmc_small")

# load the example marker genes: marker.tsv
marker_file <- system.file("extdata", "markers.tsv", package = "SGSR")

# use the export function
result_content <- ExportToSGS(
  object = pbmc_small,
  species_id = "ca4570d9c5464f10915e2def60b94780",
  track_name = "pbmc_small",
  track_type = "transcript",
  reductions = c("tsne", "umap"),
  markers.file = marker_file,
  select_group = c("groups", "RNA_snn_res.0.8","letter.idents","cluster"),
  matrix.slot = "data"
)

```

## Result
the function returns the status of loadding single cell object into SGS cellbrowser:  

1:loadding successfull!  

2:loadding failed  






