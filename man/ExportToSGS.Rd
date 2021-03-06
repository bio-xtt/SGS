% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sgsR.R
\name{ExportToSGS}
\alias{ExportToSGS}
\title{ExportToSGS}
\usage{
ExportToSGS(
  object,
  species_id,
  track_name,
  track_type,
  select_group,
  reductions = NULL,
  markers.file = NULL,
  ident.field = NULL,
  matrix.slot = "data",
  markers.n = 100
)
}
\arguments{
\item{object}{Seurat object}

\item{species_id}{the id of the species}

\item{track_name}{the name of the single cell track}

\item{track_type}{the type of the single cell track, transcript or atac}

\item{select_group}{vector of cell group names to export.}

\item{reductions}{vector of reduction names to export, defaults to all reductions.}

\item{markers.file}{path to file with marker genes.}

\item{ident.field}{name of the idents used to caculate the marker gene when there is no marker file}

\item{matrix.slot}{matrix to use ("data" or "counts"), default is 'scale.data'.}

\item{markers.n}{no marker files were supplied, FindAllMarkers is run, markers.n is used to setting the number of gene to export.}
}
\value{
This function exports Seurat object as a set of merged matrix feather file and the json format information to send post
}
\description{
this function used to add Seurat object to SGS browser
}
\examples{
\dontrun{

data("pbmc_small")

marker_file <- system.file("extdata", "markers.tsv", package = "SGSR")

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
}

}
\author{
xtt
}
