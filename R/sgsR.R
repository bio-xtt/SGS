
#' Findmatrix type
#' @description find expression matrix from seurat object like:counts,scale.data
#' @param object Seurat object
#' @param matrix.slot the name of the slot
#'
#' @return Return a matrix object from a Seurat object or show an error message
#'
findMatrix <- function(object, matrix.slot) {
  if (matrix.slot == "counts") {
    counts <- GetAssayData(object = object, slot = "counts")
  } else if (matrix.slot == "scale.data") {
    counts <- GetAssayData(object = object, slot = "scale.data")
  } else if (matrix.slot == "data") {
    counts <- GetAssayData(object = object)
  } else {
    stop("matrix.slot can only be one of: counts, scale.data, data")
  }
  return(counts)
}





#' Extract the data from Seurat object
#' Used by ExportsToSGS
#' @description  This function used to  gain the data from the seurat object
#' @param object Seurat object
#' @param reductions vector of reduction names to export, defaults to all reductions.
#' @param markers.file path to file with marker genes.
#' @param ident.field cluster.field name of the metadata field containing cell cluster
#' @param matrix.slot matrix to use, default is 'scaled data'
#' @param markers.n if no markers were supplied, FindAllMarkers is run.
#' This parameter indicates how many markers to calculate, default is 100
#' @return Return the list of merged df and list of informs
#'
#'
gain_data <- function(object,
                      reductions = NULL,
                      markers.file = NULL,
                      ident.field = NULL,
                      matrix.slot = "scale.data",
                      markers.n = 100) {

  ## look for seurat and check the version of object is the same as the version of the package
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("This script requires that Seurat (V2 or V3) is installed")
  }

  message("Seurat Version installed: ", packageVersion("Seurat"))
  message("Object was created with Seurat version ", object@version)

  objMaj <- package_version(object@version)$major
  pkgMaj <- package_version(packageVersion("Seurat"))$major


  if (objMaj != 2 && objMaj != 3 && objMaj != 4) {
    stop("can only process Seurat2, Seurat3  or Seurat4 objects, object was made with Seurat ", object@version)
  }

  if (objMaj != pkgMaj) {
    stop("The installed major version of Seurat is different from Seurat input object. You have to down- or upgrade your installed Seurat version. See the Seurat documentation.")
  }


  # compatibility layer for Seurat 2 vs 3 / 4
  # see https://satijalab.org/seurat/essential_commands.html

  if (inherits(x = object, what = "seurat")) {
    # Seurat v2 objects are called "seurat" (Paul Hoffman)
    # -> Seurat 2 data access
    idents <- object@ident # Idents() in Seurat3
    meta <- object@meta.data
    cellOrder <- object@cell.names

    ## this has changed by xtt
    if (matrix.slot == "counts") {
      counts <- object@raw.data
    } else if (matrix.slot == "scale.data") {
      counts <- object@scale.data
    } else if (matrix.slot == "data") {
      counts <- object@data
    } else {
      stop("matrix.slot can only be one of: counts, scale.data, data")
    }
    dr <- object@dr
  } else {
    # Seurat 3 / 4 functions
    idents <- Idents(object = object)
    meta <- object[[]]
    cellOrder <- colnames(x = object)
    counts <- findMatrix(object = object, matrix.slot = matrix.slot)
    if (dim(x = counts)[1] == 0) {
      message(paste0("The Seurat data slot '", matrix.slot, "' contains no data. Trying default assay."))
      defAssay <- DefaultAssay(object)
      assay <- GetAssay(object, defAssay)
      message(paste0("Default assay is ", defAssay))
      counts <- findMatrix(assay, matrix.slot)
      if (dim(x = counts)[1] == 0) {
        stop(
          "Could not find an expression matrix",
          "Please select the correct slot where the matrix is stored, possible ",
          "values are 'counts', 'scale.data' or 'data'. To select a slot, ",
          "use the option 'matrix.slot' from R or the cbImportSeurat option -s from the command line."
        )
      }
    }
    dr <- object@reductions
  }

  #### set the dir to write the save the data files
  ## this has add by xtt
  if (!requireNamespace("uuid")) install.packages("uuid")
  sc_id <- uuid::UUIDgenerate()

  if (!dir.exists("/home/sgs/data/rstudio_id")) {
    dir.create("/home/sgs/data/rstudio_id")
  }

  ################ set the sc track dir
  dir <- paste0("/home/sgs/data/rstudio_id/", sc_id)
  if (!dir.exists(paths = dir)) {
    dir.create(path = dir)
  }
  if (!dir.exists(paths = dir)) {
    stop("Output directory ", dir, " cannot be created or is a file")
  } else {
    message("the output dir is:", dir)
  }


  # Export expression matrix and deal the matrix for merge
  counts <- as.data.frame(counts)

  counts_t <- t(counts)

  ###### translate the gene name into lower
  colnames(counts_t) <- tolower(colnames(counts_t))
  exp_df_t <- data.frame(cell = rownames(counts_t), counts_t, check.names = FALSE)
  rownames(exp_df_t) <- NULL

  # Export metadata and deal the matrix for merge
  enum.fields <- c()

  meta.fields <- colnames(x = meta)

  df <- data.frame(row.names = cellOrder, check.names = FALSE)

  ### this has changed by xtt

  if (!is.null(ident.field)) {
    df <- cbind(idents, df)
    colnames(df)[1] <- ident.field
  } else {
    df$cluster <- idents
  }

  ## this has changed by xtt
  for (field in meta.fields) {
    name <- field
    df[[name]] <- meta[[field]]
    if (!is.numeric(df[[name]])) {
      enum.fields <- c(enum.fields, name)
    }
  }


  meta_df <- data.frame(cell = rownames(df), df, check.names = FALSE)

  ###translate the datatype into character
  meta_df <- as.data.frame(lapply(meta_df, as.character))

  cell_meta_column <- colnames(df)

  rownames(meta_df) <- NULL

  #### to merger the merge exp_df_t and meta_df into total_df

  total_df <- merge(exp_df_t, meta_df, by = "cell")


  ## Export cell embeddings/reductions for merge
  # Export cell embeddings/reductions
  reducNames <- reductions

  if (is.null(reducNames)) {
    reducNames <- names(dr)
    message("Using all embeddings contained in the Seurat object: ", reducNames)
  }

  for (embedding in reducNames) {
    emb <- dr[[embedding]]
    if (is.null(x = emb)) {
      message("Embedding ", embedding, " does not exist in Seurat object. Skipping. ")
      next
    }
    df <- emb@cell.embeddings
    if (ncol(x = df) > 2) {
      warning("Embedding ", embedding, " has more than 2 coordinates, taking only the first 2")
      df <- df[, 1:2]
    }
    colnames(x = df) <- c(sprintf("%s_x", embedding), sprintf("%s_y", embedding))

    df <- data.frame(cell = rownames(x = df), df, check.names = FALSE)

    rownames(df) <- NULL

    total_df <- merge(total_df, df, by = "cell")
  }

  ## cell_plots used to return plot type
  cell_plots <- reducNames


  ## set path list
  path_list <- list()

  ######## write to father
  if (!requireNamespace("feather", quietly = TRUE)) {
    message("feather has not been install,please install it and run the script!")
    # devtools::install_github("wesm/feather/R")
  } else {

    feather_file_path <- file.path(dir, "all_file.feather")
    message("the feather file path is:", feather_file_path)
    feather::write_feather(total_df, feather_file_path)
    path_list$feather_file <- feather_file_path
  }

  ##### add the expression gene name

  expression_genes_list <- list("expression_genes" = colnames(counts_t))
  exp_gene_path <- file.path(dir, "expression_gene_file.json")
  jsonlite::write_json(expression_genes_list, path = exp_gene_path)

  path_list$expression_gene_file <- exp_gene_path


  ### Export markers

  if (is.null(markers.file)) {
    ext <- "tsv"
  } else {
    ext <- tools::file_ext(markers.file)
  }
  file <- paste0("markers.", ext)
  fname <- file.path(dir, file)
  if (!is.null(markers.file)) {
    message("Copying ", markers.file, " to ", fname)
    file.copy(markers.file, fname)
    path_list$marker_gene_meta_files <- list(fname)
  }
  if (is.null(markers.file)) {
    file <- NULL
  }

  if (is.null(markers.file)) {
    if (length(levels(idents)) > 1) {
      markers.helper <- function(x) {
        partition <- markers[x, ]

        ## seurat v4 has changed
        if (objMaj == 4) {
          ord <- order(partition$p_val_adj < 0.05, -partition$avg_log2FC)
        } else {
          ord <- order(partition$p_val_adj < 0.05, -partition$avg_logFC)
        }
        res <- x[ord]
        naCount <- max(0, length(x) - markers.n)
        res <- c(res[1:markers.n], rep(NA, naCount))
        return(res)
      }

      if (.hasSlot(object, "misc") && !is.null(x = object@misc["markers"][[1]])) {
        message("Found precomputed markers in obj@misc['markers']")
        markers <- object@misc["markers"]$markers
      } else {
        message("Running FindAllMarkers(), using wilcox test, min logfc diff 0.25")
        markers <- FindAllMarkers(
          object,
          do.print = TRUE,
          print.bar = TRUE,
          test.use = "wilcox",
          logfc.threshold = 0.25
        )
      }
      message("Writing top ", markers.n, ", cluster markers to ", fname)
      markers.order <- ave(x = rownames(x = markers), markers$cluster, FUN = markers.helper)
      top.markers <- markers[markers.order[!is.na(x = markers.order)], ]

      if (!is.null(ident.field)) {
        rownames(top.markers) <- NULL
        marker_df <- top.markers[, -which(colnames(top.markers) %in% c("cluster", "gene"))]
        marker_df <- cbind(gene = top.markers$`gene`, top.markers$`cluster`, marker_df)

        colnames(marker_df)[1] <- ident.field
      } else {
        rownames(top.markers) <- NULL
        marker_df <- top.markers[, -which(colnames(top.markers) %in% c("cluster", "gene"))]
        marker_df <- cbind(gene = top.markers$`gene`, cluster = top.markers$`cluster`, marker_df)
      }



      readr::write_tsv(x = marker_df, file = fname)
      ## this add by xtt
      path_list$marker_gene_meta_files <- list(fname)
    } else {
      message("No clusters found in Seurat object and no external marker file provided, so no marker genes can be computed")
      file <- NULL
    }
  }
  if (!is.null(file)) {
    markers.string <- sprintf(
      'markers = [{"file": "%s", "shortLabel": "Seurat Cluster Markers"}]',
      file
    )
  }

  ## return the merged dataframe and path list
  return(list(sc_id, cell_meta_column, cell_plots, total_df, path_list, dir))
}






#' ExportToSGS
#' @description this function used to add Seurat object to SGS browser
#' @param object Seurat object
#' @param species_id the id of the species
#' @param track_name the name of the single cell track
#' @param track_type the type of the single cell track, transcript or atac
#' @param select_group vector of cell group names to export.
#' @param reductions vector of reduction names to export, defaults to all reductions.
#' @param markers.file path to file with marker genes.
#' @param ident.field name of the idents used to caculate the marker gene when there is no marker file
#' @param matrix.slot matrix to use ("data" or "counts"), default is 'scale.data'.
#' @param markers.n no marker files were supplied, FindAllMarkers is run, markers.n is used to setting the number of gene to export.
#' @return This function exports Seurat object as a set of merged matrix feather file and the json format information to send post
#'
#' @importFrom Seurat Project Idents GetAssayData Embeddings FetchData DefaultAssay FindAllMarkers GetAssay
#' @importFrom methods .hasSlot
#' @importFrom stats ave
#' @importFrom utils install.packages packageVersion
#' @importFrom readr write_tsv
#' @importFrom feather write_feather
#' @importFrom httr POST add_headers
#' @importFrom uuid UUIDgenerate
#' @importFrom jsonlite write_json toJSON
#' @author xtt
#' @export
#'
#' @examples
#' \dontrun{
#'
#' data("pbmc_small")
#'
#' marker_file <- system.file("extdata", "markers.tsv", package = "SGSR")
#'
#' result_content <- ExportToSGS(
#'   object = pbmc_small,
#'   species_id = "ca4570d9c5464f10915e2def60b94780",
#'   track_name = "pbmc_small",
#'   track_type = "transcript",
#'   reductions = c("tsne", "umap"),
#'   markers.file = marker_file,
#'   select_group = c("groups", "RNA_snn_res.0.8","letter.idents","cluster"),
#'   matrix.slot = "data"
#' )
#' }
#'
ExportToSGS <- function(object,
                        species_id,
                        track_name,
                        track_type,
                        select_group,
                        reductions = NULL,
                        markers.file = NULL,
                        ident.field = NULL,
                        matrix.slot = "data",
                        markers.n = 100) {

  results <- gain_data(
    object,
    reductions,
    markers.file,
    ident.field,
    matrix.slot,
    markers.n)

  path_list <- results[[5]]
  path_list$all_meta_columns <- results[[2]] ## 默认输出全部meta列
  path_list$cell_plots <- results[[3]]
  path_list$species_id <- species_id
  path_list$sc_name <- track_name
  path_list$sc_type <- track_type
  path_list$select_meta_columns <- select_group
  path_list_j <- jsonlite::toJSON(path_list, auto_unbox = TRUE)


  # gain the dir to write the conf json file
  data_dir <- results[[6]]
  conf_path <- file.path(data_dir, "conf.json")
  jsonlite::write_json(path_list, path = conf_path)


  ########## send the post request
  header <- c(
    "Content-Type" = "application/json",
    "Accept" = "application/json, text/javascript, */*; q=0.01",
    "Accept-Encoding" = "gzip, deflate, br",
    "Connection" = "keep-alive"
  )


  post_url <- "http://sgs-api:6102/api/sc/add/seurat"
  # post_url <- "http://47.74.241.105:6102/api/sc/add/seurat"

  post_body <- gsub("/home/sgs/data", "", path_list_j)

  post_result <- httr::POST(url = post_url, body = post_body, encode = "json", add_headers(.headers = header))

  post_status <- httr::status_code(post_result)

  post_content <- httr::content(post_result)

  if (post_status == "200") {
    message("the single cell track load successful!")
  } else {
    ##delete the dir
    if(dir.exists(data_dir)){
      unlink(data_dir, recursive = TRUE)
    }
    stop("the single cell track add failed")
  }

  return(post_content)
}
