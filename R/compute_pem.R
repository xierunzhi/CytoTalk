#' Compute Preferential Expression Measure
#'
#' Using all cell types in the input folder, compute the PEM; this is done per
#' gene, per cell type.
#'
#' @examples {
#' dir_in <- "~/scRNAseq-data"
#' dir_out <- "~/CytoTalk-output"
#'
#' compute_pem(dir_in, dir_out)
#' }
#'
#' @param dir_in Input directory, contains scRNAseq files
#' @param dir_out Output directory
#' @return None
#' @export
compute_pem <- function(dir_in, dir_out) {
    # format filepaths
    pattern <- "^scRNAseq_(.+)\\.csv$"
    fnames_scRNA <- list.files(dir_in, pattern = pattern)
    fpaths_scRNA <- file.path(dir_in, fnames_scRNA)
    fpath_out <- file.path(dir_out, "GeneCellTypeSpecific.txt")

    # skip if output already generated
    if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # for all files
    lst <- list()
    for (i in seq_len(length(fpaths_scRNA))) {
        # load in
        df <- suppressMessages(vroom::vroom(fpaths_scRNA[i], progress = FALSE))
        df <- tibble::column_to_rownames(df, names(df)[1])

        # exponentiate
        df <- expm1(df)

        # store in list
        lst[[i]] <- df
    }
    df <- NULL

    # rowmeans per scRNA file
    lst_rowmeans <- lapply(lst, rowMeans)
    # sums of rowmeans, per scRNA file
    lst_sums <- lapply(lst_rowmeans, sum)
    # sums of rowmeans, per gene (combined scRNA files)
    vec_gene_sums <- rowSums(do.call(cbind, lst_rowmeans))
    # overall sum of rowmeans
    total_sum <- sum(vec_gene_sums)

    # for all files
    pem_out <- list()
    for (i in seq_len(length(lst_rowmeans))) {
        pem_out[[i]] <- vector()

        # what proportion of this cell type's rowmean sum
        # accounts for the whole?
        cell_type_prop <- lst_sums[[i]] / total_sum

        # for all genes
        for (j in seq_len(length(vec_gene_sums))) {

            # scale gene sum to cell type proportion
            gene_prop <- vec_gene_sums[j] * cell_type_prop
            # scale cell type rowmean to gene proportion
            pem_out[[i]][j] <- log10(lst_rowmeans[[i]][j] / gene_prop)
        }
    }

    # join columns to dataframe
    pem_out <- as.data.frame(do.call(cbind, pem_out))

    # extract valid type names, copy over rownames
    names(pem_out) <- gsub(pattern, "\\1", fnames_scRNA)
    rownames(pem_out) <- rownames(lst[[1]])

    # write out
    pem_out <- tibble::rownames_to_column(pem_out)
    vroom::vroom_write(pem_out, fpath_out, progress = FALSE)

    # cleanup!
    rm(list = ls())
    gc()
    NULL
}
