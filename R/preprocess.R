#' @noRd
filter_file <- function(fpath_in, fpath_out, cutoff, proteins) {
    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_in)) {
        stop(sprintf("cannot find input file: %s", fpath_in))
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return(0)
    }

    # load in file
    df <- suppressMessages(vroom::vroom(fpath_in, progress = FALSE))
    # set first column to rownames
    df <- tibble::column_to_rownames(df, names(df)[1])

    # calculate threshold
    thresh <- floor(ncol(df) * cutoff)
    # filter for number of zeros in a row
    index1 <- rowSums(df != 0) > thresh
    # filter for only protein-coding genes,
    # casefold upper
    index2 <- !is.na(match(toupper(rownames(df)), proteins))

    # apply filters
    out <- df[index1 & index2, ]

    # reintegrate rownames
    out <- tibble::rownames_to_column(out)
    # write out file
    vroom::vroom_write(out, fpath_out, progress = FALSE)

    NULL
}

#' Filter Lowly Expressed Genes
#'
#' Given cell types, cutoff values, and the names of protein-coding genes,
#' filter scRNAseq data so that only protein-coding genes with expression levels
#' over the cutoff remain.
#'
#' @examples \dontrun{
#' proteins <- CytoTalk::pcg_mouse
#' type_a <- "BCells"
#' type_b <- "TCells"
#' cutoff_a <- 0.1
#' cutoff_b <- 0.1
#' dir_in <- "~/scRNAseq-data"
#' dir_out <- "~/CytoTalk-output"
#'
#' preprocess(
#'     proteins, type_a, type_b, cutoff_a, cutoff_b, dir_in, dir_out
#' )
#' }
#'
#' @param proteins Character vector of protein-coding gene names
#' @param type_a Cell type A
#' @param type_b Cell type B
#' @param cutoff_a Minimum proportion of non-zero columns, per row, cell type A
#' @param cutoff_b Minimum proportion of non-zero columns, per row, cell type B
#' @param dir_in Input directory, contains scRNAseq files
#' @param dir_out Output directory
#' @return None
#' @export
preprocess <- function(
    proteins, type_a, type_b, cutoff_a, cutoff_b, dir_in, dir_out) {
    
    # casefold upper
    proteins <- toupper(proteins)

    # filter cell A genes
    filter_file(
        file.path(dir_in, sprintf("scRNAseq_%s.csv", type_a)),
        file.path(dir_out, "TypAExp_rmRepGf_dm.txt"),
        cutoff_a, proteins
    )

    # filter cell B genes
    filter_file(
        file.path(dir_in, sprintf("scRNAseq_%s.csv", type_b)),
        file.path(dir_out, "TypBExp_rmRepGf_dm.txt"),
        cutoff_b, proteins
    )

    NULL
}
