#' @noRd
filter_indirect_edges_type <- function(dir_out, type) {
    # format file path
    fpath_in <- file.path(dir_out, sprintf("MutualInfo_Typ%s_Para.txt", type))
    fpath_out <- file.path(dir_out, sprintf("IntracellularNetwork_Type%s.txt", type))

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_in)) {
        stop(sprintf("cannot find input file: %s", fpath_in))
    } else if (file.exists(fpath_out)) {
        return()
    }

    # load in file
    mat <- suppressMessages(vroom::vroom(fpath_in, progress = FALSE))
    # set first column to rownames
    mat <- tibble::column_to_rownames(mat, names(mat)[1])
    # as matrix
    mat <- as.matrix(mat)

    # zero out diagonal
    for (i in 1:ncol(mat)) {
        mat[i, i] <- 0
    }

    # compute ARACNE gene interaction network
    mat <- as.data.frame(parmigene::aracne.m(mat, 0.15))

    # write out
    rownames(mat) <- names(mat) <- names(mat)
    mat <- tibble::rownames_to_column(mat)
    vroom::vroom_write(mat, fpath_out, progress = FALSE)
}
    
#' Filter Out Indirect Edges
#'
#' Using the ARACNE.M function, filter edges so that only direct edges remain.
#'
#' @param dir_out Output directory
#' @return NULL
#' @export
filter_indirect_edges <- function(dir_out) {
    filter_indirect_edges_type(dir_out, "A")
    filter_indirect_edges_type(dir_out, "B")
    NULL
}
