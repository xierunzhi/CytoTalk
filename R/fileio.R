#' @noRd
new_named_list <- function(mat, cell_types) {
    list(mat = mat, cell_types = cell_types)
}

#' @noRd
vroom_silent <- function(...) {
    suppressMessages(vroom::vroom(..., progress = FALSE))
}

#' @noRd
vroom_write_silent <- function(x, file, rownames=FALSE) {
    x <- as.data.frame(x)
    if (rownames) { x <- tibble::rownames_to_column(x) }
    suppressMessages(vroom::vroom_write(x, file, progress = FALSE))
}

#' @noRd
vroom_with_rownames <- function(..., row_names=1) {
    dat <- vroom_silent(...)
    tibble::column_to_rownames(dat, names(dat)[row_names])
}

#' @noRd
vroom_sparse_with_rownames <- function(..., row_names=1) {
    dat <- vroom_with_rownames(..., row_names = row_names)
    Matrix::Matrix(Matrix::as.matrix(dat), sparse = TRUE)
}

#' @rdname doc_fileio
#' @export
read_matrix_folder <- function(
    dpath, pattern=".*scRNAseq_(.+?)\\..+", auto_transform=TRUE) {

    # initial parameters
    mat <- NULL
    msg <- "not all rownames identical between input files"

    # determine filepaths
    fpaths <- dir_full(dpath, pattern = pattern)

    # read in all files
    for (fpath in fpaths) {
        cell_type <- gsub(pattern, "\\1", fpath)
        if (is.null(mat)) {
            # read in
            mat <- vroom_sparse_with_rownames(fpath)
            # start cell types vector
            cell_types <- rep(cell_type, ncol(mat))
        } else {
            # read in
            new <- vroom_sparse_with_rownames(fpath)
            # check for identical rownames
            errorifnot(identical(rownames(mat), rownames(new)), msg)
            # accumulate cell type names
            cell_types <- c(cell_types, rep(cell_type, ncol(new)))
            # combine all data
            mat <- cbind(mat, new)
        }
    }

    # check for count data
    mat <- check_count_data(mat, auto_transform)

    # return sparse matrix
    new_named_list(mat, cell_types)
}

#' @rdname doc_fileio
#' @export
read_matrix_with_meta <- function(fpath_mat, fpath_meta, auto_transform=TRUE) {
    # read in
    mat <- vroom_sparse_with_rownames(fpath_mat)
    meta <- vroom_with_rownames(fpath_meta)

    # match meta to matrix
    index <- match(rownames(meta), colnames(mat))
    # ensure index matches exactly
    errorifnot(!any(is.na(index)), "meta file does not match matrix colnames")
    # otherwise, reorder
    cell_types <- meta[index, 1]

    # check for count data
    mat <- check_count_data(mat, auto_transform)

    # return sparse matrix
    new_named_list(mat, cell_types)
}
