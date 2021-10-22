#' @noRd
now <- function() {
    format(Sys.time(), "%H:%M:%S")
}

#' @noRd
tick <- function(step, msg) {
    cat("[", step, " / 8] (", now(), ") ", msg, "\n", sep = "")
}

#' @noRd
minmax <- function(x, ...) {
    (x - min(x, ...)) / (max(x, ...) - min(x, ...))
}

#' Extract Valid Cell Types
#'
#' Searches input directory for a particular regular expression pattern,
#' filtering out invalid filenames, and finally returning a list of valid cell
#' types.
#'
#' @examples {
#' dir_in <- "~/scRNAseq-data"
#' check_valid_names(dir_in)
#' }
#'
#' @param dir_in Input directory, contains scRNAseq files
#' @return A character vector, includes only valid scRNAseq cell types
#' @export
check_valid_names <- function(dir_in) {
    # must have valid data directory
    if (!dir.exists(dir_in)) {
        stop("no input directory found")
    }

    # get filenames
    fnames <- dir(dir_in)

    # check for validity
    pattern <- "^scRNAseq_(.+)\\.csv$"
    index <- grepl(pattern, fnames)
    fnames_valid <- fnames[index]

    # extract valid type names
    type_names <- gsub(pattern, "\\1", fnames_valid)
    type_names
}
