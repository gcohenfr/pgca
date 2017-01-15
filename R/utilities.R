#' Print information about the accession dictionary
#'
#' @export
print.accession.dict <- function(x, ...) {
    cat(sprintf("Dictionary of %d accession numbers in %d protein groups",
                nrow(x$dictionary), max(x$dictionary$PG)))
    if (!is.null(x$dir)) {
        cat(sprintf(' built from %d files in directory "%s"\n', length(x$files), x$dir))
    } else if (!is.null(x$files)) {
        cat(" built from files:\n")
        cat(sprintf("\t%s\n", x$files))
    }
    cat("\n")
}

#' Write an accession dictionary to a text file
#'
#' Write the accession dictionary to a text file using the \code{\link[utils]{write.table}}
#' function. By default, a tab-separated file is written, but this can be changed
#' by changing the arguments to \code{\link[utils]{write.table}}.
#'
#' @param x an accession dictionary.
#' @param file either a character string naming a file or a \code{\link[base]{connection}}
#'      open for writing. \code{""} indicates output to the console.
#' @param ... further arguments passed to \code{\link[utils]{write.table}}.
#' @importFrom utils write.table
#' @export
save.dict <- function(x, file = stop("`file` must be specified"), ...) {
    if (class(x) != "accession.dict") {
        stop("`x` must be an accession dictionary")
    }

    write.table(x$dictionary, file = file, row.names = FALSE, ...)
}

