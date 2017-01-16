#' Print information about the PGCA dictionary
#'
#' @param x the PGCA dictionary.
#' @param ... currently not used.
#' @export
print.pgca.dict <- function(x, ...) {
    cat(sprintf("Dictionary mapping %d proteins to %d protein groups",
                nrow(x$dictionary), max(x$dictionary$pg)))
    if (!is.null(x$directory)) {
        cat(sprintf(' built from %d files in directory "%s"\n', length(x$files), x$directory))
    } else if (!is.null(x$files)) {
        cat(" built from files:\n")
        cat(sprintf("\t%s\n", x$files))
    }
    cat("\n")
}

#' Write PGCA dictionary to a text file
#'
#' Write the dictionary to a text file using the \code{\link[utils]{write.table}}
#' function. By default, a tab-separated file is written, but this can be changed
#' by changing the arguments to \code{\link[utils]{write.table}}.
#'
#' @param dict a PGCA dictionary.
#' @param file either a character string naming a file or a \code{\link[base]{connection}}
#'      open for writing. \code{""} indicates output to the console.
#' @param ... further arguments passed to \code{\link[utils]{write.table}}.
#' @importFrom utils write.table
#'
#' @seealso \code{\link{pgca}} to create the dictionary, and
#'      \code{\link{translate}} to saved translated data files.
#'
#' @export
#'
#' @examples
#' # Build accession dictionary from all files in a directory
#' dict <- pgca(
#'          system.file("extdata", package = "pgca"),
#'          col.mapping = c(gene.symbol = "Gene_Symbol")
#' )
#'
#' \dontrun{
#' # Save dictionary to a file
#' save.dict(dict, file = "dictionary.txt")
#'
#' # Change the separator string to a tab
#' save.dict(dict, file = "dictionary.txt", sep = "\t")
#' }
save.dict <- function(dict, file = stop("`file` must be specified"), ...) {
    if (class(dict) != "pgca.dict") {
        stop("`dict` must be a PGCA dictionary")
    }

    cn <- dict$col.mapping[c("group.identifier", "accession.nr", "protein.name", "gene.symbol")]
    d <- dict$dictionary
    colnames(d) <- c("PG", cn[!is.na(cn)])

    write.table(dict$dictionary, file = file, row.names = FALSE, ...)
}


#' Apply a PGCA dictionary to data files
#'
#' Apply the dictionary to the data files and write the translated files to disk.
#'
#' The dictionary is applied to the data specified in the \code{input} argument.
#' If the \code{input} argument is missing, the dictionary is applied to the files used
#' to create the dictionary. The \code{input} argument can either be a directory,
#' a character vector of file names, or a list of \code{data.frames}s.
#'
#' If the output directory \code{out.dir} is specified, the translated files will be saved
#' in this directory. Otherwise the files will be written to the same directory as the
#' input files. The function will not overwrite existing files and will fail if the
#' files already exist. Parameters \code{out.suffix} and \code{out.prefix} can be used
#' to ensure unique new file names.
#' In case the input is a list of \code{data.frame}s and no output
#' directory is specified, the \code{data.frame}s will be translated and returned as a
#' list. The function will also return the translated \code{data.frame}s as list if
#' the \code{out.dir = NULL}.
#'
#' @param dict the PGCA dictionary to use.
#' @param input input directory, files, or data (see details).
#' @param out.dir the directory to save the translated files in (see details).
#'      If \code{NULL}, the translated data frames will be returned directly.
#' @param out.suffix,out.prefix suffix and prefix that will be added to the translated files.
#' @param col.mapping the column mapping for the input files. Defaults to the same as used
#'      to build the dictionary.
#' @param out.pg.col the name of the column to store the protein group.
#'
#' @return Either a list of \code{data.frame}s or nothing (see details).
#' @seealso \code{\link{pgca}} to create the dictionary
#' @importFrom utils write.table
#' @export
#'
#' @examples
#' # Build PGCA dictionary from all files in a directory
#' dict <- pgca(
#'          system.file("extdata", package = "pgca"),
#'          col.mapping = c(gene.symbol = "Gene_Symbol")
#' )
#'
#' # Translate all files in the directory and return as a list of data.frames
#' trans <- translate(dict, input = system.file("extdata", package = "pgca"),
#'                    out.dir = NULL)
#'
#' # Translate only some files in the directory and return as a list of data.frames
#' trans <- translate(dict, input = c(
#'                        system.file("extdata", "accs_no_1947.txt", package = "pgca"),
#'                        ssystem.file("extdata", "accs_no_2007.txt", package = "pgca")
#'                    ), out.dir = NULL)
#' str(trans)
#'
#' \dontrun{
#' # Translate all files in the directory and save to another directory
#' translate(dict, input = system.file("extdata", package = "pgca"),
#'           out.dir = "translated")
#'
#' # Translate all files in the directory, save to the same directory, but
#' # add a prefix to the files
#' translate(dict, input = system.file("extdata", package = "pgca"),
#'           out.prefix = "translated_")
#' }
#'
translate <- function(dict, input, out.dir, out.suffix = "", out.prefix = "", col.mapping,
                      out.pg.col = "PGC") {
    # Validate the column mapping
    col.mapping <- if (missing(col.mapping)) {
        dict$col.mapping
    } else {
        .get.col.mapping(col.mapping)
    }

    ##
    ## Read in data
    ##
    if (missing(input)) {
        if (!is.null(dict$files)) {
            input <- dict$files
        } else {
            lc <- dict$call
            lc[[1]] <- quote(list)
            if (!is.null(lc[["col.mapping"]])) {
                lc[["col.mapping"]] <- NULL
            }
            input <- eval(lc, parent.frame())
        }
    }

    if (is.character(input)) {
        input <- normalizePath(input, mustWork = TRUE)
        input.finfo <- file.info(input, extra_cols = FALSE)

        dirs <- input[input.finfo$isdir]
        files <- input[!input.finfo$isdir]

        files <- c(files, unlist(lapply(dirs, function(x) {
            nodes <- dir(x, full.names = TRUE)
            dirs.in.dir <- list.dirs(x, recursive = FALSE)
            setdiff(nodes, dirs.in.dir)
        })))

        if (length(files) == 0L) {
            stop("`input` directory is empty")
        }

        dfs <- lapply(files, read.delim, header = TRUE, stringsAsFactors = FALSE)

        names(dfs) <- if (!is.null(names(files))) {
            names(files)
        } else {
            files
        }
    } else if (is.list(input)) {
        if (!all(sapply(input, class) == "data.frame")) {
            stop('one or more list items of `input` are not data.frames.')
        }
        dfs <- input
        files <- NULL
    } else {
        stop('`input` not recognized.')
    }

    ##
    ## Match accession nr's in files with dictionary
    ##
    dfs <- lapply(dfs, function(df) {
        df[[out.pg.col]] <- NA_integer_
        accs.matches <- match(df[[col.mapping["accession.nr"]]], dict$dictionary$accs)
        df[[out.pg.col]] <- dict$dictionary$pg[accs.matches]
        return(df)
    })

    ##
    ## If out.dir is NULL, return the data directly
    ##
    if (!missing(out.dir) && is.null(out.dir)) {
        return(dfs)
    }

    ##
    ## Determine file names
    ##
    if (missing(out.dir) && !is.null(files)) {
        out.dir <- dirname(files)
    }
    out.files <- if (!is.null(files)) {
        fcenter <- sub("\\.[a-zA-Z]+$", "", basename(files))
        fext <- sub(".+\\.([a-zA-Z]+)$", "\\1", basename(files))
        file.path(out.dir, sprintf("%s%s%s.%s", out.prefix, fcenter, out.suffix, fext))
    }

    if (any(file.exists(out.files))) {
        stop("output files already exist and will not be overwritten")
    }

    mapply(write.table, x = dfs, file = out.files, sep = "\t", row.names = FALSE)

    invisible(NULL)
}

