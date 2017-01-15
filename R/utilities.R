#' Print information about the accession dictionary
#'
#' @export
print.accession.dict <- function(x, ...) {
    cat(sprintf("Dictionary of %d accession numbers in %d protein groups",
                nrow(x$dictionary), max(x$dictionary$pg)))
    if (!is.null(x$directory)) {
        cat(sprintf(' built from %d files in directory "%s"\n', length(x$files), x$directory))
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
#' @param dict an accession dictionary.
#' @param file either a character string naming a file or a \code{\link[base]{connection}}
#'      open for writing. \code{""} indicates output to the console.
#' @param ... further arguments passed to \code{\link[utils]{write.table}}.
#' @importFrom utils write.table
#'
#' @seealso \code{\link{accession.dictionary}} to create the dictionary, and
#'      \code{\link{translate}} to saved translated data files.
#'
#' @export
#'
#' @examples
#' # Build accession dictionary from all files in a directory
#' dict <- accession.dictionary.dir("dir",
#'                  col.mapping = c(gene.symbol = "Gene_Symbol"))
#'
#' \dontrun{
#' # Save dictionary to a file
#' save.dict(dict, file = "dictionary.txt")
#'
#' # Change the separator string
#' save.dict(dict, file = "dictionary.txt", sep = "\t")
#' }
save.dict <- function(dict, file = stop("`file` must be specified"), ...) {
    if (class(dict) != "accession.dict") {
        stop("`dict` must be an accession dictionary")
    }

    d <- dict$dictionary
    colnames(d) <- c("PG", na.omit(col.mapping[c("group.identifier", "accession.nr",
                                                 "protein.name", "gene.symbol")]))

    write.table(dict$dictionary, file = file, row.names = FALSE, ...)
}


#' Apply the dictionary to data files
#'
#' Apply the dictionary to the data files and write the translated files to disk.
#'
#' The dictionary is applied to the files used to create the dictionary. If
#' the output directory \code{out.dir} is specified, the translated files will be saved
#' in this directory. Otherwise the files will be written to the same directory as the
#' original file. The function will not overwrite existing files and will fail if the
#' files already exist. Parameters \code{out.suffix} and \code{out.prefix} can be used
#' to ensure unique new file names.
#'
#' In case the dictionary was built from \code{data.frame}s directly and no output
#' directory is specified, the \code{data.frame}s will be translated and returned as a
#' list. The function will also return the translated \code{data.frame}s as list if
#' the \code{out.dir = NULL}.
#'
#' @param dict accession dictionary to use.
#' @param out.dir the directory to save the translated files in (see details).
#'      If \code{NULL}, the translated data frames will be returned directly.
#' @param out.suffix,out.prefix suffix and prefix that will be added to the translated files.
#' @param col.mapping the column mapping for the input files. Defaults to the same as used
#'      to build the dictionary.
#' @param out.pg.col the name of the column to store the protein group.
#' @param ...
#'
#' @return Either a list of \code{data.frame}s or nothing (see details).
#' @seealso \code{\link{accession.dictionary}} to create the dictionary
#' @importFrom utils write.table
#' @export
#'
#' @examples
#' # Build accession dictionary from all files in a directory
#' dict <- accession.dictionary.dir("dir",
#'                  col.mapping = c(gene.symbol = "Gene_Symbol"))
#'
#' # Translate all files in the directory and return as a list of data.frames
#' trans <- translate(dict, out.dir = NULL)
#'
#' \dontrun{
#' # Translate all files in the directory and save to another directory
#' translate(dict, out.dir = "translated")
#'
#' # Translate all files in the directory, save to the same directory, but
#' # add a prefix to the files
#' translate(dict, out.prefix = "translated_")
#' }
#'
translate <- function(dict, out.dir, out.suffix = "", out.prefix = "", col.mapping,
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
    if (!is.null(dict$files)) {
        files <- dict$files
        dfs <- lapply(files, read.delim, header = TRUE, stringsAsFactors = FALSE)

        names(dfs) <- if (!is.null(names(files))) {
            names(files)
        } else {
            files
        }
    } else {
        lc <- dict$call
        lc[[1]] <- quote(list)
        if (!is.null(lc[["col.mapping"]])) {
            lc[["col.mapping"]] <- NULL
        }
        dfs <- eval(lc, parent.frame())
        files <- NULL
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

