#' @export
print.pgcaDict <- function(x, ...) {
    cat(sprintf("Dictionary mapping %d proteins to %d protein groups",
                nrow(x$dictionary), max(x$dictionary$pg)))

    concat <- "built from"

    if (!is.null(x$files) && length(x$files) > 0L) {
        cat(
            sprintf(' built from %d file%s',
                    length(x$files),
                    ifelse(length(x$files) > 1L, "s", "")),
            x$files,
            sep="\n\t"
        )
        concat <- "and"
    }

    if (!is.null(x$directories)) {
        cat(
            switch(
                (length(x$directories) > 1L) + 1L,
                'in the directory',
                'in the directories'
            ),
            x$directories,
            sep="\n\t"
        )
    }

    df.names <- setdiff(x$df.names, x$files)
    if (!is.null(df.names) & length(df.names) > 0L) {
        cat(
            sprintf(' %s %d data frame%s', concat, length(df.names),
                    ifelse(length(df.names) > 1L, "s", "")),
            df.names,
            sep="\n\t"
        )
    }
    cat("\n")
}

#' Write a PGCA dictionary to a text file
#'
#' Write the dictionary to a text file using the
#' \code{\link[utils]{write.table}} function. By default, a tab-separated file
#' is written, but this can be changed
#' by changing the arguments to \code{\link[utils]{write.table}}.
#'
#' @param dict a PGCA dictionary.
#' @param file either a character string naming a file or a
#'      \code{\link[base]{connection}} open for writing.
#'      \code{""} indicates output to the console.
#' @param ... further arguments passed to \code{\link[utils]{write.table}}.
#' @return This function returns \code{NULL} invisibly.
#' @importFrom utils write.table
#'
#' @seealso \code{\link{pgcaDict}} to create the dictionary, and
#'      \code{\link{applyDict}} to apply the dictionary for translating
#'      data files.
#'
#' @export
#'
#' @examples
#' # Build accession dictionary from all files in a directory
#' dict <- pgcaDict(
#'          system.file("extdata", package="pgca"),
#'          col.mapping=c(gene.symbol="Gene_Symbol")
#' )
#'
#' # Save dictionary to a temporary file
#' dictOutFile <- tempfile()
#' saveDict(dict, file=dictOutFile)
#'
#' # Change the separator string to a tab
#' dictOutFile <- tempfile()
#' saveDict(dict, file=dictOutFile, sep="\t")
saveDict <- function(dict, file=stop("`file` must be specified"), ...) {
    if (class(dict) != "pgcaDict") {
        stop("`dict` must be a PGCA dictionary")
    }

    cn <- dict$col.mapping[c("group.identifier", "accession.nr", "protein.name",
                             "gene.symbol")]
    d <- dict$dictionary
    colnames(d) <- c("PG", cn[!is.na(cn)])

    write.table(dict$dictionary, file=file, row.names=FALSE, ...)
    invisible(NULL)
}


#' Apply a PGCA dictionary to data files
#'
#' Apply the dictionary to the data files and write the translated files to
#' disk.
#'
#' The dictionary is applied to the data specified  argument.
#' If no input is provided, the dictionary is applied to the
#' files used to create the dictionary. The inputs can be
#' directory names, file names, or \code{data.frames}s.
#'
#' If the output directory \code{out.dir} is specified, the translated files
#' will be saved in this directory. Otherwise the files will be written to the
#' same directory as the input files. The function will not overwrite existing
#' files and will fail if the files already exist. Parameters \code{out.suffix}
#' and \code{out.prefix} can be used to ensure unique new file names.
#' In case the input is a list of \code{data.frame}s and no output
#' directory is specified, the \code{data.frame}s will be translated and
#' returned as a list. The function will also return the translated
#' \code{data.frame}s as list if the \code{out.dir=NULL}.
#'
#' @param ... input (see details).
#' @param dict the PGCA dictionary to use.
#' @param out.dir the directory to save the translated files in (see details).
#'      If \code{NULL}, the translated data frames will be returned directly.
#' @param out.suffix,out.prefix suffix and prefix that will be added to the
#'      translated files.
#' @param col.mapping the column mapping for the input files. Defaults to the
#'      same as used to build the dictionary.
#' @param out.pg.col the name of the column to store the protein group.
#'
#' @return Either a list of \code{data.frame}s or nothing (see details).
#' @seealso \code{\link{pgcaDict}} to create the dictionary
#' @importFrom utils write.table
#' @importFrom stats setNames
#' @export
#'
#' @examples
#' # Build PGCA dictionary from all files in a directory
#' dict <- pgcaDict(
#'          system.file("extdata", package="pgca"),
#'          col.mapping=c(gene.symbol="Gene_Symbol")
#' )
#'
#' # Translate all files in the directory and return as a list of data.frames
#' trans <- applyDict(system.file("extdata", package="pgca"), dict=dict,
#'                    out.dir=NULL)
#'
#' # Translate only some files in the directory and return as a list of
#' # data.frames
#' trans <- applyDict(
#'     system.file("extdata", "BET1947_v339.txt", package="pgca"),
#'     system.file("extdata", "BET2047_v339.txt", package="pgca"),
#'     dict=dict
#' )
#' str(trans)
#'
#' # Translate all files in the directory and save to another directory
#' out.dir <- tempdir()
#' applyDict(system.file("extdata", package="pgca"), dict=dict,
#'           out.dir=out.dir)
#'
applyDict <- function(..., dict, out.dir=NULL, out.suffix="",
                      out.prefix="", col.mapping, out.pg.col="PGC") {
    if (class(dict) != "pgcaDict") {
        stop("`dict` must be a PGCA dictionary")
    }

    # Validate the column mapping
    col.mapping <- if (missing(col.mapping)) {
        dict$col.mapping
    } else {
        .getColMapping(col.mapping)
    }

    ## Get the data
    args <- list(...)
    cl <- match.call(expand.dots=TRUE)
    if (length(args) == 0L) {
        rem.args <- names(formals(pgcaDict))
        rem.args <- rem.args[rem.args != "..."]
        cl <- dict$call[!(names(dict$call) %in% rem.args)]
        cl[[1L]] <- quote(list)
        args <- eval(cl, parent.frame())
    }

    ##
    ## Get names of arguments
    ##
    # remove function name and known arguments
    cl <- cl[!(names(cl) %in% names(formals()))][-1L]
    argnames <- unlist(lapply(cl, deparse))
    if (!is.null(names(argnames))) {
        argnames <- ifelse(names(argnames) == "", argnames, names(argnames))
    }

    ##
    ## Read in data
    ##
    args <- setNames(args, argnames)
    dfs <- .toDataFrame(args)

    ##
    ## Match accession nr's in files with dictionary
    ##
    dfs <- lapply(dfs, function(df) {
        df[[out.pg.col]] <- NA_integer_
        accs.matches <- match(df[[col.mapping["accession.nr"]]],
                              dict$dictionary$accs)
        df[[out.pg.col]] <- dict$dictionary$pg[accs.matches]
        return(df)
    })

    ##
    ## If out.dir is NULL, return the data directly
    ##
    if (is.null(out.dir)) {
        return(dfs)
    }

    ##
    ## Determine file names
    ##
    df.names <- names(dfs)
    df.file.names <- unlist(lapply(dfs, function (df) {
        fn <- attr(df, "file", exact=TRUE)
        if (is.null(fn)) {
            return(NA_character_)
        }
        return(sub("\\.[^\\.]+$", "", basename(fn)))
    }))
    df.names[!is.na(df.file.names)] <- df.file.names[!is.na(df.file.names)]

    out.files <-  file.path(
        out.dir,
        sprintf("%s%s%s.txt", out.prefix, basename(df.names), out.suffix)
    )

    if (any(file.exists(out.files))) {
        stop("output files already exist and will not be overwritten")
    }

    mapply(write.table, x=dfs, file=out.files, sep="\t",
           row.names=FALSE)

    invisible(NULL)
}

