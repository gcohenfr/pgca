#' Link Protein Groups Created from MS/MS Data
#'
#' Build a dictionary for protein groups from MS/MS data.
#' Details of the algorithm can be found in Takhar et al. (Under revision).
#' "PGCA: An Algorithm to Link Protein Groups Created from MS/MS Data.".
#'
#' If the \code{group.identifier} column is logical (i.e., \code{TRUE} or
#' \code{FALSE}) or \code{master.gene.identifier} is given,
#' the \code{TRUE} accessions are assumed to be a "master gene" and the
#' data set is assumed to be in the correct order. This means all
#' \code{FALSE} values following the master gene are assumed to belong to
#' the same group.
#'
#' The \code{col.mapping} maps the column names in the data files to a specific
#' function. It nees to be a named character vector, whereas the name of each
#' item is the "function" of the given column name. The algorithm knows about
#' the following columns:
#' \describe{
#'      \item{\code{"group.identifier"}}{Column containing the group
#'                                       identifier.}
#'      \item{\code{"accession.nr"}}{Column containing the accession nr.}
#'      \item{\code{"protein.name"}}{Column containing the protein name.}
#'      \item{\code{"gene.symbol"}}{Column containing the gene symbol (if any,
#'                                  can be missing)}
#' }
#' The default column mapping is \code{c(group.identifier="N", accession.nr =
#' "Accession", protein.name="Protein_Name")}. The supplied column mapping can
#' ignore those columns that are already correct in the default map.
#' For instance, if the accession nr. is stored in column \emph{AccessionNr}
#' instead of \emph{Accession}, but the remaining columns are the same as in the
#' default mapping, specifying \code{col.mapping=c(accession.nr="AccessionNr")}
#' is sufficient.
#'
#' @param ... arbitrary number of directory names, file names, or
#'       \code{data.frame}s used as input.
#' @param col.mapping column mapping (see Details).
#' @param master.gene.identifier if given, genes with this value in the
#'      column \code{group.identifier} are considered master genes.
#'
#' @return An object of type \code{pgcaDict}.
#' @importFrom stats setNames
#' @export
#'
#' @references Takhar M, Sasaki M, Hollander Z, McManus B, McMaster W, Ng R and
#'    Cohen Freue G (Under revision). "PGCA: An Algorithm to Link Protein Groups
#'    Created from MS/MS Data." PLOS ONE.
#' @seealso \code{\link{applyDict}} to apply the dictionary to the data files
#'      and \code{\link{saveDict}} to save the dictionary itself.
#'
#' @examples
#' # Build the dictionary from all files in a directory
#' # and specifying the column "Gene_Symbol" holds the `gene.symbol`.
#' dict.dir <- pgcaDict(
#'          system.file("extdata", package="pgca"),
#'          col.mapping=c(gene.symbol="Gene_Symbol")
#' )
#'
#' # Build the dictionary from a list of files
#' dict.files <- pgcaDict(
#'      system.file("extdata", "BET1947_v339.txt", package="pgca"),
#'      system.file("extdata", "BET2007_v339.txt", package="pgca"),
#'      system.file("extdata", "BET2047_v339.txt", package="pgca"),
#'      col.mapping=c(gene.symbol="Gene_Symbol")
#' )
#'
#' # Build the dictionary from already read-in data.frames
#' dict.data <- pgcaDict(BET1947_v339, BET2047_v339,
#'                       col.mapping=c(gene.symbol="Gene_Symbol"))
#'
pgcaDict <- function(..., col.mapping, master.gene.identifier) {
    ##
    ## Get names of arguments
    ##
    cl <- match.call(expand.dots=TRUE)
    # remove function name and known arguments
    cl <- cl[!(names(cl) %in% names(formals()))][-1L]
    argnames <- unlist(lapply(cl, deparse))
    if (!is.null(names(argnames))) {
        argnames <- ifelse(names(argnames) == "", argnames, names(argnames))
    }

    args <- setNames(list(...), argnames)

    # Parse input to data.frames
    dfs <- .toDataFrame(args)

    # Remove
    directories <- unique(unlist(lapply(dfs, attr, "directory", exact=FALSE),
                                 use.names = FALSE))
    # Collect file names
    files <- unlist(lapply(dfs, function (df) {
        f <- attr(df, "file", exact=TRUE)
        if (is.null(f)) {
            f <- NA_character_
        }
        return(f)
    }), use.names = FALSE)

    # Remove duplicated files
    dfs <- dfs[is.na(files) | !duplicated(files)]
    files <- unique(files[!is.na(files)])

    # Validate column mapping
    col.mapping <- .getColMapping(col.mapping)

    ret.obj <- .buildDict(
        dfs,
        col.mapping=col.mapping,
        master.gene.identifier=master.gene.identifier
    )
    ret.obj$directories <- directories
    ret.obj$files <- files
    ret.obj$df.names <- names(dfs)
    ret.obj$call <- match.call(expand.dots=TRUE)
    return(ret.obj)
}

## Work-horse function to create the dictionary
.buildDict <- function(dfs, col.mapping, master.gene.identifier) {
    ##
    ## Validate all data frames
    ##
    dfs <- unlist(lapply(dfs, function(df) {
        if (is.list(df) && !is.data.frame(df)) {
            return(df)
        }
        return(list(df))
    }), recursive=FALSE, use.names=TRUE)

    dfs.names <- if (!is.null(names(dfs))) {
        names(dfs)
    } else {
        as.character(seq_along(dfs))
    }

    has.mgi <- TRUE
    mig.arg <- NULL
    if (missing(master.gene.identifier)) {
        has.mgi <- FALSE
    } else{
        mig.arg <- master.gene.identifier
    }

    ## Validate data frame structure
    dfs <- mapply(
        .validateInputDF,
        df=dfs,
        name=dfs.names,
        MoreArgs=list(
            col.mapping=col.mapping,
            master.gene.identifier=mig.arg
        ),
        SIMPLIFY=FALSE
    )

    ##
    ## Work with a smaller data set since we only need to have
    ## the group id and the accessions to build the dictionary
    ##
    all.accessions <- lapply(dfs, function(df) {
        accs <- df[[col.mapping["accession.nr"]]]
        uq <- !duplicated(accs)
        rbind(
            accs=accs[uq],
            prot=df[[col.mapping["protein.name"]]][uq],
            gene=df[[col.mapping["gene.symbol"]]][uq],
            img=df[[col.mapping["img"]]][uq]
        )
    })
    all.accessions <- do.call(cbind, all.accessions)
    all.accessions <- all.accessions[ , sort.list(all.accessions["accs", ],
                                                  method="radix")]

    has.mgi <- (col.mapping["img"] %in% rownames(all.accessions))

    if (has.mgi) {
        # `img` is set to TRUE whenever at least one occurence was "TRUE"
        any.img <- tapply(all.accessions["img", ] == "TRUE",
                          all.accessions["accs", ], any)
    }

    all.accessions <- all.accessions[ , !duplicated(all.accessions["accs", ])]
    if (has.mgi) {
        all.accessions["img", ] <- any.img[all.accessions["accs", ]]
    }

    stripped <- lapply(dfs, function(df) {
        accs <- factor(df[[col.mapping["accession.nr"]]],
                       levels=all.accessions["accs", ])
        matrix(
            c(as.integer(df[[col.mapping["group.identifier"]]]),
              as.integer(accs),
              rep.int(NA_integer_, nrow(df))),
            ncol=3L,
            byrow=FALSE,
            dimnames=list(
                NULL,
                c("gid", "accs", "pg")
            )
        )
    })

    ##
    ## Merge groups in individual files according to overlapping Acc#s
    ##
    merged <- lapply(stripped, function(x) {
        accs.cont <- tabulate(x[ , "accs"], nbins=max(x[ , "accs"]))
        accs.ids <- which(accs.cont > 1L)
        if (length(accs.ids) > 0L) {
            for (k in accs.ids) {
                rep.n <- x[x[, "accs"] == k, "gid"]
                x[x[ , "gid"] %in% rep.n, "gid"] <- min(rep.n)
            }
        }

        x <- x[sort.list(x[ , "gid"], method="radix"), ]
        x <- cbind(x, rank=x[ , "gid"])
        x[duplicated(x[, "rank"]), "rank"] <- 0L

        return(x)
    })

    state.env <- new.env(size=5L)
    state.env$dict <- merged[[1L]][!duplicated(merged[[1L]][, "accs"]), ]
    state.env$pg.counter <- max(state.env$dict[ , "gid"])

    state.env$dict[, "pg"] <- state.env$dict[, "gid"]

    update.dict <- function(x, state.env) {
        for (i in seq_len(nrow(x))) {
            #check if the acc# exists already
            accs.matches <- which(state.env$dict[ , "accs"] == x[i, "accs"])
            if (length(accs.matches) == 0L) {
                # check if rank == 0, if so: it is an isoform, otherwise it is
                # a new PG
                if (x[i, "rank"] == 0L) {
                    # (i - 1) > 0 because the first row is always a new protein,
                    # thus rank != 0
                    x[i, "pg"] <- x[i - 1L, "pg"]
                } else {
                    state.env$pg.counter <- state.env$pg.counter + 1L
                    x[i, "pg"] <- state.env$pg.counter
                }
            } else {
                x[i, "pg"] <- state.env$dict[accs.matches[1L], "pg"]
            }
            #end replacing NA with PG translation
        }

        # to avoid new isoforms to be called a new PG:
        pos.rank <- which(x[ , "rank"] > 0L)
        # add total number of rows + 1 (CORRECT IN PREVIOUS VERSIONS!!!!) at
        # the end
        pos.rank <- c(pos.rank, nrow(x) + 1L)

        for(m in seq_len(length(pos.rank) - 1L)) {
            # If a group is merged, then we use the minimum PG# in that group
            # and we need to change the PG# in ALL that set and in the previous
            # version of the dictionary.
            pos.seq <- seq.int(pos.rank[m], pos.rank[m + 1L] - 1L)
            new.pg <- min(x[pos.seq, "pg"])
            pg.mismatch <- which(x[pos.seq, "pg"] != new.pg)
            if (length(pg.mismatch) > 0L) {
                pg.mismatch <- x[pos.seq[pg.mismatch], "pg"]

                for (pgm in pg.mismatch) {
                    x[x[ , "pg"] == pgm, "pg"] <- new.pg

                    # change current dictionary as well: some Acc# may not be
                    # present in this set but in the dictionary, thus, the old
                    # PG# will not be changed there.
                    state.env$dict[state.env$dict[ , "pg"] == pgm, "pg"] <-
                        new.pg
                }
            }
        }

        # combine new set with old dictionary
        state.env$dict <- rbind(state.env$dict, x)

        # re-order to and find duplicated acc#
        state.env$dict <- with(state.env, {
            dict[order(dict[, "accs"], -dict[, "rank"]), ]
        })
        dupl.accs <- duplicated(state.env$dict[, "accs"])
        state.env$dict <- state.env$dict[!dupl.accs, ]

        #re-order
        state.env$dict <- with(state.env, {
            dict[order(dict[, "pg"], -dict[, "rank"]), ]
        })
        return(NULL)
    }

    lapply(merged[-1L], update.dict, state.env=state.env)

    dict <- state.env$dict
    remove(state.env)

    dict.df <- data.frame(
        dict[ , c("pg", "gid")],
        t(all.accessions[ , dict[ , "accs"]]),
        stringsAsFactors=FALSE
    )

    if (has.mgi) {
        dict.df$gid <- (dict.df$img == TRUE)
        dict.df$img <- NULL
        colnames(dict.df)[2L] <- col.mapping[["group.identifier"]]
    }

    return(structure(list(
        dictionary=dict.df,
        col.mapping=col.mapping,
        call=match.call(expand.dots=TRUE)
    ), class="pgcaDict"))
}

## Helper function to validate a single input data frame
## This functions checks that
##   - each item is a data.frame
##   - each data.frame has the correct columns
##   - each column has the correct type
.validateInputDF <- function(df, name, col.mapping,
                             master.gene.identifier = NULL) {
    col.mapping.na <- col.mapping[!is.na(col.mapping) &
                                      names(col.mapping) != "img"]

    has.mgi <- !is.null(master.gene.identifier)

    if (!is.data.frame(df)) {
        stop(sprintf("Data set '%s' must be a data.frame", name))
    }

    present <- col.mapping.na %in% colnames(df)
    if (!all(present)) {
        stop(sprintf("Data set '%s' does not contain the column%s: '%s'",
                     name,
                     ifelse(sum(!present) > 1L, "s", ""),
                     paste(col.mapping.na[!present], collapse="', '")))
    }

    if (!is.character(df[[col.mapping["accession.nr"]]]) &&
        !is.factor(df[[col.mapping["accession.nr"]]])) {
        warning(sprintf("Accession Nr. column '%s' in data set '%s' is not a character or a factor",
                        col.mapping["accession.nr"], name))
    }
    if (!is.character(df[[col.mapping["accession.nr"]]])) {
        df[[col.mapping["accession.nr"]]] <-
            as.character(df[[col.mapping["accession.nr"]]])
    }

    # If the `master.gene.identifier` is given, the group.identifier
    # is converted to logical by comparing to the given value
    if (has.mgi) {
        df[[col.mapping["group.identifier"]]] <-
            df[[col.mapping["group.identifier"]]] == master.gene.identifier
    }

    # If the group identifier is a logical, we assume a "is master-gene"
    # type of column
    if (is.logical(df[[col.mapping["group.identifier"]]])) {
        if (!isTRUE(df[[col.mapping["group.identifier"]]][1L])) {
            stop(sprintf("The first 'is master-gene' value in column '%s' in data set '%s' must be `TRUE`",
                         col.mapping["group.identifier"], name))
        }

        df[[col.mapping["img"]]] <- df[[col.mapping["group.identifier"]]]

        df[[col.mapping["group.identifier"]]] <-
            cumsum(df[[col.mapping["group.identifier"]]])
    } else if (is.character(df[[col.mapping["group.identifier"]]])) {
        df[[col.mapping["group.identifier"]]] <-
            as.factor(df[[col.mapping["group.identifier"]]])
    }

    # If the gene.symbol or protein.name is a factor, convert it to
    # character (otherwise the integer factor level would be converted to
    # character)
    if (is.factor(df[[col.mapping["protein.name"]]])) {
        df[[col.mapping["protein.name"]]] <-
            as.character(df[[col.mapping["protein.name"]]])
    }
    if (is.factor(df[[col.mapping["gene.symbol"]]])) {
        df[[col.mapping["gene.symbol"]]] <-
            as.character(df[[col.mapping["gene.symbol"]]])
    }

    return(df)
}

## Helper function to convert a list of arguments
## (can be data.frames, directory names, or file names)
## into a flat list of data.frames
#' @importFrom utils read.delim
.toDataFrame <- function(input.list) {
    fun.envir <- environment()
    dfs <- lapply(input.list, function (arg) {
        if (is.data.frame(arg)) {
            # data.frames are returned directly
            return(arg)
        } else if (is.list(arg)) {
            # lists are recursively parsed
            return(do.call(.toDataFrame, arg, envir = fun.envir))
        } else if (is.character(arg)) {
            # character strings can be either directory
            # names or file names

            # Check if `arg` points to an existing path
            arg <- path.expand(arg)
            if (!file.exists(arg)) {
                stop(sprintf("File/directory %s does not exist", arg))
            }

            finfo <- file.info(arg)
            files <- dir(arg, full.names=TRUE)

            # Check if `arg` points to a (populated) directory
            if (isTRUE(finfo$isdir)) {
                if (length(files) == 0L) {
                    stop(sprintf("Directory %s is empty", dir))
                }
                # Check if any of the items in the directory are directories itself and
                # remove them
                dirs.in.dir <- list.dirs(arg, recursive=FALSE)
                files <- setdiff(files, dirs.in.dir)
            } else {
                files <- arg
            }

            dfs <- lapply(files, read.delim, header=TRUE, stringsAsFactors=FALSE,
                          check.names=FALSE)

            names(dfs) <- if (!is.null(names(files))) {
                names(files)
            } else {
                 files
            }

            if (isTRUE(finfo$isdir)) {
                attr(dfs[[1L]], "directory") <- arg
            }

            for (i in seq_along(dfs)) {
                attr(dfs[[i]], "file") <- files[i]
            }

            return(dfs)
        }
    })

    dfs <- .flattenDFList(dfs)

    return(dfs)
}

## Helper function to flatten a (possibly nested)
## list of data.frames into a flat list of data.frames
#' @importFrom stats setNames
.flattenDFList <- function(obj) {
    out <- NULL
    obj.names <- names(obj)

    for (i in seq_along(obj)) {
        if (is.data.frame(obj[[i]])) {
            out <- c(out, setNames(
                list(obj[[i]]),
                obj.names[i]
            ))
        } else if (is.list(obj[[i]])) {
            out <- c(out, .flattenDFList(obj[[i]]))
        }
    }

    return(out)
}

## Helper function to merge the user-supplied column mapping
## with the default column mapping and validate the mapping.
.getColMapping <- function(col.mapping) {
    default.col.mapping <- c(
        group.identifier="N",
        accession.nr="Accession",
        protein.name="Protein_Name",
        gene.symbol=NA_character_
    )

    if (!missing(col.mapping)) {
        override <- match(names(col.mapping), names(default.col.mapping))
        default.col.mapping[override] <- col.mapping
    }

    if (anyNA(default.col.mapping[c("group.identifier", "accession.nr",
                                    "protein.name")])) {
        stop("The columns `group.identifier`, `accession`, and ",
             "`protein.name` must be present")
    }

    # Find unique column name for the "is master gene" identifier"
    img.col.name <- "img"
    while (img.col.name %in% default.col.mapping) {
        img.col.name <- sprintf("img_%05d", sample.int(99999L, 1L))
    }
    default.col.mapping <- c(default.col.mapping, "img"=img.col.name)

    return (default.col.mapping)
}

