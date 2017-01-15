#' Create a dictionary for the accession groups
#'
#' If the \code{group.identifier} column is logical (i.e., \code{TRUE} or \code{FALSE}),
#' the \code{TRUE} accessions are assumed to be a "master gene" and the
#' data set is assumed to be in the correct order (i.e., all \code{FALSE} values following
#' the master gene are assumed to belong to the same group).
#'
#' The \code{col.mapping} maps the column names in the data files to a specific function.
#' It nees to be a named character vector, whereas the name of each item is the "function"
#' of the given column name. The accession dictionary knows about the following columns:
#' \describe{
#'      \item{\code{"group.identifier"}}{Column containing the group identifier.}
#'      \item{\code{"accession.nr"}}{Column containing the accession nr.}
#'      \item{\code{"protein.name"}}{Column containing the protein name.}
#'      \item{\code{"gene.symbol"}}{Column containing the gene symbol (if any, can be missing)}
#' }
#' The default column mapping is \code{c(group.identifier = "N", accession.nr = "Accession",
#' protein.name = "Protein_Name")}. The supplied column mapping can miss those columns that are
#' already correct in the default map. For instance, if the accession nr. is stored in column
#' \emph{AccessionNr} instead of \emph{Accession}, but the remaining columns are the same as in the
#' default mapping, specifying \code{col.mapping = c(accession.nr = "AccessionNr")} is
#' sufficient.
#'
#' @param ... arbitrary number of \code{data.frame}s or filenames.
#' @param col.mapping column mapping (see Details).
#' @export
#' @seealso \code{\link{translate}} to apply the dictionary to the data files and
#'      \code{\link{save.dict}} to save the dictionary itself.
#' @examples
#' # Build accession dictionary from all files in a directory
#' dict.dir <- accession.dictionary.dir("dir",
#'                  col.mapping = c(gene.symbol = "Gene_Symbol"))
#'
#' # Build accession dictionary from a list of files
#' dict.files <- accession.dictionary.files(
#'      "accs_no_1947.txt",
#'      "accs_no_2007.txt",
#'      "accs_no_2047.txt",
#'      col.mapping = c(gene.symbol = "Gene_Symbol")
#' )
#'
#' # Build dictionary from already read-in data.frames
#' df.1947 <- read.delim("accs_no_1947.txt")
#' df.2007 <- read.delim("accs_no_2007.txt")
#' dict.data <- accession.dictionary(df.1947, df.2007,
#'                  col.mapping = c(gene.symbol = "Gene_Symbol"))
#'
accession.dictionary <- function(..., col.mapping) {
    # Validate the column mapping
    col.mapping <- .get.col.mapping(col.mapping)
    col.mapping.na <- col.mapping[!is.na(col.mapping)]

    ##
    ## Validate all data frames
    ##
    dfs <- list(...)
    dfs <- unlist(lapply(dfs, function(df) {
        if (is.list(df) && !is.data.frame(df)) {
            return(df)
        }
        return(list(df))
    }), recursive = FALSE, use.names = TRUE)

    # Check that
    #   - each item is a data.frame
    #   - each data.frame has the correct columns
    #   - each column has the correct type (we convert logical group identifiers!)
    dfs <- mapply(function(df, name) {
        if (!is.data.frame(df)) {
            stop(sprintf("Data set '%s' must be a data.frame", name))
        }

        present <- col.mapping.na %in% colnames(df)
        if (!all(present)) {
            stop(sprintf("Data set '%s' does not contain the column%s: '%s'",
                         name,
                         ifelse(sum(!present) > 1L, "s", ""),
                         paste(col.mapping.na[!present], collapse = "', '")))
        }

        if (!is.character(df[[col.mapping["accession.nr"]]]) &&
            !is.factor(df[[col.mapping["accession.nr"]]])) {
            warning(sprintf("Accession Nr. column '%s' in data set '%s' is not a character or a factor",
                            col.mapping["accession.nr"], name))
        }
        if (is.character(df[[col.mapping["accession.nr"]]])) {
            df[[col.mapping["accession.nr"]]] <- as.character(df[[col.mapping["accession.nr"]]])
        }

        # If the group identifier is a logical, we assume a "is master-gene" type
        # of column
        if (is.logical(df[[col.mapping["group.identifier"]]])) {
            if (!isTRUE(df[[col.mapping["group.identifier"]]][1L])) {
                stop("The first 'is master-gene' value in column '%s' in data set '%s' must be `TRUE`",
                     col.mapping["group.identifier"], name)
            }
            df[[col.mapping["group.identifier"]]] <- cumsum(df[[col.mapping["group.identifier"]]])
        } else if (is.character(df[[col.mapping["group.identifier"]]])) {
            df[[col.mapping["group.identifier"]]] <- as.factor(df[[col.mapping["group.identifier"]]])
        }

        return(df)
    }, dfs, names(dfs), SIMPLIFY = FALSE)

    ##
    ## Work with a smaller data set since we only need to have
    ## the group id and the accessions to build the dictionary
    ##
    all.accessions <- lapply(dfs, function(df) {
        accs <- df[[col.mapping["accession.nr"]]]
        uq <- !duplicated(accs)
        rbind(
            accs = accs[uq],
            prot = df[[col.mapping["protein.name"]]][uq],
            gene = df[[col.mapping["gene.symbol"]]][uq]
        )
    })
    all.accessions <- do.call(cbind, all.accessions)
    all.accessions <- all.accessions[ , sort.list(all.accessions["accs", ], method = "radix")]
    all.accessions <- all.accessions[ , !duplicated(all.accessions["accs", ])]

    stripped <- lapply(dfs, function(df) {
        accs <- factor(df[[col.mapping["accession.nr"]]],
                       levels = all.accessions["accs", ])
        matrix(
            c(as.integer(df[[col.mapping["group.identifier"]]]),
              as.integer(accs),
              rep.int(NA_integer_, nrow(df))),
            ncol = 3L,
            byrow = FALSE,
            dimnames = list(
                NULL,
                c("gid", "accs", "pg")
            )
        )
    })

    ##
    ## Merge groups in individual files according to overlapping Acc#s
    ##
    merged <- lapply(stripped, function(x) {
        accs.cont <- tabulate(x[ , "accs"], nbins = max(x[ , "accs"]))
        accs.ids <- which(accs.cont > 1L)
        if (length(accs.ids) > 0L) {
            for (k in accs.ids) {
                rep.n <- x[x[, "accs"] == k, "gid"]
                x[x[ , "gid"] %in% rep.n, "gid"] <- min(rep.n)
            }
        }

        x <- x[sort.list(x[ , "gid"], method = "radix"), ]
        x <- cbind(x, rank = x[ , "gid"])
        x[duplicated(x[, "rank"]), "rank"] <- 0L

        return(x)
    })

    state.env <- new.env(size = 5L)
    state.env$dict <- merged[[1L]][!duplicated(merged[[1L]][, "accs"]), ]
    state.env$pg.counter <- max(state.env$dict[ , "gid"])

    state.env$dict[, "pg"] <- state.env$dict[, "gid"]

    update.dict <- function(x, state.env) {
        for (i in seq_len(nrow(x))) {
            #check if the acc# exists already
            accs.matches <- which(state.env$dict[ , "accs"] == x[i, "accs"])
            if (length(accs.matches) == 0L) {
                # check if rank == 0, if so: it is an isoform, otherwise it is a new PG
                if (x[i, "rank"] == 0L) {
                    # (i - 1) > 0 because the first row is always a new protein, thus rank != 0
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
        # add total number of rows + 1 (CORRECT IN PREVIOUS VERSIONS!!!!) at the end
        pos.rank <- c(pos.rank, nrow(x) + 1L)

        for(m in seq_len(length(pos.rank) - 1L)) {
            # If a group is merged, then we use the minimum PG# in that group and
            # we need to change the PG# in ALL that set and in the previous version
            # of the dictionary.
            pos.seq <- seq.int(pos.rank[m], pos.rank[m + 1L] - 1L)
            new.pg <- min(x[pos.seq, "pg"])
            pg.mismatch <- which(x[pos.seq, "pg"] != new.pg)
            if (length(pg.mismatch) > 0L) {
                pg.mismatch <- x[pos.seq[pg.mismatch], "pg"]

                for (pgm in pg.mismatch) {
                    x[x[ , "pg"] == pgm, "pg"] <- new.pg

                    # change current dictionary as well: some Acc# may not be present in this
                    # set but in the dictionary, thus, the old PG# will not be changed there.
                    state.env$dict[state.env$dict[ , "pg"] == pgm, "pg"] <- new.pg
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

    lapply(merged[-1L], update.dict, state.env = state.env)

    dict <- state.env$dict
    remove(state.env)

    dict.df <- data.frame(
        dict[ , c("pg", "gid")],
        t(all.accessions[ , dict[ , "accs"]]),
        stringsAsFactors = FALSE
    )

    return(structure(list(
        dictionary = dict.df,
        col.mapping = col.mapping,
        call = match.call(expand.dots = TRUE)
    ), class = "accession.dict"))
}

#' @describeIn accession.dictionary Read in data from all files in the given directory
#' @param dir path to the directory.
#' @export
accession.dictionary.dir <- function(dir, col.mapping) {
    # Check parameter `dir`
    if (length(dir) != 1L || anyNA(dir) || !is.character(dir)) {
        stop("`dir` must be a single path")
    }

    # Check if `dir` points to an existing path
    dir <- path.expand(dir)
    if (!file.exists(dir)) {
        stop(sprintf("Directory %s does not exist", dir))
    }

    # Check if `dir` points to a populated directory
    files <- dir(dir, full.names = TRUE)
    if (length(files) == 0L) {
        stop(sprintf("Directory %s is empty", dir))
    }

    # Check if any of the items in the directory are directories itself and remove them
    dirs.in.dir <- list.dirs(dir, recursive = FALSE)
    files <- setdiff(files, dirs.in.dir)

    ret.obj <- accession.dictionary.files(files, col.mapping = col.mapping)
    ret.obj$directory <- dir
    ret.obj$call <- match.call(expand.dots = FALSE)
    return(ret.obj)
}

#' @describeIn accession.dictionary Read in data from the given files
#' @importFrom utils read.delim
#' @export
accession.dictionary.files <- function(..., col.mapping) {
    # Check all files
    files <- c(...)

    if (length(files) == 0L || anyNA(files) || !is.character(files)) {
        stop("A list of file names (characters) or vectors of file names must be supplied.")
    }

    missing.files <- !file.exists(files)
    if (sum(missing.files) == 1L) {
        stop(sprintf("File '%s' does not exist", files[missing.files]))
    } else if (sum(missing.files) > 1L) {
        stop(sprintf("Files '%s' do not exist", paste(files[missing.files], collapse = "', ")))
    }

    # Read in all files
    dfs <- lapply(files, read.delim, header = TRUE, stringsAsFactors = FALSE)

    names(dfs) <- if (!is.null(names(files))) {
        names(files)
    } else {
        files
    }

    ret.obj <- accession.dictionary(dfs, col.mapping = col.mapping)
    ret.obj$files <- files
    ret.obj$call <- match.call(expand.dots = TRUE)
    return(ret.obj)
}


#' Helper function to merge the user-supplied column mapping
#' with the default column mapping and validate the mapping.
.get.col.mapping <- function(col.mapping) {
    default.col.mapping <- c(
        group.identifier = "N",
        accession.nr = "Accession",
        protein.name = "Protein_Name",
        gene.symbol = NA_character_
    )

    if (!missing(col.mapping)) {
        override <- match(names(col.mapping), names(default.col.mapping))
        default.col.mapping[override] <- col.mapping
    }

    if (anyNA(default.col.mapping[c("group.identifier", "accession.nr", "protein.name")])) {
        stop("The columns `group.identifier`, `accession`, and `protein.name` must be present")
    }

    return (default.col.mapping)
}

