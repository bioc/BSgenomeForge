### =========================================================================
### forgeBSgenomeDataPkgFromTwobitFile()
### -------------------------------------------------------------------------
###
### Create a BSgenome data package from a 2bit file.
###


.get_seqnames_from_twobit_file <- function(filepath)
{
    if (!isSingleString(filepath) || filepath == "")
        stop(wmsg("'filepath' must be a single (non-empty) string"))
    if (!file.exists(filepath))
        stop(wmsg(filepath, ": no such file or directory"))
    twobitfile <- rtracklayer::TwoBitFile(filepath)
    ans <- try(seqlevels(twobitfile), silent=TRUE)
    if (inherits(ans, "try-error"))
        stop(wmsg("Failed to read file ", filepath, ". Is this a 2bit file?"))
    ans
}

### Same checks as check_circ_seqs() plus the "cannot be empty" check.
.check_seqnames <- function(seqnames)
{
    if (is.null(seqnames))
        return()
    if (!is.character(seqnames))
        stop(wmsg("'seqnames' must be NULL or a character vector"))
    if (length(seqnames) == 0L)
        stop(wmsg("'seqnames' cannot be empty"))
    ## "primary key" constraints (see GenomeInfoDb:::is_primary_key)
    if (anyNA(seqnames))
        stop(wmsg("'seqnames' cannot contain NA's"))
    if (!all(nzchar(seqnames)))
        stop(wmsg("'seqnames' cannot contain empty strings"))
    if (anyDuplicated(seqnames))
        stop(wmsg("'seqnames' cannot contain duplicate values"))
}

.make_pkgtitle <- function(organism, provider, genome)
{
    paste0("Full genomic sequences for ", organism,
           " (genome assembly ", genome, " from ", provider, ")")
}

.make_pkgdesc <- function(organism, provider, genome)
{
    paste0("Full genomic sequences for ", organism, " as ",
           "provided by ", provider, " (genome assembly ", genome, "). ",
           "The sequences are stored in DNAString objects.")
}

.sort_twobit_file1 <- function(origfile, destfile, seqnames)
{
    dna <- rtracklayer::import.2bit(origfile)
    m <- match(seqnames, names(dna))
    ## 'seqnames' should have been checked already so this should never happen,
    ## but just in case...
    bad_idx <- which(is.na(m))
    if (length(bad_idx) != 0L) {
        errmsg <- c("sequence name(s) in 'seqnames' not in file ",
                    origfile, ": ", paste(seqnames[bad_idx], collapse=", "))
        stop(wmsg(errmsg))
    }
    dna <- dna[m]
    rtracklayer::export.2bit(dna, destfile)
    destfile
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### forgeBSgenomeDataPkgFromTwobitFile()
###

forgeBSgenomeDataPkgFromTwobitFile <- function(filepath,
    organism, provider, genome,
    pkg_maintainer, pkg_author=NA,
    pkg_version="1.0.0",
    pkg_license="Artistic-2.0",
    seqnames=NULL,
    circ_seqs=NULL,  # must be a subset of 'seqnames' when 'seqnames' not NULL
    destdir=".")
{
    if (!isSingleString(organism) || organism == "")
        stop(wmsg("'organism' must be a single (non-empty) string"))
    if (!isSingleString(provider) || provider == "")
        stop(wmsg("'provider' must be a single (non-empty) string"))
    if (!isSingleString(genome) || genome == "")
        stop(wmsg("'genome' must be a single (non-empty) string"))
    if (!isSingleString(pkg_maintainer) || pkg_maintainer == "")
        stop(wmsg("'pkg_maintainer' must be a single (non-empty) string"))
    if (identical(pkg_author, NA)) {
        pkg_author <- pkg_maintainer
    } else if (!isSingleString(pkg_author) || pkg_author == "") {
        stop(wmsg("'pkg_author' must be a single (non-empty) string"))
    }
    if (!isSingleString(pkg_version) || pkg_version == "")
        stop(wmsg("'pkg_version' must be a single (non-empty) string"))
    if (!isSingleString(pkg_license) || pkg_license == "")
        stop(wmsg("'pkg_license' must be a single (non-empty) string"))
    .check_seqnames(seqnames)   # shallow check
    check_circ_seqs(circ_seqs)  # shallow check
    if (!isSingleString(destdir) || destdir == "")
        stop(wmsg("'destdir' must be a single (non-empty) string"))

    all_seqnames <- .get_seqnames_from_twobit_file(filepath)

    ## Deep 'seqnames' check.
    if (!is.null(seqnames)) {
        bad_seqnames <- setdiff(seqnames, all_seqnames)
        if (length(bad_seqnames) != 0L) {
            errmsg <- c("sequence name(s) in 'seqnames' not in file ",
                        filepath, ": ", paste(bad_seqnames, collapse=", "))
            stop(wmsg(errmsg))
        }
    }

    ## Deep 'circ_seqs' check.
    if (!is.null(circ_seqs)) {
        if (is.null(seqnames)) {
            valid_circ_seqs <- all_seqnames
            where <- c("file ", filepath)
        } else {
            valid_circ_seqs <- seqnames
            where <- "'seqnames'"
        }
        bad_circ_seqs <- setdiff(circ_seqs, valid_circ_seqs)
        if (length(bad_circ_seqs) != 0L) {
            errmsg <- c("sequence name(s) in 'circ_seqs' not in ", where, ": ",
                        paste(bad_circ_seqs, collapse=", "))
            stop(wmsg(errmsg))
        }
    }

    if (!is.null(seqnames)) {
        sorted_twobit_file <- file.path(tempdir(), "single_sequences.2bit")
        filepath <- .sort_twobit_file1(filepath, sorted_twobit_file, seqnames)
    }

    organism <- format_organism(organism)
    abbr_organism <- abbreviate_organism_name(organism)
    pkgname <- make_pkgname(abbr_organism, provider, genome)
    pkgtitle <- .make_pkgtitle(organism, provider, genome)
    pkgdesc <- .make_pkgdesc(organism, provider, genome)
    check_pkg_maintainer(pkg_maintainer)
    biocview <- organism2biocview(organism)
    circ_seqs <- build_Rexpr_as_string(circ_seqs)

    ## Create the package.
    origdir <- system.file("pkgtemplates", "2bit_BSgenome_datapkg",
                           package="BSgenomeForge")
    symValues <- list(BSGENOMEOBJNAME=abbr_organism,
                      PKGTITLE=pkgtitle,
                      PKGDESCRIPTION=pkgdesc,
                      PKGVERSION=pkg_version,
                      PKGAUTHOR=pkg_author,
                      PKGMAINTAINER=pkg_maintainer,
                      PKGLICENSE=pkg_license,
                      ORGANISM=organism,
                      PROVIDER=provider,
                      GENOME=genome,
                      ORGANISMBIOCVIEW=biocview,
                      CIRCSEQS=circ_seqs)
    pkg_dir <- createPackage(pkgname, destdir, origdir, symValues,
                             unlink=TRUE, quiet=FALSE)[[1]]
    if (!is.null(seqnames)) {
        move_file_to_datapkg(filepath, pkg_dir)
    } else {
        to <- file.path(pkg_dir, "inst", "extdata", "single_sequences.2bit")
        stopifnot(file.copy(filepath, to))
    }

    invisible(pkg_dir)
}

