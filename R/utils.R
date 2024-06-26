### =========================================================================
### Some low-level utilities
### -------------------------------------------------------------------------
###
### Nothing in this file is exported.
###


### TODO: Move this to S4Vectors (or BiocBaseUtils).
load_package_gracefully <- function(package, ...)
{
    if (!requireNamespace(package, quietly=TRUE))
        stop("Could not load package ", package, ". Is it installed?\n\n  ",
             wmsg("Note that ", ..., " requires the ", package, " package. ",
                  "Please install it with:"),
             "\n\n    BiocManager::install(\"", package, "\")")
}

### If 'collapse' is FALSE (the default), returns a logical vector parallel
### to 'x'. Otherwise, returns a single logical TRUE or FALSE.
### TODO: Move this to Biostrings (put it in same file as replaceAmbiguities).
hasAmbiguities <- function(x, other.than.N=FALSE, collapse=FALSE)
{
    stopifnot(is(x, "DNAStringSet"),
              isTRUEorFALSE(other.than.N),
              isTRUEorFALSE(collapse))
    ambig_letters <- setdiff(names(IUPAC_CODE_MAP), DNA_BASES)
    if (other.than.N)
        ambig_letters <- setdiff(ambig_letters, "N")
    af <- alphabetFrequency(x, collapse=collapse)
    if (collapse) {
        ## 'af' is a named integer vector.
        stopifnot(all(ambig_letters %in% names(af)))
        ambig_count <- sum(af[ambig_letters])
    } else {
        ## 'af' is an integer matrix with 1 row per sequence in 'x'.
        stopifnot(all(ambig_letters %in% colnames(af)))
        ambig_count <- rowSums(af[ , ambig_letters, drop=FALSE])
        if (!is.null(names(x)))
            names(ambig_count) <- names(x)
    }
    ambig_count != 0L
}

### Returns TRUE if 'assembly_accession' is a GenBank accession, or FALSE
### if it's a RefSeq accession, or an error if it's none.
is_GenBank_accession <- function(assembly_accession)
{
    stopifnot(isSingleString(assembly_accession))
    if (assembly_accession == "")
        stop(wmsg("'assembly_accession' cannot be an empty string"))
    ## Make sure that 'assembly_accession' looks either like a GenBank
    ## or like a RefSeq accession. We only look at its first 5 characters.
    if (grepl("^GCA_[0-9]", assembly_accession))
        return(TRUE)
    if (grepl("^GCF_[0-9]", assembly_accession))
        return(FALSE)
    stop(wmsg("malformed assembly accession: ", assembly_accession))
}

drop_rows_with_NA_accns <- function(chrominfo, accession_col)
{
    accns <- chrominfo[ , accession_col]
    drop_idx <- which(is.na(accns))
    if (length(drop_idx) != 0L)
        chrominfo <- S4Vectors:::extract_data_frame_rows(chrominfo, -drop_idx)
    chrominfo
}

### 'character_class' must be a string containing a set of character
### specifications like in [0-9\\s] but **without** the square brackets
### e.g. "0-9\\s".
### The function is vectorized with respect to 'x'.
### Returns a 2-col matrix with one row per string in 'x'.
### Sanity check:
###   x <- c("", " ", "a", "2", "  ", " a", " 2", "a ", "aa", "a2", "2 ",
###          "2a", "22", "a23", "a 3", "a3 ", "2 3 ", " a  1 bb 2  33 ")
###   prefix_suffix <- .split_suffix(x, "0-9\\s")
###   reconstructed <- paste0(prefix_suffix[ , 1], prefix_suffix[ , 2])
###   stopifnot(identical(x, reconstructed))
.split_suffix <- function(x, character_class)
{
    pattern <- sprintf("^.*[^%s]([%s]*)$", character_class, character_class)
    suffix <- sub(pattern, "\\1", x, perl=TRUE)
    prefix <- substr(x, 1L, nchar(x) - nchar(suffix))
    colnames <- c("prefix", "suffix")
    matrix(c(prefix, suffix), ncol=2L, dimnames=list(NULL, colnames))
}

format_organism <- function(organism)
{
    ## Remove leading and trailing whitespaces (like strip() in Python).
    organism <- gsub("^\\s*|\\s*$", "", organism)
    ## Replace multiple whitespaces with single space.
    organism <- gsub("\\s+", " ", organism)
    first_letter <- substring(organism, 1, 1)
    other_letters <- substring(organism, 2)
    paste0(toupper(first_letter), tolower(other_letters))
}

abbreviate_organism_name <- function(organism)
{
    ## If 'organism' has a numeric suffix, we handle it separately.
    ## This is an attempt at doing something sensible with organism names
    ## like "Torque teno virus 1" where we want .abbreviate_organism_name()
    ## to return "Tvirus1".
    ## Note that we allow whitespaces in the numeric suffix but we'll remove
    ## them before adding the suffix back to the abbreviated organism name.
    prefix_suffix <- .split_suffix(organism, "0-9\\s")
    organism <- prefix_suffix[ , "prefix"]

    ## Abbreviate 'organism' e.g.
    ##   "Homo sapiens" -> "Hsapiens"
    ##   "Canis lupus familiaris" -> "Cfamiliaris"
    ##   "Torque teno virus" -> "Tvirus"
    parts <- strsplit(organism, "\\s+")[[1]]
    if (length(parts) <= 1) {
        abbr_organism <- parts
    } else {
        first_letter <- substr(head(parts, 1), 1, 1)
        last_part <- tail(parts, 1)
        abbr_organism <- paste0(first_letter, last_part)
    }

    ## Remove whitespaces from the numeric suffix.
    suffix <- gsub("\\s", "", prefix_suffix[ , "suffix"])

    ## Add numeric suffix back.
    paste0(abbr_organism, suffix)
}

check_pkg_maintainer <- function(pkg_maintainer)
{
    if (!(isSingleString(pkg_maintainer) && nzchar(pkg_maintainer)))
        stop(wmsg("'pkg_maintainer' must be a single (non-empty) string"))
    pattern <- "\\<[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,}\\>"
    if (!grepl(pattern, pkg_maintainer, ignore.case=TRUE))
        stop(wmsg("please provide a valid email address"))
}

make_pkgname <- function(abbr_organism, provider, genome)
{
    part4 <- gsub("[^0-9a-zA-Z.]", "", genome)
    paste("BSgenome", abbr_organism, provider, part4, sep=".")
}

organism2biocview <- function(organism)
{
    parts <- strsplit(organism, " +")[[1]]
    first_part <- head(parts, 1)
    last_part <- tail(parts, 1)
    paste0(first_part, "_", last_part)
}

build_Rexpr_as_string <- function(seqnames)
{
    stopifnot(is.character(seqnames))
    if (length(seqnames) == 0L)
        return("character(0)")
    paste0('c', '(', paste0('"', seqnames, '"', collapse=","), ')')
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helper functions used by .get_circ_seqs_from_NCBI() and
### .get_circ_seqs_from_UCSC()
###

check_circ_seqs <- function(circ_seqs)
{
    if (is.null(circ_seqs))
        return()
    if (!is.character(circ_seqs))
        stop(wmsg("'circ_seqs' must be NULL or a character vector"))
    ## "primary key" constraints (see GenomeInfoDb:::is_primary_key)
    if (anyNA(circ_seqs))
        stop(wmsg("'circ_seqs' cannot contain NA's"))
    if (!all(nzchar(circ_seqs)))
        stop(wmsg("'circ_seqs' cannot contain empty strings"))
    if (anyDuplicated(circ_seqs))
        stop(wmsg("'circ_seqs' cannot contain duplicate values"))
}

get_circ_seqs_for_registered_assembly_or_genome <-
    function(assembly_or_genome, seqnames, is_circular,
             circ_seqs=NULL, what="assembly")
{
    known_circ_seqs <- seqnames[is_circular]
    if (is.null(circ_seqs))
        return(known_circ_seqs)
    if (setequal(circ_seqs, known_circ_seqs))
        return(circ_seqs)
    msg <- c("This ", what, " is registered in the GenomeInfoDb package ")
    if (length(known_circ_seqs) == 0L) {
        msg <- c(msg, "and it has no known circular sequences.")
    } else {
        in1string <- paste0("\"", known_circ_seqs, "\"", collapse=", ")
        msg <- c(msg, "which means that its circular sequences are known, ",
                      "so you are not required to specify them. However, ",
                      "if you do specify them, then they must match the ",
                      "known ones. The known circular sequences for ",
                      "registered ", what, " ", assembly_or_genome, " ",
                      "are: ", in1string)
    }
    stop(wmsg(msg))
}

.assembly_has_no_assembled_molecules <- function(circ_seqs, what)
{
    if (length(circ_seqs) == 0L)
        return(character(0))
    stop(wmsg("This ", what, " contains no assembled molecules ",
              "so cannot have circular sequences."))
}

### 'is_assembled' will be set to:
### - a logical vector parallel to 'seqnames' when the function is called
###   by .get_circ_seqs_from_NCBI();
### - NULL when the function is called by .get_circ_seqs_from_UCSC().
get_circ_seqs_for_unregistered_assembly_or_genome <-
    function(assembly_or_genome, seqnames, is_assembled,
             circ_seqs=NULL, what="assembly")
{
    if (!is.null(is_assembled) && !any(is_assembled))
        return(.assembly_has_no_assembled_molecules(circ_seqs, what))
    if (is.null(circ_seqs))
        stop(wmsg("This ", what, " is not registered in the ",
                  "GenomeInfoDb package so I don't know what its ",
                  "circular sequences are (if any). Please provide ",
                  "their names in a character vector passed to ",
                  "the 'circ_seqs' argument (set 'circ_seqs' to ",
                  "character(0) if the ", what, " has no circular ",
                  "sequences)."))
    ## The sequence names in 'circ_seqs' must belong to the assembly or genome.
    if (anyNA(match(circ_seqs, seqnames)))
        stop(wmsg("'circ_seqs' contains sequence names that ",
                  "do not belong to the specified ", what, " ",
                  "(", assembly_or_genome, ")"))
    if (!is.null(is_assembled)) {
        ## The sequence names in 'circ_seqs' must be names of **assembled**
        ## molecules.
        assembled_molecules <- seqnames[is_assembled]
        if (!all(circ_seqs %in% assembled_molecules))
            stop(wmsg("all the sequence names in 'circ_seqs' must be ",
                      "names of assembled molecules"))
    }
    circ_seqs
}

