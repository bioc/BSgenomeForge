
.get_circ_seqs_from_NCBI <- function(assembly_accession, chrominfo,
                                     organism=NULL, circ_seqs=NULL)
{
    NCBI_assembly_info <- BSgenomeForge:::.extract_NCBI_assembly_info(
                                                   assembly_accession,
                                                   chrominfo,
                                                   organism=organism,
                                                   circ_seqs=circ_seqs)
    NCBI_assembly_info$circ_seqs
}

test_that(".get_circ_seqs_from_NCBI() works for an unregistered assembly
          with NO assembled molecules", {
    ## NCBI assembly Acidianus infernus (no circular sequences).
    chrominfo <- getChromInfoFromNCBI("GCA_009729545.1")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_009729545.1", chrominfo,
                                          organism="it doesn't matter")
    expect_equal(circ_seqs, character(0))

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_009729545.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs=character(0))
    expect_equal(circ_seqs, character(0))

    regexp <- "cannot[\\s]+have[\\s]+circular[\\s]+sequences"
    expect_error(.get_circ_seqs_from_NCBI("GCA_009729545.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs="MT"),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

test_that(".get_circ_seqs_from_NCBI() works for an unregistered assembly
          with assembled molecules", {
    ## NCBI assembly ASM836960v1 (Vibrio cholerae, a g-proteobacteria).
    ## 3 sequences: "1", "2", "unnamed". All of them are circular.
    ## See CP043554.1, CP043556.1, and CP043555.1, in NCBI Nucleotide
    ## database at https://www.ncbi.nlm.nih.gov/nuccore/
    chrominfo <- getChromInfoFromNCBI("GCA_008369605.1")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_008369605.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs=c("1", "2"))
    expect_equal(circ_seqs, c("1", "2"))

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_008369605.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs="unnamed")
    expect_equal(circ_seqs, "unnamed")

    regexp <- "not[\\s]+registered"
    expect_error(.get_circ_seqs_from_NCBI("GCA_008369605.1", chrominfo,
                                          organism="it doesn't matter"),
                 regexp, ignore.case=TRUE, perl=TRUE)

    regexp <- "sequence[\\s]+names[\\s]+that[\\s]+do[\\s]+not[\\s]+belong"
    expect_error(.get_circ_seqs_from_NCBI("GCA_008369605.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs="MT"),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

test_that(".get_circ_seqs_from_NCBI() works for an unregistered assembly
          with assembled molecules and NO circular sequences", {
    ## NCBI assembly ASM369215v1 (Arthrobacter phage DrManhattan, a virus).
    ## Only 1 sequence: "NC_048093.1". It's an assembled molecule but it's
    ## not circular. See MH834610.1 in NCBI Nucleotide database at
    ## https://www.ncbi.nlm.nih.gov/nuccore/
    chrominfo <- getChromInfoFromNCBI("GCA_003692155.1")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_003692155.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs=character(0))
    expect_equal(circ_seqs, character(0))

    regexp <- "not[\\s]+registered"
    expect_error(.get_circ_seqs_from_NCBI("GCA_003692155.1", chrominfo,
                                          organism="it doesn't matter"),
                 regexp, ignore.case=TRUE, perl=TRUE)

    regexp <- "sequence[\\s]+names[\\s]+that[\\s]+do[\\s]+not[\\s]+belong"
    expect_error(.get_circ_seqs_from_NCBI("GCA_003692155.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs="MT"),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

test_that(".get_circ_seqs_from_NCBI() throws an error when both assembled
          and unassembled molecules are specified", {
    ## NCBI assembly ASM2578275v1 (Saccharomyces cerevisiae, baker's yeast).
    ## Only 1 circular sequence, "MT". See CP088317.1 in NCBI Nucleotide
    ## database at https://www.ncbi.nlm.nih.gov/nuccore/
    chrominfo <- getChromInfoFromNCBI("GCA_025782755.1")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_025782755.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs="MT")
    expect_equal(circ_seqs, "MT")

    regexp <- "must[\\s]+be[\\s]+names[\\s]+of[\\s]+assembled[\\s]+molecules"
    expect_error(.get_circ_seqs_from_NCBI("GCA_025782755.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs="CP088301.1"),
                 regexp, ignore.case=TRUE, perl=TRUE)

    regexp <- "must[\\s]+be[\\s]+names[\\s]+of[\\s]+assembled[\\s]+molecules"
    expect_error(.get_circ_seqs_from_NCBI("GCA_025782755.1", chrominfo,
                                          organism="it doesn't matter",
                                          circ_seqs=c("MT", "CP088301.1")),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

test_that(".get_circ_seqs_from_NCBI() works for a registered assembly
          with one assembled molecule", {
    ## NCBI assembly Felis_catus_9.0 (1 circular sequence: "MT").
    chrominfo <- getChromInfoFromNCBI("GCF_000181335.3")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCF_000181335.3", chrominfo)
    expect_equal(circ_seqs, "MT")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCF_000181335.3", chrominfo,
                                          circ_seqs="MT")
    expect_equal(circ_seqs, "MT")

    regexp <- "circular[\\s]+sequences[\\s]+are[\\s]+known"
    expect_error(.get_circ_seqs_from_NCBI("GCF_000181335.3", chrominfo,
                                          circ_seqs="chrA1"),
                 regexp, ignore.case=TRUE, perl=TRUE)

    regexp <- "circular[\\s]+sequences[\\s]+are[\\s]+known"
    expect_error(.get_circ_seqs_from_NCBI("GCF_000181335.3", chrominfo,
                                          circ_seqs=c("chrA1", "MT")),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

test_that(".get_circ_seqs_from_NCBI() works for a registered assembly
          with assembled molecules and more than one circular sequence", {
    ## NCBI assembly for Vitis vinifera (2 circular sequences: "MT"
    ## and "Pltd").
    chrominfo <- getChromInfoFromNCBI("GCF_000003745.3")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCF_000003745.3", chrominfo)
    expect_equal(circ_seqs, c("MT", "Pltd"))

    circ_seqs <- .get_circ_seqs_from_NCBI("GCF_000003745.3", chrominfo,
                                          circ_seqs=c("MT", "Pltd"))
    expect_equal(circ_seqs, c("MT", "Pltd"))

    regexp <- "circular[\\s]+sequences[\\s]+are[\\s]+known"
    expect_error(.get_circ_seqs_from_NCBI("GCF_000003745.3", chrominfo,
                                          circ_seqs="1"),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

test_that(".get_circ_seqs_from_NCBI() works for a registered assembly
          with NO circular sequences", {
    ## NCBI assembly for Betacoronavirus (no circular sequences).
    chrominfo <- getChromInfoFromNCBI("GCA_009948555.1")

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_009948555.1", chrominfo)
    expect_equal(circ_seqs, character(0))

    circ_seqs <- .get_circ_seqs_from_NCBI("GCA_009948555.1", chrominfo,
                                          circ_seqs=character(0))
    expect_equal(circ_seqs, character(0))

    regexp <- "no[\\s]+known[\\s]+circular[\\s]+sequences"
    expect_error(.get_circ_seqs_from_NCBI("GCA_009948555.1", chrominfo,
                                          circ_seqs="MT"),
                 regexp, ignore.case=TRUE, perl=TRUE)

    regexp <- "no[\\s]+known[\\s]+circular[\\s]+sequences"
    expect_error(.get_circ_seqs_from_NCBI("GCA_009948555.1", chrominfo,
                                          circ_seqs="MN996531.1"),
                 regexp, ignore.case=TRUE, perl=TRUE)
})

