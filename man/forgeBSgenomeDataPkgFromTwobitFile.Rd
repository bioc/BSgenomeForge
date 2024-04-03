\name{forgeBSgenomeDataPkgFromTwobitFile}

\alias{forgeBSgenomeDataPkgFromTwobitFile}

\title{Create a BSgenome data package from a 2bit file}

\description{
  The \code{forgeBSgenomeDataPkgFromTwobitFile} function allows the user
  to create a BSgenome data package from a 2bit file.
}

\usage{
forgeBSgenomeDataPkgFromTwobitFile(filepath, organism, provider, genome,
             pkg_maintainer, pkg_author=NA,
             pkg_version="1.0.0", pkg_license="Artistic-2.0",
             seqnames=NULL,
             circ_seqs=NULL,
             destdir=".")
}

\arguments{
  \item{filepath}{
    A single string containing a path to a 2bit file.
  }
  \item{organism}{
    The full name of the organism e.g. \code{"Escherichia coli"},
    \code{"Mus musculus"}, \code{"Saccharomyces cerevisiae"}, etc...
  }
  \item{provider}{
    The provider of the sequence data stored in the 2bit file e.g.
    \code{"NCBI"}, \code{"UCSC"}, \code{"BDGP"}, \code{"FlyBase"},
    \code{"WormBase"}, etc... Will be used to make part 3 of the package name.
  }
  \item{genome}{
    A single string specifying the name of the genome assembly
    e.g. \code{"ASM1484v1"}, \code{"mm39"}, or \code{"sacCer3"}.
    Will be used to make part 4 of the package name.
  }
  \item{pkg_maintainer}{
    A single string containing the name and email address of the package
    maintainer (e.g \code{"Jane Doe, <janedoe@gmail.com>"}).
  }
  \item{pkg_author}{
    A single string containing the name of the package author. When
    unspecified, this takes the value of \code{pkg_maintainer}.
  }
  \item{pkg_version}{
    The version of the package. Set to \code{"1.0.0"} by default.
  }
  \item{pkg_license}{
    The license of the package. This must be the name of a software license
    used for free and open-source packages. Set to \code{"Artistic-2.0"}
    by default.
  }
  \item{seqnames}{
    NULL (the default), or a character vector containing a subset of the
    sequence names stored in the 2bit file. Use this to select and/or
    reorder the sequences that will go in the BSgenome data package.

    By default (i.e. when \code{seqnames} is NULL), all the sequences in
    the file are selected and put in the BSgenome data package in the order
    that they appear in the 2bit file.

    Use \code{seqlevels(rtracklayer::TwoBitFile(filepath))} to see the
    sequence names stored in the 2bit file.
  }
  \item{circ_seqs}{
    NULL (the default), or a character vector providing the names of the
    sequences stored in the 2bit file that are known to be circular.
    Note that if \code{seqnames} is supplied then \code{circ_seqs} must
    be a subset of it.

    Set to \code{circ_seqs} to \code{character(0)} if none of the sequences
    in the 2bit file (or in \code{seqnames}) is known to be circular.
    
    By default (i.e. if \code{circ_seqs} is NULL), then the circular sequences
    are guessed based on their names. Because this is not really reliable,
    a warning is issued.
  }
  \item{destdir}{
    A single string containing the path to the directory where the
    BSgenome data package is to be created. This directory must already
    exist. Note that, by default, the package will be created in the
    current directory (\code{"."}).
  }
}

\details{
  This function is intended for Bioconductor users who want to forge a BSgenome
  data package from a 2bit file.
}

\value{
  The path to the created package as an invisible string.
}

\author{
  Hervé Pagès
}

\seealso{
  \itemize{
    \item The \code{\link{forgeBSgenomeDataPkgFromNCBI}} function for
          creating a BSgenome data package from an NCBI assembly
          (similar to \code{forgeBSgenomeDataPkgFromTwobitFile} but
          slightly more convenient).

    \item The \code{\link{forgeBSgenomeDataPkgFromUCSC}} function for
          creating a BSgenome data package from a UCSC genome
          (similar to \code{forgeBSgenomeDataPkgFromTwobitFile} but
          slightly more convenient).

    \item The \code{\link{fastaTo2bit}} function to convert a FASTA
          file to the 2bit format.
  }
}

\examples{
## Download the FASTA file containing the genomic sequences for
## for Escherichia coli assembly ASM1484v1 (see
## https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000014845.1/):
fasta_path <- downloadGenomicSequencesFromNCBI("GCF_000014845.1",
                                               destdir=tempdir())

## Use fastaTo2bit() to convert the file to 2bit:
twobit_path <- tempfile(fileext=".2bit")
fastaTo2bit(fasta_path, twobit_path)

## All the DNA sequences in Escherichia coli are circular:
circ_seqs <- seqlevels(rtracklayer::TwoBitFile(twobit_path))

## Note that seqlevels(rtracklayer::TwoBitFile(.)) is equivalent to
## names(rtracklayer::import.2bit(.)) but a lot more efficient because
## it doesn't load the sequences.

## Create a BSgenome data package from the 2bit file:
forgeBSgenomeDataPkgFromTwobitFile(
    filepath=twobit_path,
    organism="Escherichia coli",
    provider="NCBI",
    genome="ASM1484v1",
    pkg_maintainer="Jane Doe <janedoe@gmail.com>",
    circ_seqs=circ_seqs,
    destdir=tempdir()
)
}

\keyword{utilities}
