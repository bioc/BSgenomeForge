\name{forgeBSgenomeDataPkgFRomNCBI}

\alias{forgeBSgenomeDataPkgFromNCBI}

\title{Create a BSgenome data package from an NCBI assembly}

\description{
  The \code{forgeBSgenomeDataPkgFRomNCBI} function allows the user to
  create a BSgenome data package from an NCBI assembly.
}

\usage{
forgeBSgenomeDataPkgFromNCBI(assembly_accession,
             pkg_maintainer, pkg_author=NA,
             pkg_version="1.0.0", pkg_license="Artistic-2.0",
             organism=NULL, circ_seqs=NULL, destdir=".")
}

\arguments{
  \item{assembly_accession}{
    A single string containing a GenBank assembly accession (e.g.
    \code{"GCA_009729545.1"}) or a RefSeq assembly accession (e.g.
    \code{"GCF_000857545.1"}). Alternatively, if the assembly
    is registered in the \pkg{GenomeInfoDb} package (see
    \code{?GenomeInfoDb::\link{registered_NCBI_assemblies}}), the
    assembly name (e.g. \code{"mLoxAfr1.hap2"}) can be supplied
    instead of its GenBank or RefSeq accession.
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
  \item{organism}{
    The full name of the organism e.g. \code{"Homo sapiens"},
    \code{"Felis catus"}, \code{"Loxodonta africana"},
    \code{"Acidianus infernus"}, etc... Only needs to be specified
    if the assembly is \emph{not} registered in the \pkg{GenomeInfoDb}
    package (see \code{?GenomeInfoDb::\link{registered_NCBI_assemblies}}).
  }
  \item{circ_seqs}{
    NULL (the default), or a character vector containing the names of the
    circular sequences in the assembly. This only needs to be specified if
    the assembly is \emph{not} registered in the \pkg{GenomeInfoDb} package
    (if the assembly is registered then its circular sequences are known so
    there's no need to specify \code{circ_seqs}).

    Notes:
    \itemize{
      \item You can use \code{registered_NCBI_assemblies()} to get the list
            of assemblies that are registered in the \pkg{GenomeInfoDb}
            package.
      \item Only assembled molecules can be circular. To see the
            list of assembled molecules for a given assembly, call
            \code{getChromInfoFromNCBI(assembly_accession,
            assembled.molecules.only=TRUE)$SequenceName}.
      \item If the assembly is not registered and does not have
            circular sequences, then \code{circ_seqs} must be set
            to \code{character(0)}.
    }
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
  data package from an NCBI assembly. It typically makes use of the
  \code{downloadGenomicSequencesFromNCBI} utility function to download the
  compressed FASTA file that contains the genomic sequences of the assembly,
  and stores it in the working directory. However, if the file already
  exists in the working directory, then it is used and not downloaded again.
}

\value{
  The path to the created package as an invisible string.
}

\author{
  Atuhurira Kirabo Kakopo
}

\seealso{
  \itemize{
    \item The \code{\link[GenomeInfoDb]{registered_NCBI_assemblies}} and
          \code{\link[GenomeInfoDb]{getChromInfoFromNCBI}} functions
          defined in the \pkg{GenomeInfoDb} package.

    \item The \code{\link{downloadGenomicSequencesFromNCBI}} function that
          \code{forgeBSgenomeDataPkgFromNCBI} uses internally to download
          the genomic sequences from NCBI.

    \item The \code{\link{fastaTo2bit}} function that
          \code{forgeBSgenomeDataPkgFromNCBI} uses internally to convert
          the file downloaded by \code{downloadGenomicSequencesFromNCBI}
          from FASTA to 2bit.

    \item The \code{\link{forgeBSgenomeDataPkgFromUCSC}} function for
          creating a BSgenome data package from a UCSC genome.
  }
}

\examples{
## ---------------------------------------------------------------------
## EXAMPLE 1
## ---------------------------------------------------------------------

## Create a BSgenome data package for NCBI assembly GCA_009729545.1
## (organism Acidianus infernus):
forgeBSgenomeDataPkgFromNCBI(assembly_accession="GCA_009729545.1",
                             pkg_maintainer="Jane Doe <janedoe@gmail.com>",
                             organism="Acidianus infernus",
                             destdir=tempdir())

## ---------------------------------------------------------------------
## EXAMPLE 2
## ---------------------------------------------------------------------

## Create a BSgenome data package for NCBI assembly GCF_000857545.1
## (organism Torque teno virus 1):
forgeBSgenomeDataPkgFromNCBI(assembly_accession="GCF_000857545.1",
                             pkg_maintainer="Jane Doe <janedoe@gmail.com>",
                             organism="Torque teno virus 1",
                             circ_seqs="NC_002076.2",
                             destdir=tempdir())
}

\keyword{utilities}
