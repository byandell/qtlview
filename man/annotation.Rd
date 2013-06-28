\name{read.annotation}
\alias{read.annotation}
\title{Create trait.annotation object}
\description{
Creates a trait.annotation object that associates a_gene_id with Symbol.
}
\usage{
read.annotation(filename, update.names = NULL, drop.extra = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{File name where gene annotations are found.}
  \item{update.names}{Named character vector to translate names in file to
    \code{trait.annotation} column names if not \code{NULL}. See details.}
  \item{drop.extra}{Drop extra names in file if \code{TRUE}.}
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
The file is read and any name changes are performed. The names of
\code{update.names} must correspond (partial match) to columns in the
returned object. The values of \code{update.names} are the (partial)
match to column names in the file identified by \code{filename}.

The \code{Chromosome_Position} column, if missing, is created by
averaging the\code{Start Coordinate} and \code{End Coordinate} columns
and rescaling from base pairs to Mb.
}
\value{
A data frame with entries:
    \item{\code{a_substance_id}}{Substance ID number.}
    \item{\code{accession_code}}{Accession code.}
    \item{\code{a_gene_id}}{Gene ID number.}
    \item{\code{Symbol}}{Gene symbol.}
    \item{\code{Alternate_Symbols}}{Alternate gene symbols.}
    \item{\code{Name}}{Gene name.}
    \item{\code{Alternate_Names}}{Alternate gene names, comma separated.}
    \item{\code{Description}}{Gene description.}
    \item{\code{Chromosome}}{Chromosome name.}
    \item{\code{Strand}}{Strand direction (+1 or 1).}
    \item{\code{Chromosome_Position}}{Chromosome position in Mb.}
}
\author{Brian S. Yandell (byandell@wisc.edu)}
\examples{
\dontrun{
## Do this once when creating your cross library.
my.file.name <- "~/Rlib/eqtl/data/Rosetta/F2.Annotation.mm8.20080208.csv"
my.data.lib <- "~/Rlib/eqtl/B6BTBR07"

trait.annotation <- read.annotation(my.file.name,
  update.names = c(Symbol = "gene1", Alternate_Symbols = "gene",
    Name = "Gene_Name", Alternate_Names = "AlternateNames",
    Description = "ProteomeShortDescription",
    Chromosome = "CHROMOSOME", Start_Coordinate = "STARTCOORD",
    End_Coordinate = "ENDCOORD", Strand = "STRAND"))

save("trait.annotation",
  file = file.path(my.data.lib, "data", "trait.annotation.RData"),
  compress = TRUE)
}
}
\keyword{utilities}
