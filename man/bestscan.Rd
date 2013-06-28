\name{make.bestscan}
\alias{make.bestscan}
\title{Make 1-D scans from best patterns}
\description{
  Takes best model per trait and builds a matrix of 1-D scans.
}
\usage{
make.bestscan(cross, best, pattern, scan.type = "BIM",
  cross.name = deparse(substitute(cross)),
  maps = get(paste(cross.name, "maps", sep = ".")),
  trait.annotation = get(paste(cross.name, "traitnames", sep = ".")),
  attenuate = 1, tolerance = 1e-05, verbose = TRUE)
}
\arguments{
  \item{cross}{Object created by \code{\link[qtl]{read.cross}}.}
  \item{best}{Data frame with trait name as first column and
    \code{chrom} and \code{locus} as two other columns. Last column
    contains the QTL-specific score (e.g. explained
    \code{variance}). See details.}
  \item{pattern}{Data frame with trait name as first column,
    and columns for \code{pattern} and \code{posterior}.}
  \item{scan.type}{Type of scan as character string.}
  \item{cross.name}{Name of cross object as character string.}
  \item{maps}{Object of class \code{\link{read.maps}}.}
  \item{trait.annotation}{Annotation object for traits.}
  \item{attenuate}{Degree of attenuation of peaks.}
  \item{tolerance}{Tolerance level (drop scans below this level).}
  \item{verbose}{Show steps if \code{TRUE}.}
}
\details{
  This converts results from an object that has the best pattern per trait
  to a matrix with 1-D scans, one trait per row. This is meant to be
  saved for use in multtrait, ideally filtered through a web page.

  The \code{best} object has one row per 
  QTL for each trait; traits with no QTL can be included with \code{NA}
  for \code{chrom} and \code{locus}. Other columns in best will depend on
  the method of construction, but the last column should contain a
  QTL-specific score, which is used for plotting. For instance, using
  qb.best from suggested package R/qtlbim includes a column for the posterior
  probability for that QTL (\code{n.qtl}), and a last column with the
  QTL-explained \code{variance}.
  
  The \code{pattern} object is constructed as best as possible from
  \code{best} if not provided. The \code{trait.annotation}, created with
  \code{\link{read.annotation}}, is used for trait positions if provided
  (see \code{\link{aug.scanone}}).
  
  The \code{attenuate} determines the degree of attenuation down from
  peak. Default is an exponential decay corresponding to genetic linkage
  in cM.
}
\author{Brian S. Yandell}
\seealso{\code{\link[qtl]{read.cross}}, \code{\link{read.maps}},
  \code{\link{read.annotation}}}
\examples{
\dontrun{
tissue.bestscan <- make.bestscan(cross, tissue.best, tissue.pattern)
}
}
\keyword{utilities}
