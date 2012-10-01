\name{cistrans}
\alias{cistrans}
\alias{plot.cistrans}
\alias{print.cistrans}
\alias{summary.cistrans}
\title{Create and plot cis-trans map}
\description{
  Creates cistrans object that can be plotted with different color schemes
  to accentuate strenght of LOD peak and proximity to physical location.
}
\usage{
cistrans(x, filename, peak.chr, trans.chr, min.lod = 0,
  main, sep = "\t", use.density = TRUE, \dots)
\method{plot}{cistrans}(x, cex, main = attr(x, "main"),
  col.score = c("score","cis","trans"), cis.only,
  xlim, ylim, use.cM = FALSE, jitter = FALSE,
  peak.chr = attr(x, "peak.chr"), trans.chr = attr(x, "trans.chr"), \dots)
\method{summary}{cistrans}(object, \dots)
\method{print}{cistrans}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{data frame with LOD peak for traits (rows) and chromosomes (columns)}
  \item{filename}{name of file with \code{x} if it is to be read in}
  \item{cross.name}{name of cross (used by \code{maps} and \code{trait.annotation})}
  \item{maps}{list of \code{\link[qtl]{pull.map}} objects created by
    \code{\link{read.maps}}; default is associated with \code{cross.name}}
  \item{peak.chr}{names of chromosome(s) to show peaks}
  \item{trans.chr}{names of chromosome(s) to include traits}
  \item{min.lod}{minimum LOD to keep}
  \item{main}{title of plot (default taken from levels of \code{x$Tissue}}
  \item{sep}{separation character if reading in from \code{filename}}
  \item{use.density}{use density in plot of cumulative score if \code{TRUE}}
  \item{trait.annotation}{trait annotation data frame; default is
    associated with \code{cross.name}}
  \item{use.annot}{use annotation if \code{TRUE}}
  \item{\dots}{additional arguments passed to some subroutines}
}
\details{
  The \code{trait.annotation} has \code{a_gene_id} column that is matched
  with \code{x$a.gene.id} to then use \code{Chromosome} and
  \code{Chromosome_Position} (in Mb) for plotting.
  }
\value{
  Data frame of class \code{cistrans} with added components for summary
  and plot methods.
}
\author{Brian S. Yandell}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
}
\keyword{utilities}
