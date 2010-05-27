\name{listof.complod}
\alias{listof.complod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Composite LOD plot across multiple groups}
\description{
Run multtrait and composite LOD on a single plot.
}
\usage{
listof.complod(file.groups, dir = NULL, scan.type = "LOD", ...)
\method{plot}{listof.complod}(x, lodcolumn = 1, ylim, ylab,
  col, lty, pre.main = "", main, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file.groups}{Vector of character strings with file names to be
    used for \code{\link{multtrait}}.}
  \item{dir}{ ~~Describe \code{dir} here~~ }
  \item{scan.type}{Scan type: \code{LOD} and \code{LPD} use
    \code{support.lod} in \code{ylab}.}
  \item{x}{Object of class \code{listof.complod} created by listof.complod.}
  \item{lodcolumn}{Value 1 or 2 corresponding to composite LOD or (if
    present) support interval.}
  \item{ylim,ylab,col,lty}{Plot arguments}
  \item{pre.main}{Prefix for default main label}
  \item{main}{Main label: default has pre.main followed by color legend.}
  \item{\dots}{Additional arguments for routines called.}
}
\details{
  This is a kludge. Current use is to use the web page for multtrait to
  create the text files, then use those files again, manually, to create
  a combined plot. Text files are saved in
  \url{http://attie.wisc.edu/lab/data/data_for_intranet/r_scripts}.
}
\author{Brian S. Yandell}
\seealso{\code{\link{multtrait}}}
\examples{
\dontrun{
## Extract data for multtrait as tab-delimited files.
## Now merge the plots.
tmp <- listof.complod(c(a="scana.txt", b="scanb.txt"))
plot(tmp)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
