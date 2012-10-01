\name{read.perms}
\alias{read.perms}
\alias{summary.read.perms}
\alias{plot.read.perms}
\title{Create list of scanone perm objects.}
\description{
Use cross and file of physical locations to create a list of R/qtl map objects.
}
\usage{
read.perms(cross, filenames, row.names = NULL, perms,
  method = "hk", model = "normal", \dots)
\method{summary}{read.perms}(object, \dots)
\method{plot}{read.perms}(x, \dots)
}
\arguments{
  \item{cross}{Objects of class \code{cross}. See
    \code{\link[qlt]{read.cross}}.}
  \item{filenames}{Character vector of file names containing
    permutations. Needs to have names from \code{c("both","female","male").}
    Files should be CSV saved from call to \code{\link[qtl]{scanone}} with a large
    value (e.g. 10000) of \code{n.perm}. Usually will have two columns
    labeled \code{A} and \code{X} for autosomes and X chromosome, resp.}
  \item{row.names}{Argument for \code{\link[utils]{read.csv}}. Should be
    NULL if CSV file saved without row names, otherwise set to 1.}
  \item{perms}{Named list of data frames containing permutations. Each
    element should be a dataframe or \code{scanoneperm} object from
    \code{\link[qtl]{scanone}}. Names must be from
    \code{c("both","female","male")}.}
  \item{object,x}{Object of class \code{read.perms}.}
  \item{\dots}{Arguments to pass to \code{plot} or \code{summary} method.}
}
\details{
This collects together permutations for male, female and both sexes
combined to be used for high throughput QTL scans. The data frames by
sex are coerced to \code{scanoneperm} format.
}
\seealso{
\code{\link[B6BTBR07]{B6BTBR07.perms}}
}
\examples{
\dontrun{
example(B6BTBR07.perms)
}
}
\keyword{utilities}
