\name{add.phenos}
\alias{add.phenos}
\title{Create genetic and physical list of maps object.}
\description{
Add phenotypes to cross object.
}
\usage{
add.phenos(cross, newdata = NULL, index = NULL)
}
\arguments{
  \item{cross}{Objects of class \code{cross}. See
    \code{\link[qlt]{read.cross}}.}
  \item{newdata}{Data frame with one phenotype per column. Row names should correspond to \code{index} in \code{cross} object.}
  \item{index}{Numerical index or character name of phenotype in \code{cross} that corresponds to row names of \code{newdata}.}
}
\details{
Returns the \code{cross} object with added phenotypes.
No change if \code{newdata} is \code{NULL}.
Assumes \code{newdata} are in same individual order as \code{cross} if \code{index} is \code{NULL}; otherwise, matches row names of \code{newdata} to the phenotype identified by \code{index}.
}
\seealso{
\code{\link{snp.record}}
}
\examples{
data(hyper)
newdata <- data.frame(x = rnorm(nind(hyper)))
cross <- add.phenos(hyper, newdata)
summary(cross)

## Use of index. First need to create one.
newdata <- data.frame(myindex = seq(nind(cross)))
cross <- add.phenos(cross, newdata)
newdata <- data.frame(y = rnorm(100))
row.names(newdata) <- sample(seq(nind(cross)), 100)
cross <- add.phenos(cross, newdata)
}
\keyword{utilities}
