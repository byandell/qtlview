\name{mybim}
\alias{mybim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mybim(traitname, log10 = FALSE, n.iter = 3000, cross = B6BTBR07, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traitname}{ ~~Describe \code{traitname} here~~ }
  \item{log10}{ ~~Describe \code{log10} here~~ }
  \item{n.iter}{ ~~Describe \code{n.iter} here~~ }
  \item{cross}{ ~~Describe \code{cross} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (traitname, log10 = FALSE, n.iter = 3000, cross = B6BTBR07, 
    ...) 
{
    require(qtlbim)
    pheno.col <- find.pheno(cross, mytrait(traitname, cross))
    cross <- subset(cross, chr = (sapply(cross$geno, class) == 
        "A"))
    cross$pheno$newtrait <- if (log10) 
        log10(cross$pheno[[pheno.col]])
    else cross$pheno[[pheno.col]]
    pheno.col <- find.pheno(cross, "newtrait")
    cross$pheno$sexcov <- unclass(cross$pheno[[grep("sex", tolower(names(cross$pheno)))]])
    not.miss <- !is.na(cross$pheno[[pheno.col]])
    if (any(not.miss)) 
        cross$pheno$sexcov <- (cross$pheno$sexcov - mean(cross$pheno$sexcov[not.miss]))/sqrt(var(cross$pheno$sexcov[not.miss]))
    assign("cross", cross, pos = 1)
    qb.mcmc(cross, pheno.col = pheno.col, n.iter = n.iter, interval = rep(30, 
        19), fixcov = find.pheno(cross, "sexcov"), intcov = TRUE, 
        ...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
