\name{mytrait}
\alias{mytrait}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mytrait(traitnames, cross = B6BTBR07)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traitnames}{ ~~Describe \code{traitnames} here~~ }
  \item{cross}{ ~~Describe \code{cross} here~~ }
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
function (traitnames, cross = B6BTBR07) 
{
    if (missing(traitnames)) 
        return(names(cross$pheno))
    tmp <- names(cross$pheno)[pmatch(tolower(make.names(traitnames)), 
        tolower(names(cross$pheno)), duplicates.ok = TRUE)]
    tmp1 <- is.na(tmp)
    if (any(tmp1)) {
        tmp3 <- sapply(strsplit(traitnames[tmp1], "wk"), paste, 
            collapse = "")
        tmp2 <- pmatch(tolower(make.names(paste("wk", tmp3, sep = ""))), 
            tolower(names(cross$pheno)), nomatch = 0)
        tmp2 <- pmax(pmatch(tolower(make.names(paste("X", tmp3, 
            sep = ""))), tolower(names(cross$pheno)), nomatch = 0), 
            tmp2)
        tmp2[tmp2 == 0] <- NA
        tmp2 <- names(cross$pheno)[tmp2]
        tmp[tmp1] <- tmp2
        if (any(is.na(tmp))) {
            print(traitnames)
            print(tmp)
            stop("traitnames do not match names(cross$pheno)")
        }
    }
    tmp
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
