\name{plot.snp.record}
\alias{plot.snp.record}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.snp.record(x = snp.record[[diagname]], diagname = deparse(substitute(x)),
  chr = c(1:19, "X"), title, map ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{diagname}{ ~~Describe \code{diagname} here~~ }
  \item{chr}{ ~~Describe \code{chr} here~~ }
  \item{title}{ ~~Describe \code{title} here~~ }
  \item{map}{ ~~Describe \code{title} here~~ }
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
function (x = snp.record[[diagname]], diagname = deparse(substitute(x)), 
    chr = c(1:19, "X"), title = {
        if (length(chr) == 1) 
            paste("chr ", chr, ":", diagname)
        else diagname
    }, ...) 
{
    if (is.character(x) & length(x) == 1) {
        diagname <- x
        x <- snp.record[[diagname]]
    }
    if (is.null(x)) 
        stop("invalid x or diagname")
    require(qtl)
    map <- pull.map(B6BTBR07)
    markers <- unlist(lapply(map, names))
    scanone.obj <- data.frame(chr = ordered(rep(names(map), sapply(map, 
        length)), names(map)), pos = unlist(map))
    tmp <- rep(NA, nrow(scanone.obj))
    tmp2 <- match(snp.record$Build36.SNP, markers, nomatch = 0)
    tmp[tmp2] <- x[tmp2 > 0]
    scanone.obj[[diagname]] <- tmp
    class(scanone.obj) <- c("scanone", "data.frame")
    plot(scanone.obj, chr = chr, ...)
    if (length(chr) == 1) {
        rug(cM.map[[chr]], 0.02, side = 3)
        p <- myapprox(Mb.map[[chr]], cM.map[[chr]])
        axis(3, p$y, p$x)
        mtext("Mb", 3, 1, at = -5, adj = 1)
    }
    title(title, line = 2)
    invisible()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
