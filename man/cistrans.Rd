\name{cistrans}
\alias{cistrans}
\alias{plot.cistrans}
\alias{print.cistrans}
\alias{summary.cistrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
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
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{filename}{ ~~Describe \code{filename} here~~ }
  \item{peak.chr}{ ~~Describe \code{peak.chr} here~~ }
  \item{trans.chr}{ ~~Describe \code{trans.chr} here~~ }
  \item{min.lod}{ ~~Describe \code{min.lod} here~~ }
  \item{main}{ ~~Describe \code{main} here~~ }
  \item{sep}{ ~~Describe \code{sep} here~~ }
  \item{use.density}{ ~~Describe \code{use.density} here~~ }
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
function (x = read.table(filename, header = TRUE, sep = sep), 
    filename, peak.chr = peak.chrs, trans.chr = trans.chrs, min.lod = 0, 
    main = paste(levels(x$Tissue), collapse = ", "), sep = "\t", 
    use.density = TRUE, ...) 
{
    getreal <- function(x, realx = c(1:19, "X")) {
        if (is.null(x)) 
            NULL
        else {
            x <- levels(as.ordered(x))
            realx[realx \%in\% x]
        }
    }
    peak.chrs <- getreal(x$peak.chr)
    trans.chrs <- getreal(x$trans.chr)
    maxlods <- c("maxlod", "maxlod.chr", "maxlod.pos.Mb", "maxlod.pos.cM")
    peaks <- c("peak.score", "peak.chr", "peak.pos.Mb", "peak.pos.cM")
    tmp <- match(maxlods, names(x), nomatch = 0)
    names(x)[tmp] <- peaks[tmp > 0]
    x$peak.chr <- ordered(x$peak.chr, peak.chrs)
    is.trans <- !is.null(x$trans.chr)
    if (is.trans) 
        is.trans <- !all(is.na(x$trans.chr))
    if (is.trans) 
        x$trans.chr <- ordered(x$trans.chr, trans.chrs)
    tmp <- names(cM.map)
    peak.chr <- tmp[tmp \%in\% peak.chr]
    if (is.trans) {
        trans.chr <- tmp[tmp \%in\% trans.chr]
        if (!length(peak.chr) | !length(trans.chr)) 
            stop("must select at least on peak chr and one trans chr")
    }
    else {
        if (!length(peak.chr)) 
            stop("must select at least on peak chr")
    }
    tmp <- (x$peak.chr \%in\% peak.chr)
    if (is.trans) 
        tmp <- tmp & (x$trans.chr \%in\% trans.chr)
    if (!any(tmp)) 
        stop("no transcripts selected")
    x <- x[tmp, ]
    tmp <- x$peak.score >= min.lod
    if (!any(tmp)) 
        stop(paste("no transcripts above minimum score of", min.lod))
    x <- x[tmp, ]
    drop.cis <- FALSE
    if (drop.cis & is.trans) {
        dist.cis <- 20
        tmp <- as.character(x$peak.chr) == as.character(x$trans.chr)
        tmp <- tmp & (abs(x$peak.pos.Mb - x$trans.pos.Mb) <= 
            dist.cis)
        if (any(tmp)) 
            x <- x[!tmp, ]
    }
    if (is.null(x$peak.pos.cM)) 
        x$peak.pos.cM <- mypos(x, "Mb", "peak")
    if (is.null(x$trans.pos.cM) & is.trans) 
        x$trans.pos.cM <- mypos(x, "Mb", "trans")
    x <- cumscore(x, ...)
    class(x) <- c("cistrans", "data.frame")
    if (use.density) 
        attr(x, "dens") <- cumscore.dens(x, ...)
    attr(x, "main") <- main
    attr(x, "peak.chr") <- peak.chr
    attr(x, "trans.chr") <- trans.chr
    attr(x, "is.trans") <- is.trans
    x
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
