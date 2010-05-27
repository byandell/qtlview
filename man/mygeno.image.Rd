\name{mygeno.image}
\alias{mygeno.image}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mygeno.image(x, chr, reorder = TRUE, id = "MouseNum",
  sex = c("both", "male", "female"), genotypes, xlim,
  use.cM = FALSE, reorder.by.genotype, recomb.only = FALSE,
  keep.missing = TRUE, \dots)
\method{plot}{mygeno.image}(x, xlim, use.cM = FALSE, equal.spacing = FALSE,
  zscale = TRUE, normal.score = TRUE, main = "", \dots)
\method{summary}{mygeno.image}(object, \dots)
\method{print}{mygeno.image}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class cross (see \code{\link[qtl]{read.cross}})}
  \item{chr}{vector of chromosome names for object \code{x}}
  \item{reorder}{character vector of traits to align with genotype}
  \item{id}{name used to match individual identifiers}
  \item{sex}{choice of both, male or female sex}
  \item{genotypes}{character vector of genotype names for summary (see below)}
  \item{xlim}{range on chr to consider (ignored if length(chr) > 1)}
  \item{use.cM}{use cM map if \code{TRUE}; otherwise use Mb map}
  \item{reorder.by.genotype}{reorder by average genotype if \code{TRUE} (see below)}
  \item{recomb.only}{only show recombinant individuals if \code{TRUE}}
  \item{equal.spacing}{use equal spacing for markers if \code{TRUE}}
  \item{keep.missing}{keep individuals with missing genotypes if
    \code{TRUE}}
  \item{zscale}{include Z gray scale if \code{TRUE}}
  \item{normal.score}{use normal scores of phenotypes if \code{TRUE}}
  \item{main}{main label for plot}
  \item{\dots}{additional arguments to \code{\link[qtl]{subset.cross}}}
}
\details{
  Produces an image map of genotypes. Individuals are sorted by average
  genotype (e.g. number of B alleles) if \code{reorder} is \code{TRUE}.
  Genotypes are color-coded as blue (AA), green (AB), red (BB). If
  non-segregating markers are included, regions spanned by them are
  colored lightgray. Gray scale on right is mean genotype across markers
  (black = AA, white = BB).

  If names of phenotypes are provided with \code{reorder}, these are
  visually shown next to genotype bar on right as smooth and raw
  forms. In addition, a second plot of genotype vs. phenotype shows the
  relationship in detail.
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
function (x, chr, reorder = TRUE, id = "MouseNum", sex = c("both", 
    "male", "female"), genotypes = genos, xlim = range(Map), 
    use.cM = FALSE, reorder.by.genotype = any(reorder > 0), recomb.only = FALSE, 
    keep.missing = TRUE, ...) 
{
    sex <- match.arg(sex)
    if (!inherits(x, "cross")) 
        stop("Input should have class \"cross\".")
    cross <- subset(x, chr = chr, ...)
    if (sex != "both") {
        sexpgm <- getsex(cross)
        cross <- subset(cross, ind = (sexpgm$sex == (sex == "male")))
    }
    chr <- names(cross$geno)
    ids <- cross$pheno[[id]]
    Geno <- pull.geno(cross)
    cM.map <- pull.map(cross)
    if (length(chr) == 1) {
        Map <- if (use.cM) 
            cM.map[[chr]]
        else Mb.map[[chr]]
        xlim <- xlim.range(xlim, Map)
        mean.geno <- apply(Geno[, (Map >= xlim[1] & Map <= xlim[2])], 
            1, mean, na.rm = TRUE)
    }
    else mean.geno <- apply(Geno, 1, mean, na.rm = TRUE)
    ordset <- NULL
    if (is.character(reorder)) 
        reorder <- find.pheno(cross, reorder)
    if (is.numeric(reorder)) {
        if (reorder < 1 || reorder > nphe(cross)) 
            stop("reorder should be an integer between 1 and", 
                nphe(cross))
        ordset <- as.matrix(cross$pheno[, reorder, drop = FALSE])
    }
    else ordset <- NULL
    if (reorder.by.genotype) {
        o <- order(mean.geno)
        Geno <- Geno[o, ]
        ids <- ids[o]
        mean.geno <- mean.geno[o]
        if (!is.null(ordset)) 
            ordset <- ordset[o, , drop = FALSE]
    }
    if (recomb.only) {
        keep <- apply(Geno[, (Map >= xlim[1] & Map <= xlim[2])], 
            1, function(x) !all(x == mean(x, na.rm = TRUE)))
        keep[is.na(keep)] <- FALSE
        Geno <- Geno[keep, , drop = FALSE]
        ids <- ids[keep]
        mean.geno <- mean.geno[keep]
        if (!is.null(ordset)) 
            ordset <- ordset[keep, , drop = FALSE]
    }
    if (!keep.missing) {
        keep <- apply(Geno[, (Map >= xlim[1] & Map <= xlim[2])], 
            1, function(x) !any(is.na(x)))
        Geno <- Geno[keep, , drop = FALSE]
        ids <- ids[keep]
        mean.geno <- mean.geno[keep]
        if (!is.null(ordset)) 
            ordset <- ordset[keep, , drop = FALSE]
    }
    class(Geno) <- c("mygeno.image", "matrix")
    attr(Geno, "cM") <- cM.map
    attr(Geno, "Mb") <- Mb.map[chr, drop = FALSE]
    attr(Geno, "chr") <- chr
    attr(Geno, "id") <- ids
    attr(Geno, "geno") <- mean.geno
    attr(Geno, "ordset") <- ordset
    attr(Geno, "xlim") <- xlim
    attr(Geno, "recomb.only") <- recomb.only
    attr(Geno, "keep.missing") <- keep.missing
    genos <- c("B6", "Het", "BTBR", "same")
    tmp <- length(genotypes)
    if (tmp >= 4) 
        genotypes <- genotypes[1:4]
    else genotypes <- c(genotypes, genos[(tmp + 1):4])
    attr(Geno, "genotypes") <- genotypes
    col <- c("blue", "green", "red", "lightgray")
    names(col) <- c(1:3, "NA")
    attr(Geno, "col") <- col
    Geno
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
