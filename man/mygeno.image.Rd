\name{mygeno.image}
\alias{mygeno.image}
\title{image of genotypes}
\description{
  image of genotypes
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
