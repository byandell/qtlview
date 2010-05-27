\name{read.maps}
\alias{read.maps}
\title{Create genetic and physical list of maps object.}
\description{
Use cross and file of physical locations to create a list of R/qtl map objects.
}
\usage{
read.maps(cross, filename, chr.valid, keep.nonseg = TRUE,
  drop.extra = TRUE, reset.chr = TRUE, interp.loc = TRUE,
  genotypes = c("A","H","B"), verbose = TRUE, \dots)
\method{plot}{read.maps}(x, main = "genetic vs. physical maps", \dots)
\method{summary}{read.maps}(object, \dots)
}
\arguments{
  \item{cross}{Objects of class \code{cross}. See
    \code{\link[qlt]{read.cross}}.}
  \item{filename}{Name of file containing physical locations and
    detailed descriptions of markers.}
  \item{chr.valid}{Names of chromosomes to keep (default is all in
    \code{cross} object).}
  \item{keep.nonseg}{Keep non-segregating markers if \code{TRUE}
    (default). May be vector of marker names, or logical vector with
    names corresponding to markers.}
  \item{drop.extra}{Drop extra markers in \code{cross} object if not
    found in file when \code{TRUE} (default). If \code{FALSE}, then
    physical locations for extra markers will be interpolated.}
  \item{reset.chr}{Reset \code{chr} value for markers with missing
    values in physical location file using \code{cross} object if
    \code{TRUE} (default). If \code{FALSE}, markers 
    with missing \code{chr} values are dropped.}
  \item{interp.loc}{Interpolate physical location if missing using
    \code{cross} object when \code{TRUE} (default).}
  \item{genotypes}{Genotype names to be used in plots and summaries.}
  \item{verbose}{Verbose listing if \code{TRUE} (default).}
  \item{x,object}{Object of class \code{read.maps}.}
  \item{main}{Main title for plot.}
  \item{\dots}{Arguments passed to \code{plot} or \code{summary} method.}
}
\details{
  Supplied by NCBI via Ron Zhu.
  The object created by this routine has both physical (Mb) and
  genetic (cM) distances along each chromosome. It can have markers that
  are part of map, as well as markers that are the same for parents
  (suggesting non-segregating regions of nearly same genotype).

  Routine assumes physical location file has columns identified as
  \code{snp} = marker ID;
  \code{chr} = chromosome name (must match \code{cross} object);
  \code{loc} = physical location in Morgans;
  \code{orient} = orientation on chromosome (ignored for now).
}
\seealso{
\code{\link{snp.record}}
}
\examples{
\dontrun{
## Keep only SNPs with decent records.
keep <- snp.record$Keep
names(keep) <- snp.record$Build36.SNP
## Create new maps object.
B6BTBR07.maps <-
  read.maps(B6BTBR07, "http://attie.wisc.edu/lab/data/data_for_intranet/snps/snp_4853.txt", keep = keep)
}
}
\keyword{utilities}
