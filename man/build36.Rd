\name{build36}
\alias{build36}
\docType{data}
\title{NCBI build 36 data for Affymetrix 5K Mouse Kit}
\description{
NCBI Mouse Build 36 data with best available position and source.
}
\usage{data(build36)}
\format{
  A data frame with 4853 observations on the following 14 variables.
  \describe{
    \item{\code{ID}}{Order of entry returned from NCBI.}
    \item{\code{snp}}{Build 36 \code{rs} number.}
    \item{\code{updated}}{Date of last update: \code{12/22/2005}
      or \code{4/19/2006}.}
    \item{\code{snpold}}{Build 32 \code{rs} number.}
    \item{\code{alleles}}{Target allele pair or triplet from \code{A/C}
      \code{A/C/G} \code{A/G} \code{A/G/T} \code{A/T} \code{C/G}
      \code{C/G/T} \code{C/T} \code{G/T}.}
    \item{\code{het}}{Heterozygosity score.}
    \item{\code{sehet}}{Standard error of \code{he}.}
    \item{\code{assemply}}{One of \code{C57BL/6J} \code{Celera}
      \code{build.32}. One SNP was not found in build 36 at all.}
    \item{\code{chr}}{Chromosome name: 1 to 19 or X}
    \item{\code{chr.pos}}{Chromosome position in Mb (build 36).}
    \item{\code{accession}}{Accession number (usually \code{NT_nnn}).}
      \item{\code{ctg.start}}{Base pair start of SNP.}
    \item{\code{ctg.end}}{Base pair end of SNP (usually same as \code{ctg.start}.}
    \item{\code{orient}}{Orientation as given \code{+} or reversed \code{-}.}
  }
}
\references{
\url{http://www.ncbi.nlm.nih.gov/projects/genome/seq/NCBIContigInfo.html}
}
\examples{
data(build36)
dim(build36)
}
\keyword{datasets}
