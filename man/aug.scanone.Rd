\name{aug.scanone}
\alias{aug.scanone}
\alias{plot.aug.scanone}
\alias{summary.aug.scanone}
\alias{plot.summary.aug.scanone}
\alias{plot.summary.aug.scanone}
\title{Plot 1-D scan of LOD and/or means}
\description{
Profiles of one or more phenotypes. If only one phenotype, in
addition profile the means by genotype.
}
\usage{
aug.scanone(traitnames = mytrait(),
                         cross = B6BTBR07,
                         sex = sexes,
                         method = "ehk",
                         log10 = rep(FALSE, length(traitnames)),
                         log.offset = 1,
                         lod.error = 100,
                         category = B6BTBR07.source,
                         ...)
\method{plot}{aug.scanone}(x, chr = levels(x$chr), traitnames = names(x)[-(1:2)],
  col.scheme = c("redblue", "cm", "gray", "heat", "terrain", "topo"),
  gamma = 0.6, allow.neg = FALSE, max.names = 50, zscale = TRUE,
  main = "", threshold.level = 0.05, max.lod = 20, category = NULL, \dots)
\method{summary}{aug.scanone}(object, chr = levels(object$chr),
  threshold.level = 0.05, mean.peaks = FALSE, category = NULL, \dots)
\method{print}{summary.aug.scanone}(x, digits = 2, \dots)
\method{plot}{summary.aug.scanone}(x, chr = dimnames(x$lod)[[2]],
  threshold.level = 0.05, max.lod = 20, max.names = 100,
  by = c("chr","trait","phenotype"), scale = c("cM","Mb"),
  cex = 2, pch = 3, \dots)
}
\details{
  \code{aug.scanone} creates multiple scanone's using
  \code{\link[qtl]{scanone}}. The plot uses ideas from
  \code{\link[qtl]{plot.scantwo}}. The \code{summary} method produces a
  large list, which can itself be plotted.
}
\seealso{\code{\link{myplot}}}
\examples{
multtrait.plot(cross.name="B6BTBR07",
  category="rbm", ## Later this will allow for tissues, modules.
  traitnames=mytrait(c("il.18","mpo")),
  chr=c(1:19,"X"),
  col.scheme=c("redblue", "cm", "gray", "heat", "terrain",
        "topo"),
  threshold.level=0.05,    ## Drop traits that have max below threshold.
  max.names=100,           ## Include names if number of traits < max.names.
  max.lod = 20)            ## Truncate lod at max.lod for color scheme.
}
\keyword{ ~kwd1 }
