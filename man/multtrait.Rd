\name{multtrait}
\alias{multtrait}
\alias{plot.multtrait}
\alias{print.multtrait}
\alias{summary.multtrait}
\title{Plot 1-D scan of LOD and/or means}
\description{
  Profiles of one or more phenotypes using raw data or pre-computed
  profiles.
  If only one phenotype, in addition profile the means by genotype.
}
\usage{
multtrait(traitnames = NULL, cross.name, cross = get(cross.name),
  category, filename = NULL, chr="", main = "",
  summary.plot = c("total", "chr","none","phenotype"),
  threshold.level=0.05, threshold.lod,
  cluster = TRUE, support.lod = 1.5,
  n.clust,
  heatmap = FALSE,
  scan.type = c("LOD","LPD","BF","RAW","MOM","PAT"),
  sex = c("both","male","female"),
  ylab = c("symbol","a_gene_id","symbol.a_gene_id","none","on"),
  \dots)  
\method{plot}{multtrait}(x, heatmap = TRUE, rescale = TRUE, use.cM = FALSE,
  summary.plot, \dots)
\method{summary}{multtrait}(object, \dots)
\method{print}{multtrait}(x, \dots)
}
\details{
  \code{multtrait} is intended as the Intranet interface.
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
