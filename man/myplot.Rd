\name{myplot}
\alias{myplot}
\alias{transcript.plot}
\title{Plot 1-D scan of LOD and/or means}
\description{
Profiles of one or more phenotypes. If only one phenotype, in
addition profile the means by genotype.
}
\usage{
transcript.plot(cross.name, ...)
myplot(traitnames, chr = c(1:19, "X"), add = TRUE,
 title, col = c("black", "blue", "red", "green", "purple", "magenta"),
 sex, log10 = FALSE, previous = FALSE,
 model = c("normal", "binary", "2part", "np"),
 method = "hk", means = c("add", "only", "none"),
 snpgrid = 30, threshold = 0.05, n.perm = 0, ylim,
 cross.name = "B6BTBR07", legend = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{traitnames}{Character string of trait names. Can be partial match.}
  \item{chr}{Chromosomes to profile.}
  \item{add}{Add all traits to one plot if \code{TRUE}.}
  \item{title}{Title for LOD profile plot.}
  \item{col}{Colors for LOD profile plot.}
  \item{sex}{Character string corresponding to \code{traitnames}. See details.}
  \item{log10}{Transform by \code{log10} if \code{TRUE}.}
  \item{previous}{Plot on previous plot if \code{TRUE}.}
  \item{model}{Model argument for \code{\link[qtl]{scanone}}.}
  \item{method}{Method argument for \code{\link[qtl]{scanone}}.}
  \item{newdata}{New data frame to add phenotypes to cross. See details.}
  \item{covariates}{Character string of covariate names. See details.}
  \item{means}{Add means profile if \code{"add"}; only show means
    profile if \code{"only"}; no means profile if \code{"none"}.}
  \item{means.by.sex}{Separate profiles of means by sex if \code{TRUE}
    and \code{sex="both"}.}
  \item{lods.by.sex}{Add male and female LOD profiles if
    \code{sex="both"} and only one trait.}
  \item{title.means}{Title for means profile plot(s).}
  \item{snpgrid}{Number of ticks on SNP Mb grid.}
  \item{threshold}{Threshold p-value for permutation test.}
  \item{n.perm}{Number of permutations; none if \code{n.perm} = 0.}
  \item{ylim}{Limits for vertical axis on LOD profile; default to range
    of values across traits.}
  \item{cross.name}{Name of cross.}
  \item{legend}{Include legend if \code{TRUE}.}
  \item{\dots}{Further arguments to \code{\link[qtl]{scanone}}.}
}
\details{
  \code{transcript.plot} is a version that is used explicitly for
  transcript gene mapping, or eQTL.
  
The \code{sex} option can be "both", "male", "female" or "ignore". With
"both", Sex is included as an interacting covariate. If in addition
there is only one trait, then plots for male and female are overlaid if
\code{lods.by.sex} is \code{TRUE}.

New phenotypes can be added by passing a data frame with
\code{newdata}. This data frame must have \code{MouseNum} as the first
column, with entries \code{Mousemmmm} as IDs. Other columns must be
labeled with valid R names (no special characters). Those names may be
used as phenotypes or covariates. Covariates except Sex are used as
additive covariates only.
}
\seealso{B6BTBR07}
\examples{
myplot("il.18",11)
myplot(c("il.18","mpo"))
}
\keyword{ ~kwd1 }
