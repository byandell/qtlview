#######################################################
myrecover <- function(traitname, log10 = FALSE, n.iter = 3000,
                      cross, output.dir, ...)
{
  require(qtlbim)

  pheno.col <- find.pheno(cross, mytrait(traitname, cross$pheno))

  ## Drop X chr.
  cross <- subset(cross, chr = (sapply(cross$geno, class) == "A"))
  
  ## Create new trait, transform if desired.
  cross$pheno$newtrait <- if(log10)
    log10(cross$pheno[[pheno.col]])
  else
    cross$pheno[[pheno.col]]
  pheno.col <- find.pheno(cross, "newtrait")

  ## Set up sex covariate centered on 0 w.r.t. missing values for trait.
  cross$pheno$sexcov <-
    unclass(cross$pheno[[grep("sex", tolower(names(cross$pheno)))]])
  not.miss <- !is.na(cross$pheno[[pheno.col]])
  if(any(not.miss))
    cross$pheno$sexcov <-
      (cross$pheno$sexcov -
       mean(cross$pheno$sexcov[not.miss])) /
         sqrt(var(cross$pheno$sexcov[not.miss]))

  ## qb.mcmc has problems--qb.valid.phenoData.
  ## check updates I missed.  Do diff.
  assign("cross", cross, pos = 1)

  ## Recover qb object and make sure args set properly.
  qb <- qb.recover(cross, "newtrait", output.dir = output.dir,
                   interval = rep(30, 19),
                   fixcov = find.pheno(cross, "sexcov"),
                   intcov = TRUE, ...)
  qb$args$interval <- rep(30, 19)
  qb$args$fixcov <- find.pheno(cross, "sexcov")
  qb$args$intcov <- TRUE
  qb
}
#######################################################
mylegacy <- function(qbObject,
                     traitname, log10 = FALSE, n.iter = 3000,
                     cross, ...)
{
  require(qtlbim)

  if(missing(traitname))
    stop("need traitname")
  
  pheno.col <- find.pheno(cross, mytrait(traitname, cross$pheno))

  ## Drop X chr.
  cross <- subset(cross, chr = (sapply(cross$geno, class) == "A"))

  ## Create new trait, transform if desired.
  cross$pheno$newtrait <- if(log10)
    log10(cross$pheno[[pheno.col]])
  else
    cross$pheno[[pheno.col]]
  pheno.col <- find.pheno(cross, "newtrait")

  ## Set up sex covariate centered on 0 w.r.t. missing values for trait.
  cross$pheno$sexcov <-
    unclass(cross$pheno[[grep("sex", tolower(names(cross$pheno)))]])
  not.miss <- !is.na(cross$pheno[[pheno.col]])
  if(any(not.miss))
    cross$pheno$sexcov <-
      (cross$pheno$sexcov -
       mean(cross$pheno$sexcov[not.miss])) /
         sqrt(var(cross$pheno$sexcov[not.miss]))

  ## qb.mcmc has problems--qb.valid.phenoData.
  ## check updates I missed.  Do diff.
  assign("cross", cross, pos = 1)
  
  qb <- qb.legacy(qbObject)
  qb
}
