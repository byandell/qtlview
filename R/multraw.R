###########################################################  
## Routine multtrait assumes cross, perms and maps are available.
## Default names are cross.name, cross.name.perms, cross.name.maps.

multraw <- function(traitnames = NULL,
                    cross = get(cross.name),
                    cross.name = deparse(substitute(cross)),
                    chr = "",
                    main = "",
                    sex = sexes[1],
                    method = "hk",
                    filename = NULL,
                    newdata = NULL,
                    covariates = NULL,
                    threshold.level = 0.05,
                    lods.by.sex = TRUE,
                    perms = myget(cross.name, "perms"),
                    maps = myget(cross.name, "maps"),
                    density = FALSE,
                    category = "clinical",
                    scan.type = c("LOD","LPD","BF"),
                    ylab = c("symbol","a_gene_id","symbol.a_gene_id","none"),
                    step = 2, off.end = 0, error.prob = 0.0001, 
                    map.function = "haldane",
                    stepwidth = "variable",
                    trait.annotation = NULL,
                    ...)
{
  scan.type <- match.arg(scan.type)
  ylab <- match.arg(ylab)

  newdata <- get.raws(filename, newdata)

  if(!length(traitnames)) {
    if(is.null(newdata)) {
      if(!length(traitnames))
        stop("need to supply trait names to profile raw data")
      traitnames <- mytrait(traitnames, cross$pheno)
    }
    else
      traitnames <- names(newdata)[-1]
  }

  chr <- find.chr(chr, names(cross$geno))

  cross <- add.raws(cross, newdata, traitnames, covariates)

  n.profiles <- n.traits <- length(traitnames)
  
  sexes <- c("both","male","female","ignore")
  sex <- sexes[pmatch(sex, sexes, nomatch = 1)]
  sex <- array(sex, n.traits)

  sexpgm <- getsex(cross)

  ## Transform phenotypes.
  trans.pheno <- cross$pheno[,traitnames, drop = FALSE]

  ## Make sure we have genoprob; compute if missing.
  if (!("prob" %in% names(cross$geno[[1]])))
    cross <- calc.genoprob(cross, step, off.end, error.prob, map.function, stepwidth)

  for(i in seq(n.traits)) {
    if(is.logical(trans.pheno[[i]]))
      trans.pheno[[i]] <- as.numeric(trans.pheno[[i]])
    if(is.factor(trans.pheno[[i]]))
      trans.pheno[[i]] <- as.numeric(trans.pheno[[i]]) - 1
    trans.pheno[[i]] <- normal.trans(trans.pheno[[i]])

    ## Make sure sex is set properly.
    sextbl <- table(sexpgm$sex[!is.na(trans.pheno[[i]])])
    if(!sum(sextbl))
      stop(paste("all missing values for trait", traitnames[i]))
    if(sex[i] == "both") {
      if("0" %in% names(sextbl)) {
        ## Only females.
        if(!("1" %in% names(sextbl))) {
          sex[i] <- "female"
          warning(paste("only females for trait", traitnames[i]))
        }
      }
      else {
        ## Only males.
        sex[i] <- "male"
        warning(paste("only males for trait", traitnames[i]))
      }
    } else {
      if(sex[i] == "male" & !("1" %in% names(sextbl)))
        stop(paste("no males for trait", traitnames[i]))
      if(sex[i] == "female" & !("0" %in% names(sextbl)))
        stop(paste("no females for trait", traitnames[i]))
    }
  }

  lods.by.sex <- lods.by.sex & (n.traits == 1) & (sex[1] == "both")
  if(lods.by.sex & !density) {
    n.profiles <- 3
    traitnames <- rep(traitnames, 3)
    sex <- c("both","male","female")
  }

  sex <- array(sex, n.profiles)

  ## Need to modify this to allow for X chromosome threshold.
  threshold.lod <- rep(0, n.profiles)
  threshold.level <- threshold.level[1]
  if(threshold.level == 1)
    threshold.level <- 0

  ## Shut off warning.
  tmpopt <- options(warn = -2)
  on.exit(options(tmpopt))

  sums <- list()
  for(i in seq(traitnames)) {
    traitname <- traitnames[i]

    ## Set up cross for sex.
    switch(sex[i],
           both =, ignore = {
             crosssex <- cross
             crosssex$pheno[[traitname]] <- trans.pheno[[traitname]]
           },
           male = {
             crosssex <- subset(cross, ind = (sexpgm$sex == 1))
             crosssex$pheno[[traitname]] <-
               trans.pheno[sexpgm$sex == 1, traitname]
           },
           female = {
             crosssex <- subset(cross, ind = (sexpgm$sex == 0))
             crosssex$pheno[[traitname]] <-
               trans.pheno[sexpgm$sex == 0, traitname]
           })

    ## Covariates (should be handled outside this routine).
    covlist <- mycov(sexpgm, sex[i], category, covariates, crosssex)
  
    pheno.col <- find.pheno(crosssex, traitname)

    ## Profile LOD.
    sums[[i]] <- scanone(crosssex, chr, pheno.col, method = method,
                         intcov = covlist$intcov, addcov = covlist$addcov)

    ## Need to modify this to allow for X chromosome threshold.
    threshold.lod[i] <- threshold.perm(perms, sex[i], threshold.level)
  }

  mynames <- paste(traitnames, sex, sep = ".")
  names(threshold.lod) <- names(sums) <- names(sex) <- mynames

  ## Add trait positions.
  trait.position <-
    find.trait.position(traitnames, maps, trait.annotation, ...)

  out <- list(threshold.level = threshold.level,
              threshold.lod = threshold.lod, sex = sex,
              sums = sums, traitnames = traitnames,
              n.traits = n.traits, trans.pheno = trans.pheno,
              tissue.batch = covlist$tissue.batch, cross = cross,
              main = main, ylab = ylab, maps = maps,
              trait.position = trait.position)

  class(out) <- c("multraw", "list")
  out
}
###########################################################
## This should be removed ultimately.
mycov <- function(sexpgm, sex, category, covariates, crosssex)
{
  ## Interacting covariate: sex.
  intcov <- if(sex == "both")
    as.matrix(sexpgm$sex)
  else
    NULL
  
  ## Additive covariate
  addcov <- intcov
  
  ## Add batch if category not clinical.
  if(category != "clinical") {
    tissue.batch <- paste(category, "batch", sep = ".")
    form <- formula(paste("~", tissue.batch))
    batchcov <- model.matrix(form, crosssex$pheno)[, -1]
    tmp <- matrix(NA, nind(crosssex), ncol(batchcov))
    tmp[!is.na(crosssex$pheno[[tissue.batch]]), ] <- batchcov
    addcov <- cbind(addcov, tmp)
  }
  else
    tissue.batch <- NULL
  
  ## Add covariates if any.
  if(!is.null(covariates)) {
    covariates <- mytrait(covariates, crosssex$pheno)
    addcov <- cbind(addcov, crosssex$pheno[, covariates])
  }
  list(intcov = intcov, addcov = addcov, tissue.batch = tissue.batch)
}
###########################################################  
print.multraw <- function(x, ...) print(summary(x, ...), ...)
###########################################################  
summary.multraw <- function(object, traitnames = object$traitnames,
                            chr = levels(object$sums[[1]]$chr),
                            abbrev = 20,
                            ...)
{
  if(!all(traitnames %in% object$traitnames))
    stop(paste("traitnames match with object:",
               paste(object$traitnames, collapse = ",")))

  sums <- lapply(object$sums, sumone.scanone, chr = chr)
  
  sumall <- sums[[1]][, 1, drop = FALSE]
  if(length(unique(object$traitnames)) == 1) {
    for(i in seq(length(object$sex))) {
      names(sums[[i]]) <-
        paste(names(sums[[i]]), object$sex[i], sep = ".")
      sumall <- cbind(sumall, sums[[i]][, -1])
    }
  }
  else {
    ## Abbreviate traitnames. First remove annoying "."s.
    trnames <- traitnames
    trnames <- sapply(strsplit(trnames, ".", fixed = TRUE),
                      function(x) paste(x, collapse = ""))
    if(abbrev > 0)
      trnames <- abbreviate(trnames, abbrev)
    names(trnames) <- traitnames
  
    for(i in match(traitnames, object$traitnames)) {
      names(sums[[i]]) <-
        paste(names(sums[[i]]), trnames[object$traitnames[i]], sep = ".")
      sumall <- cbind(sumall, sums[[i]][, -1])
    }
  }
  class(sumall) <- c("summary.multraw", class(sums[[1]]))
  attr(sumall, "threshold.level") <- object$threshold.level
  attr(sumall, "threshold.lod") <- object$threshold.lod
  attr(sumall, "sex") <- object$sex
  attr(sumall, "trait.position") <- object$trait.position

  sumall
}
###########################################################  
print.summary.multraw <- function(x, ...)
{
  threshold.level <- attr(x, "threshold.level")
  if (threshold.level > 0) {
    cat(paste("Permutation thresholds at ", threshold.level,
              "% level\n", sep = ""))

    ## Threshold by sex.
    sex <- attr(x, "sex")
    tmp <- !duplicated(sex)
    cat("Threshold values:",
        paste(sex[tmp],
              round(attr(x, "threshold.lod")[tmp], 2),
              sep = " = ", collapse = ", "), "\n\n")
  }

  cat("Trait positions (if known):\n")
  print(attr(x, "trait.position"))
  cat("\n")

  NextMethod(x, ...)
}
###########################################################  
plot.multraw <- function(x, chr = "",
                         add = TRUE, 
                         title = mytitle,
                         col = c("black","blue","red","green","purple","magenta"),
                         previous = FALSE,
                         covariates = NULL,
                         means = c("add","only","none"),
                         means.by.sex = (length(chr) > 1),
                         title.means = mytitle.means,
                         ylim = ylims,
                         legend = TRUE,
                         density = FALSE,
                         use.cM = FALSE,
                         ...)
{
  chr <- find.chr(chr, levels(x$sums[[1]]$chr))
                  
  n.profiles <- length(x$traitnames)
  
  means <- match.arg(means)
  show.means <- (means != "none") & (x$n.traits == 1)
  legend <- legend & add

  col <- array(col, n.profiles)

  means.by.sex <- means.by.sex & show.means & (x$sex[1] == "both")

  ## Use traitnames or no labels.
  tmp <- mylabels(x$traitnames, ylab = x$ylab)
  if(is.logical(tmp)) {
    mynames <- ""
    legend <- FALSE
  }
  else
    mynames <- paste(tmp, "(", x$sex, "/", col, ")", sep = "")

  ## Default title for LOD profile.
  if(x$main == "")
    mytitle <- paste(mynames, collapse = ", ")
  else
    mytitle <- x$main

  if(density & x$n.traits) { ## Do density instead.
    require(lattice)

    ## Limit to 5 traits.
    traitnames <- x$traitnames[seq(min(x$n.traits, 5))]

    ## Extract pheno from cross.
    pheno <- x$cross$pheno
    
    ## Make traitnames legit names.
    tmp <- match(traitnames, names(pheno))
    traitnames <- make.names(traitnames)
    names(pheno)[tmp] <- traitnames

    ## Make sex as Female, Male.
    sexpgm <- getsex(x$cross)
    pheno$sex <- factor(c("Female","Male")[1 + sexpgm$sex])

    ## Plot densities.
    form <- formula(paste("~", paste(traitnames, collapse = "+"), "| sex"))
    plot(densityplot(form, pheno, ref = TRUE, outer = TRUE))
    return(invisible())
  }

  ## Default title and colors for means profiles.
  sex.col <- c("blue", "green", "red")
  sex.names <-  paste(c("B6", "Het", "BTBR"), sex.col, sep = "=")
  mytitle.means <- paste(sex.names, collapse = ", ")
  if(any(x$sex == "both") & !means.by.sex)
    mytitle.means <- paste(mytitle.means, "female=solid, male=dashed",
                           sep = "; ")

  old.mar <- par("mar")
  old.las <- par("las")
  old.mfrow <- par("mfrow")
  on.exit(par(las = old.las, mar = old.mar, mfrow = old.mfrow))
  new.mar <- c(3.1,4.1,1.6 + 2 * (length(chr) == 1),0.6)
  par(mar = new.mar)

  n.plots <- 1 + ((!add) * (n.profiles - 1)) + ((means == "add") & show.means) +
    means.by.sex
  if(legend) {
    layout(cbind(seq(n.plots),n.plots+seq(n.plots)), c(6,1))
  }
  else
    par(mfrow = c(n.plots, 1))

  ## Set up cM to Mb translation.  
  if(length(chr) == 1) {
    base <- ifelse(use.cM, "cM", "Mb")
    p <- qm.approx(x$maps, base, chr)
    if(!use.cM) { ## use Mb.
      for(i in seq(x$traitnames))
        x$sums[[i]]$pos <- qm.approx(maps, "cM", chr, x$sums[[i]]$pos)$y
    }
  }


  ## Set vertical limits to include all peaks.
  ylims <- range(0, sapply(x$sums, function(x) max(x)$lod))

  for(i in seq(x$traitnames)) {
    traitname <- x$traitnames[i]
  
    plot(x$sums[[i]], add = previous | (add & i > 1),
         col = col[i], ylim = ylim)

    threshold.lines(x$sums[[i]], x$threshold.lod[i], col = col[i], ...)

    if(i == 1 | !add)
      add.rug(chr, title, x$maps, p, use.cM)
  }

  one <- x$sums[[1]]

  ## This is very specific to mouse crosses.
  if(show.means) {
    phenotypes <- x$cross$pheno

    ## Transform trait using normal scores.
    phenotypes$trans.pheno <- x$trans.pheno[, 1]

    ## Remove batch effect if category is tissue.
    if(!is.null(x$tissue.batch)) {
      form <- formula(paste("trans.pheno",
                            x$tissue.batch, sep = "~"))
      tmp <- !is.na(phenotypes$trans.pheno)
      phenotypes$trans.pheno[tmp] <- resid(lm(form, phenotypes))
    }
    
    make.mean <- function(cross, chr, phenotype) {
      lapply(cross$geno[chr],
             function(w, x) apply(w$prob, 2:3, function(w, x)
                                  weighted.mean(x,w, na.rm=TRUE),
                                  phenotype))
    }
    if(x$sex[1] == "ignore") {
      tmp <- make.mean(x$cross, chr, phenotype)
      one$AA <- unlist(lapply(tmp, function(x) x[,1]))
      one$AB <- unlist(lapply(tmp, function(x) x[,2]))
      one$BB <- unlist(lapply(tmp, function(x)
                              if(ncol(x) == 3) x[,3] else rep(NA, nrow(x))))
    }
    if(x$sex[1] %in% c("female","both")) {
      tmp <- make.mean(subset(x$cross, ind = (getsex(x$cross)$sex == 0)),
                       chr, phenotype)
      one$AA.F <- unlist(lapply(tmp, function(x) x[,1]))
      one$AB.F <- unlist(lapply(tmp, function(x) x[,2]))
      one$BB.F <- unlist(lapply(tmp, function(x)
                                if(ncol(x) == 3) x[,3] else rep(NA, nrow(x))))
    }
    if(x$sex[1] %in% c("male","both")) {
      tmp <- make.mean(subset(x$cross, ind = (getsex(x$cross)$sex == 1)),
                       chr, phenotype)
      one$AA.M <- unlist(lapply(tmp, function(x) x[,1]))
      one$AB.M <- unlist(lapply(tmp, function(x) x[,2]))
      one$BB.M <- unlist(lapply(tmp, function(x)
                                if(ncol(x) == 3) x[,3] else rep(NA, nrow(x))))
    }
    switch(x$sex[1],
           both = {
             if(means.by.sex) {
               names(one)[4] <- "female"
               names(one)[7] <- "male"
             }
             else
               names(one)[4] <- "female=solid, male=dashed"
           },
           ignore = {
             names(one)[4] <- "means"
           },
           male =, female = {
             names(one)[4] <- x$sex[1]
           })
    
    plot(one, lodcol = 2:4, col = sex.col,
         ylim = range(one[, 4:ncol(one)], na.rm = TRUE))

    add.rug(chr, title.means, x$maps, p, use.cM)

    if (ncol(one) > 6) {
      if (means.by.sex) {
        plot(one, lodcol = 5:7, col = sex.col,
             ylim = range(one[, 4:ncol(one)], na.rm = TRUE))
        add.rug(chr, title.means, x$maps, p, use.cM)
      }
      else
        plot(one, lodcol = 5:7, col = sex.col,
             lty = 2, add = TRUE)
    }
  }

  ## Add legend if desired.
  if(legend) {
    new.mar[c(2,4)] <- 0
    par(mar = new.mar)
    old.xaxt <- par("xaxt")
    old.yaxt <- par("yaxt")
    old.bty <- par("bty")
    
    on.exit(par(xaxt = old.xaxt, yaxt = old.yaxt, bty = old.bty), add = TRUE)
    par(xaxt = "n", yaxt = "n")
    plot(c(1,2), c(0,n.profiles), type = "n", xlab = "", ylab = "")
    for(i in seq(n.profiles)) {
      lines(c(1,2), rep(i-1,2), col = col[i], lwd = 2)
      text(1, i-0.5, mynames[i], col = col[i], adj = 0)
    }
    if(add & n.plots > 1) {
      for(j in seq(n.plots - 1)) {
        plot(c(1,2), c(0,3), type = "n", xlab = "", ylab = "")
        for(i in seq(3)) {
          lines(c(1,2), rep(i-1,2), col = sex.col[i], lwd = 2)
          text(1, i-0.5, sex.names[i], col = sex.col[i], adj = 0)
        }
      }
    }
  }
  invisible()
}
###########################################################  
is.raws <- function(filename = NULL)
{
  if(is.null(filename))
    return(TRUE)

  ## This assumes files have traits in rows.
  new.names <- scan(filename, "", sep = ",", nlines = 1, quiet = TRUE)
  is.names <- new.names[grep("^Mouse", new.names)]
  if(length(is.names))
    return(TRUE)

  is.names <- new.names[grep("^rs", new.names)]
  if(length(is.names))
    return(FALSE)

  ## Otherwise we need to look at first row.
  new.names <- scan(filename, "", sep = ",",
                    nlines = 1, skip = 1, n = 1, quiet = TRUE)
  if(pmatch("Mouse", new.names, nomatch = 0))
    return(TRUE)
  if(pmatch("rs", new.names, nomatch = 0))
    return(FALSE)

  stop("cannot determine if file contains LOD profiles or raw data")
  invisible()
}
###########################################################  
get.raws <- function(filename = NULL, newdata = NULL)
{
  if(!is.null(filename)) {
    newdata <- read.csv(filename, header = TRUE)
    new.names <- names(newdata)
    mouse.names <- new.names[grep("^Mouse", new.names)]

    if(length(mouse.names) > 1) {
      ## Need to transpose.

      ## Check if there is a separate tissue column.
      tmp <- match("Tissue", new.names)
      if(!is.na(tmp)) {
        tissues <- newdata[[tmp]]
        newdata <- newdata[, -tmp]
        new.names <- names(newdata)
      }
      else
        tissues <- NULL

      if(length(mouse.names) < length(new.names)) {
        ## First column has MouseNum.
        new.names <- as.character(newdata[[1]])

        ## Add tissue to newdata names.
        if(!is.null(tissues))
          new.names <- paste(tissues, new.names, sep = ".")
        
        tmp <- t(newdata[, mouse.names])
        dimnames(tmp)[[2]] <- new.names
        newdata <- data.frame(MouseNum = mouse.names)
        newdata <- cbind(newdata, tmp)
        
      }
      else {
        mouse.names <- new.names
        new.names <- row.names(newdata)

        ## Add tissue to newdata names.
        if(!is.null(tissues))
          new.names <- paste(tissues, new.names, sep = ".")
        
        newdata <- t(newdata)
        dimnames(newdata) <- list(NULL, new.names)
        newdata <- cbind(data.frame(MouseNum = mouse.names), newdata)
      }
    }
  }
  newdata
}  
###########################################################  
add.phenos <- function(cross, newdata = NULL, index = NULL)
{
  if(!is.null(newdata)) {
    if(any(names(newdata) %in% names(cross$pheno)))
      warning("some cross phenotypes overwritten with new data")

    if(is.null(index)) {
      n.ind <- nind(cross)
      if(nrow(newdata) != n.ind)
        stop(paste("newdata must have number of rows (",
                   nrow(newdata), ") as cross individuals (",
                   n.ind, ")", sep = ""))

      ## Must assume newdata in same order as cross here.
      for(i in names(newdata))
        cross$pheno[[i]] <- newdata[[i]]
    }
    else {
      ## The row.names of newdata must be index.

      mat <- match(row.names(newdata), cross$pheno[[index]], nomatch = 0)
      if(!any(mat > 0))
        stop("no row names of newdata match index")

      tmp <- rep(NA, nind(cross))
      for(i in names(newdata)) {
        tmp[mat] <- newdata[mat > 0, i]
        cross$pheno[[i]] <- tmp
      }
    }
  }
  cross
}
###########################################################  
add.raws <- function(cross, newdata = NULL,
                     traitnames = names(newdata),
                     covariates = NULL,
                     index = "MouseNum")
{
  ## Allow user to pass a data frame
  if(!is.null(newdata)) {
    ## Expand trait names if abbeviated.
    tmp3 <- mytrait(c(traitnames,covariates), newdata, warn = FALSE)
    tmp3 <- tmp3[!is.na(tmp3)]

    ## Add new phenotypes and/or covariates.
    if(length(tmp3)) {
      ## First column of newdata must have MouseNum as "Mousemmmm".
      row.names(newdata) <- newdata[[1]]
      newdata <- newdata[, tmp3]
      cross <- add.phenos(cross, newdata, index)
    }
    else
      warning("no traitnames or covariates match newdata")
  }
  cross
}
