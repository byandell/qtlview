################################################################
## Create aug.scanone object.
## This is augmented scanone object.
aug.scanone <- function(traitnames = mytrait(pheno = cross$pheno),
                        sex = sexes,
                        method = "hk",
                        lod.error = 100,
                        cross = get(cross.name),
                        cross.name = deparse(substitute(cross)),
                        normal.scores = TRUE,
                        intcov = NULL,
                        addcov = NULL,
                        trait.annotation = NULL,
                        maps = myget(cross.name, "maps"),
                        scan.type = c("LOD","LPD","2logBF"),
                        ...)
{
  sexes <- c("both","male","female","ignore")
  sex <- match.arg(sex, sexes)

  ## Allow for qb.scanone or scanone
  scan.type <- match.arg(scan.type)
  if(scan.type == "LOD") {
    if(is.null(attr(cross$geno[[1]]$prob, "step")))
      stop("Need to first run calc.genoprob")

    scan.type <- "lod"
    scan.prog <- scanone
  }
  else {
    if(is.null(attr(cross$geno[[1]]$prob, "step")))
      stop("Need to first run qb.genoprob")

    method <- scan.type
    scan.prog <- scanone.qb
  }

  n.traits <- length(traitnames)

  ## Set up data frame to store LODs.
  lods <- pull.pseudomarkers(cross, traitnames = traitnames, ...)
  
  switch(sex,
         both =, ignore = {
         },
         male = {
           cross <- subset(cross, ind = (cross$pheno$Sex == "Male"))
         },
         female = {
           cross <- subset(cross, ind = (cross$pheno$Sex == "Female"))
         })
  if(sex == "both")
    intcov <- cbind(intcov, model.matrix(~Sex, cross$pheno)[,-1, drop = FALSE])

  ## Here assume addcov and intcov are distinct.
  addcov <- cbind(addcov, intcov)

  normal.scores <- array(normal.scores, n.traits)
  
  for(i in seq(n.traits)) {
    traitname <- traitnames[i]
    cat(traitname, "...\n")
    if(is.logical(cross$pheno[[traitname]]))
      cross$pheno[[traitname]] <-
        as.numeric(cross$pheno[[traitname]])
    if(is.factor(cross$pheno[[traitname]]))
      cross$pheno[[traitname]] <-
        as.numeric(cross$pheno[[traitname]]) - 1
    if(normal.scores[i])
      cross$pheno[[traitname]] <- normal.trans(cross$pheno[[traitname]])

    pheno.col <- find.pheno(cross, traitname)

    if(!all(is.na(cross$pheno[[pheno.col]]))) {
      ## Profile LOD.
      tmp <- scan.prog(cross, , pheno.col, method = method,
                     intcov = intcov, addcov = addcov, ...)
      lods[[traitname]] <- tmp[[scan.type]]
    }
  }

  ## Recode really large lods with NA.
  lods.attr <- attributes(lods)
  tmp <- lods[, traitnames]
  tmp[tmp > lod.error & !is.na(tmp)] <- NA
  lods[, -(1:2)] <- tmp
  ## Recover attributes that got trashed.
  for(i in names(lods.attr)) {
    if(is.null(attr(lods, i)))
       attr(lods, i) <- lods.attr[[i]]
  }

  attr(lods, "method") <- method

  ## Add trait positions.
  trait.position <-
    find.trait.position(traitnames, maps, trait.annotation, ...)
  attr(lods, "trait.position") <- trait.position
    
  class(lods) <- c("aug.scanone", "scanone", "data.frame")
  
  lods
}
################################################################
scanone.qb <- function(cross, chr, pheno.col, method = "LPD",
                       intcov = intcov, addcov = addcov,
                       ...,
                       interval = rep(30, length(cross$geno)),
                       verbose = FALSE)
{
  ## Wrapper for qb.scanone.
  require(qtlbim)

  ## Cannot handle X chr yet. Fake it.
  ## Need to create rows for X that are all set to zero.
  ## Assume X chr is at the end of list.
  is.X <- sapply(cross$geno, class) == "X"
  if(any(is.X)) {
    loci <- pull.pseudomarkers(cross, ...)
    chr.names <- names(cross$geno)
    loci <- loci[loci$chr %in% chr.names[is.X], ]
    for(i in chr.names[is.X])
      cross$geno[[i]] <- NULL
    ## And see below for the wrap.
  }
  
  ## Covariates are handled differently from scanone.
  ## Assume intcov is subset of addcov from calling routine.
  ## Also assume addcov columns are numeric already.
  
  ## Set up sex covariate centered on 0 w.r.t. missing values for trait.
  for(i in seq(ncol(addcov))) {
    addcov.name <- paste("covar", i, sep = "")
    if(match(addcov.name, names(cross$pheno)))
      stop(paste("Some phenotype is named", addcov.name))
    
    ## Center covariates to improve mixing.
    cross$pheno[[addcov.name]] <- (addcov[[i]] - mean(addcov[[i]], na.rm = TRUE)) /
      sd(addcov[[i]], na.rm = TRUE)
  }
  fixcov <- find.pheno(cross, paste("covar", seq(ncol(addcov)), sep = ""))
  intcovs <- !is.na(match(dimnames(addcov)[[2]], dimnames(intcov)[[2]]))

  ## Call qtlbim routines. These take some time!
  tmp <- qb.mcmc(cross, traitname, normal.scores[i], cross.name = cross.name,
                 cross = cross, fixcov = fixcov, intcov = intcovs, interval = interval,
                 ..., verbose = verbose)
  out <- qb.scanone(tmp, type = method)
  
  if(any(is.X)) {
    ## Append rows with 0 values.
    loci <- cbind(loci, matrix(0, nrow(loci), ncol(out) - 2))
    names(loci)[-(1:2)] <- names(out)[-(1:2)]
    loci$chr <- as.character(loci$chr)
    out$chr <- as.character(out$chr)
    out <- rbind(out, loci)
    out$chr <- ordered(out$chr, chr.names)
  }
  out
}
################################################################
max.aug.scanone <- function(object,
                             chr = levels(object$chr),
                             traitnames = names(object)[-(1:2)],
                             threshold.lod = 0,
                             mean.peaks = FALSE,
                             ...)
{
  object.class <- c("scanone", "data.frame")
  maps <- attr(object, "maps")
  trait.position <- attr(object, "trait.position")

  ## Reduce to selected chromosomes.
  object <- object[object$chr %in% chr, ]

  traitnames <- names(object)[match(traitnames, names(object), nomatch = 0)]
  object <- object[, c(TRUE, TRUE, names(object)[-(1:2)] %in% traitnames)]
  
  ## Drop traits below threshold.lod.
  if(any(threshold.lod > 0)) {
    tmp <- object[, -(1:2)]
    object[, -(1:2)] <- tmp[, threshold.pass(tmp, threshold.lod, object$chr)]
  }
  else { ## Drop any that are all NA.
    tmp <- apply(object, 2, function(x) !all(is.na(x)))
    object <- object[, tmp, drop = FALSE]
  }
  if(!nrow(object))
    return(NULL)
 
  ## Set class to scanone.
  class(object) <- object.class
  
  n.chr <- length(chr)
  n.traits <- ncol(object) - 2
  traitnames <- names(object)[-(1:2)]
  out <- data.frame(chr = ordered(rep(0, n.traits), levels(object$chr)),
                    pos.cM = rep(0, n.traits),
                    pos.Mb = rep(0, n.traits),
                    lod = rep(0, n.traits))
  row.names(out) <- traitnames

  for(i in traitnames) {
    tmp <- max(object[ , c("chr","pos",i)], ...)
    if(nrow(tmp) > 1) {
      ## Multiple ties to max lod.
      tmp2 <- tmp[!duplicated(tmp$chr),]
      tmpfn <-  {
        if(mean.peaks)
          mean
        else
          function(x) ifelse(length(x) == 1, x, median(x))
      }
      tmp2$pos <- unlist(tapply(tmp$pos, tmp$chr, tmpfn))[tmp2$chr]
      tmp <- tmp2
    }
    if(nrow(tmp)) {
      tmp$pos <- round(tmp$pos, 3)

      ## Add position in Mb.
      tmp$pos.Mb <- rep(NA, nrow(tmp))
      tmp <- tmp[, c(1:2,4,3)]
      for(j in levels(tmp$chr)) {
        jj <- tmp$chr == j
        if(any(jj))
          tmp$pos.Mb[jj] <- qm.approx(maps, "cM", j, tmp$pos[jj])$y
      }
      out[i, ] <- tmp[1,]
    }
  }

  ## Add trait positions in cM and Mb.
  add.trait.position(out, trait.position)
}
################################################################
print.aug.scanone <- function(x, ...) print(summary(x, ...), ...)
################################################################
summary.aug.scanone <- function(object,
                                chr = levels(object$chr),
                                traitnames = names(object)[-(1:2)],
                                threshold.level = 0,
                                threshold.lod = 0,
                                mean.peaks = FALSE,
                                hc = NULL, n.clust = 0, comp.lod = NULL,
                                summary.plot = "none",
                                long = FALSE,
                                col = c("black","blue","red","green","purple","magenta"),
                                heatmap = FALSE,
                                ...)
{
  ylab <- object$ylab
  if(long) {
    object.class <- c("scanone", "data.frame")

    ## Reduce to selected chromosomes.
    object <- object[object$chr %in% chr, ]

    traitnames <- names(object)[match(traitnames, names(object),
                                      nomatch = 0)]
    object <-
      object[, c(TRUE, TRUE, names(object)[-(1:2)] %in% traitnames)]
    
    ## Drop in traits below threshold.lod.
    if(any(threshold.lod > 0)) {
      tmp <- object[, -(1:2)]
      object[, -(1:2)] <- tmp[, threshold.pass(tmp, threshold.lod, object$chr)]
    }
    else { ## Drop any that are all NA.
      tmp <- apply(object, 2, function(x) !all(is.na(x)))
      object <- object[, tmp, drop = FALSE]
    }
    
    ## Set class to scanone.
    class(object) <- object.class
    
    n.chr <- length(chr)
    n.traits <- ncol(object) - 2
    traitnames <- names(object)[-(1:2)]
    
    out <- list(lod = matrix(NA, n.traits, n.chr))
    dimnames(out$lod) <- list(traitnames, chr)
    out$pos <- out$lod

  
    for(i in traitnames) {
      tmp <- sumone.scanone(object[ , c("chr","pos",i)], chr = chr, ...)

      out$pos[i, as.character(tmp$chr)] <- tmp$pos
      out$lod[i, as.character(tmp$chr)] <- tmp[[i]]
    }
    out$loci <- object[,1:2]
  }
  else {
    maps <- attr(object, "maps")
    
    ## Print max lod by trait.
    maxout <- max(object, chr = chr, threshold.lod = 0, ...)
    tmp2 <- nrow(maxout)
    if(!is.null(hc)) {
      ## Order traits for subset that are plotted.
      ## Order other traits by decreasing lod.
      is.selected <- match(rev(hc$labels[hc$order]), dimnames(maxout)[[1]])
      if(n.clust > 1) {
        maxout <- cbind(maxout, cluster = 0)
        
        ## Assign cluster numbers from n.clust on down.
        tmp <- NULL
        for(i in 2:n.clust) {
          tmp2 <- rev(cutree(hc, i)[hc$order])
          tmp2 <- i + 1 - c(unclass(ordered(tmp2, unique(tmp2))))
          tmp <- paste(tmp, tmp2, sep = ifelse(i < 10, "", "."))
        }
        maxout[is.selected, "cluster"] <- tmp
        maxout$cluster <- factor(maxout$cluster)
      }
    }
    else
      ##***Modify for X chromosome.
      is.selected <- seq(tmp2)[maxout$lod > threshold.lod[1]]

    ## Add phenotype name as column for printing.
    ## Use index as row.names to track sorting
    ## of those selected above threshold.
    maxout <- cbind(phenotype = row.names(maxout), maxout)
    row.names(maxout) <- NULL

    out <- list(selected = maxout[is.selected, ],
                threshold.level = threshold.level,
                threshold.lod = threshold.lod,
                trait.position = object$trait.position,
                maps = maps)

    tmp2 <- nrow(maxout)
    if(length(is.selected) < tmp2) {
      tmp2 <- seq(tmp2)[-is.selected]
      if(!is.null(hc))
        tmp2 <- tmp2[order(-maxout[tmp2,"lod"])]
      out$rest <- maxout[tmp2, ]
    }

    out$comp.lod <- comp.lod
  }
  attr(out, "long") <- long

  col <- array(col, length(traitnames))

  attr(out, "color") <- col
  attr(out, "heatmap") <- heatmap
  attr(out, "ylab") <- ylab

  class(out) <- c("summary.aug.scanone", "list")
  out
}
################################################################
plot.summary.aug.scanone <- function(x,
                                      chr = dimnames(x$lod)[[2]],
                                      threshold.lod = 0,
                                      max.lod = 20,
                                      max.names = 100,
                                      main = "",
                                      by = c("chr","trait","phenotype"),
                                      scale = c("cM","Mb"),
                                      cex = 2, pch = 3,
                                      use.cM = FALSE,
                                      ...)
{
  ylab.choice <- attr(x, "ylab")
  
  by <- match.arg(by)
  if(length(chr) == 1)
    by <- "trait"
  scale <- match.arg(scale)

  n.chr <- length(chr)

  maps <- x$maps

  ## Subset to chr.
  x$lod <- as.matrix(x$lod[, chr, drop = FALSE])

  tmp <- threshold.pass(x$lod, threshold.lod, chr, 1)
  x$lod <- as.matrix(x$lod[tmp,, drop = FALSE])
  pos <- as.matrix(x$pos[tmp, chr])
  n.traits <- sum(tmp)
  if(n.traits == 0)
    stop("no phenotypes pass threshold.lod")
  
  ## Limit lod to (0, max.lod).
  x$lod[is.na(x$lod)] <- 0
  x$lod[x$lod < 0] <- 0
  max.lod <- min(max(x$lod, na.rm = TRUE), max.lod)
  x$lod[x$lod > max.lod] <- max.lod
  traitnames <- dimnames(x$lod)[[1]]

  ## Translate to Mb if desired.
  if(scale == "Mb" & by == "chr") {
    for(i in chr)
      pos[,i] <- qm.approx(maps, "cM", i, pos[,i])$y
  }
  pos <- c(pos)
  
  xlab <- paste("Position (", scale, ")", sep = "")
  if(by == "chr") {
    n.y <- n.chr
    n.z <- n.traits
    reps <- rep(n.z, n.y)
    ats <- seq(n.y)
    labels <- chr
    ylab <- "Chromosome"
    tmpar <- par(mar = c(3, 4, 2 + 2 * (length(chr) == 1), 0) + 0.1)
#    on.exit(par(tmpar))
    plot(range(pos, na.rm = TRUE), c(0, n.y) + 0.5, type = "n",
         xlab = "", ylab = ylab, yaxt = "n")
    mtext(xlab, 1, 2)
    abline(h = seq(0, n.y) + 0.5, col = "gray", lty = 3, lwd = 2)
    tmp <- rep(seq(n.y), reps)
    points(pos, jitter(tmp), cex = cex * c(x$lod) / max.lod, pch = pch)
    axis(side = 2, at = ats, labels = labels, las = 1)
  }
  else { ## By phenotype = trait.
    ## This only makes sense now for length(chr) == 1.
    ## Want to show scale like for plot.aug.scanone.
    n.y <- n.traits
    n.z <- n.chr
    reps <- n.z
    add.names <- n.traits <= max.names
    ats <- if(add.names)
      seq(n.traits)
    else
      pretty(seq(n.traits))
    labels <- mylabels(traitnames, max.names, ylab.choice)
    ylab <- ifelse(add.names, "", "Phenotype Index")
    tmpar <- par(mar = c(3, 5 + 3 * add.names, 2 + 2 * (n.chr == 1), 0) + 0.1)
    on.exit(par(tmpar))
    if(n.chr == 1) {
      plot(range(maps$cM.map[[chr]]), c(0, n.y) + 0.5, type = "n",
           xlab = "", ylab = ylab, yaxt = "n")
      abline(h = seq(0, n.y) + 0.5, col = "gray", lty = 3, lwd = 2)
      tmp <- rep(seq(n.y), reps)
      points(pos, jitter(tmp), cex = cex * c(x$lod) / max.lod, pch = pch)
      axis(side = 2, at = ats, labels = labels, las = 1)

      rug(maps$cM.map[[chr]], 0.02, quiet = TRUE, side = 1)
      add.rug(chr, "", maps)
    }
    else { ## n.chr > 1
      n.loci <- nrow(x$loci)
      xlab <- "Chromosome"
      loci <- x$loci[x$loci$chr %in% chr,]
      chr.offset <- unlist(tapply(loci$pos, loci$chr, max))
      chr.offset[is.na(chr.offset)] <- 0
      chr.offset <- c(0, cumsum(chr.offset + 5)[chr])
      plot(range(chr.offset), c(0, n.y) + 0.5, type = "n",
           xlab = "", ylab = ylab, yaxt = "n", xaxt = "n")
      mtext(xlab, 1, 2)
      abline(h = seq(0, n.y) + 0.5, col = "gray", lty = 3, lwd = 2)
      pos <- pos + chr.offset[rep(seq(n.chr), rep(n.traits, n.chr))]
      tmp <- rep(seq(n.y), reps)
      points(pos, jitter(tmp), cex = cex * c(x$lod) / max.lod, pch = pch)
      axis(side = 2, at = ats, labels = labels, las = 1)
      abline(v = chr.offset - 2.5)
      axis(1, at = (chr.offset[-1] + chr.offset[-n.chr - 1]) / 2,
           labels = chr)
    }
  }
  invisible(title(main))
}
###############################################################################
print.summary.aug.scanone <- function(x, digits = 3, ...)
{
  if(attr(x, "long")) {
    cat("LOD peak positions\n")
    print(signif(x$pos, digits))
    cat("\nLOD peak values\n")
    invisible(print(signif(x$lod, digits)))
  }
  else {
    tmpopt <- options(width = 133)
    on.exit(options(tmpopt))

    heatmap <- attr(x, "heatmap")
    if(!heatmap)
      x$selected$color <- attr(x, "color")[as.numeric(row.names(x$selected))]
    print(x$selected)
    
    if(any(x$threshold.lod > 0))
      cat("\nThreshold:", paste(round(x$threshold.lod, 2), collapse = ","),
          "level:", x$threshold.level, "\n")
    
    if(!is.null(x$rest)) {
      cat("\n")

      if(!heatmap)
        x$rest$color <- attr(x, "color")[as.numeric(row.names(x$rest))]
      print(x$rest)
    }

    if(!is.null(x$comp.lod)) {
      cat("\nSummary of LOD average across all traits:\n\n")
      print(mysum.scanone(x$comp.lod, threshold.lod = 0, ...))
    }
  }
}
###############################################################################
plot.aug.scanone <- function(x,
                             chr = levels(x$chr),
                             traitnames = names(x)[-(1:2)],
                             col.scheme = c("redblue", "gray", "cm", "heat",
                               "terrain", "topo"),
                             gamma = 0.6, allow.neg = FALSE,
                             max.names = 100,
                             zscale = ifelse(add.names, "value", "gray"),
                             main = "",
                             threshold.lod = 0,
                             max.lod = 10,
                             rescale = TRUE,
                             cluster = TRUE,
                             add.position = TRUE,
                             cex = cexs,
                             lwd = lwds,
                             support.lod = 1.5,
                             hc = NULL,
                             n.clust = 0,
                             beta = 1,
                             use.cM = FALSE,
                             heatmap = TRUE,
                             ylab = c("symbol","a_gene_id","symbol.a_gene_id","none"),
                             lod.name = "LOD",
                             ...) 
{
  ylab.choice <- match.arg(ylab)

  if(!heatmap)
    return(myplot.scanone(x, threshold.lod, main, chr = chr,
                          lod.name = lod.name, ...))
  if(is.null(support.lod))
    support.lod <- 0

  ## Code borrowed from plot.scantwo in R/qtl. (But Yandell wrote part of it.)
  col.scheme <- match.arg(col.scheme)
  
  cols <- switch(col.scheme,
                 gray = {
                   if (gamma <= 0) 
                     rev(gray(seq(0, 1, len = 256)))
                   else rev(gray(log(seq(1, exp(gamma), len = 256))/gamma))
                 }, 
                 heat = heat.colors(256),
                 terrain = terrain.colors(256), 
                 topo = topo.colors(256), cm = cm.colors(256),
                 redblue = rev(rainbow(256, 
                   start = 0, end = 2/3, gamma = gamma)))

  traitnames <- names(x)[match(traitnames, names(x), nomatch = 0)]

  lod <- x[, traitnames, drop = FALSE]
  
  ## Make values positive.
  if (!allow.neg && any(!is.na(lod) & lod < -1e-06)) {
    u <- !is.na(lod) & lod < 0
    n <- sum(u)
    warning(n, " LOD scores <0, set to 0")
    lod[u] <- 0
  }
  if (any(!is.na(lod) & lod == Inf)) {
    u <- !is.na(lod) & lod == Inf
    n <- sum(u)
    warning(n, " LOD scores =Inf, set to NA")
    lod[u] <- NA
  }

  ## Reduce to traits with max(lod) > threshold.lod.
  lod <- lod[, threshold.pass(lod, threshold.lod, x$chr), drop = FALSE]

  if(ncol(lod) == 0) {
    stop(paste("No traits pass threshold of", threshold.lod))
    return(myplot.scanone(x,
                          threshold.lod = threshold.lod, main = main,
                          chr = chr, ...))
  }
  
  n.traits <- ncol(lod)
  if(n.traits == 0)
    stop("no phenotypes with LOD above threshold.lod")

  traitnames <- dimnames(lod)[[2]]
  labels <- mylabels(traitnames, length(traitnames), ylab.choice)

  lod <- lod[x$chr %in% chr,]
  
  if(is.logical(rescale))
    rescale <- "support"
  else
    rescale <- match.arg(rescale, c("support","peaks","none"))
  
  if(rescale != "none") {
    ## Rescale to 10 + LOD - max(LOD).
    ## Perhaps do this per chromosome above threshold?
    if(rescale == "support") {
      if(support.lod > 0)
        max.lod <- 3 * support.lod
      tmpfn <- function(x) {
        maxx <- max(x, na.rm = TRUE)
        ## Uses autosome threshold.
        if(maxx > threshold.lod[1])
          pmax(0, max.lod + x - maxx)
        else
          rep(0, length(x))
      }
      tmpfn2 <- function(x, chr) unlist(tapply(x, chr, tmpfn))
      plod <- apply(as.matrix(lod), 2, tmpfn2, x$chr[x$chr %in% chr])
    }
    else { ## rescale == "peaks"
      tmpfn <- function(x) {
        maxx <- max(x, na.rm = TRUE)
        x / maxx
      }
      plod <- apply(as.matrix(lod), 2, tmpfn)
    }
  }
  else {
    max.lod <- min(max(as.matrix(lod), na.rm = TRUE), max.lod)
    plod <- pmin(as.matrix(lod), max.lod)
  }

  if(!is.null(hc)) {
    lod <- lod[, hc$order]
    plod <- plod[, hc$order]
    if(!add.position)
      n.clust <- 0
    if(n.clust > 1) {
      clust <- cutree(hc, n.clust)[hc$order]
      clust <- seq(clust)[!duplicated(clust)][-1] - 0.5
    }
    else
      n.clust <- 0
  }
  else
    n.clust <- 0

  maps <- attr(x, "maps")

  ## Set up cM to Mb translation.  
  if(length(chr) == 1) {
    if(use.cM)
      p <- qm.approx(maps, "Mb", chr)
    else ## use Mb.
      p <- qm.approx(maps, "cM", chr)
  }

  ## Make room for Z scale.
  add.names <- n.traits <= max.names
  dots <- list(...)
  if(is.logical(zscale))
    zscale <- "gray"
  if (zscale == "gray") {
    if ("layout" %in% names(dots)) 
      layout(dots[["layout"]][[1]], dots[["layout"]][[2]])
    else layout(cbind(1, 2), c(ifelse(rescale != "none", 40, 10), 1))
    lo <- ifelse(allow.neg, -max.lod, 0)
  }
  if(add.names & ylab.choice != "none")
    par(cex = 0.75)
  tmpar <- par(mar = {
    if ("mar1" %in% names(dots)) 
      dots[["mar1"]]
    else
      c(4,
        5 + round(add.names * max(1 + nchar(labels)) / 5),
        3 + (length(chr) == 1),
        (zscale == "value") * 3) + 0.1
  })

  ylab <- ifelse(add.names, "", "Phenotype Index")

  ## Set up to plot trait position if it is gene transcript.
  trait.position <- attr(x, "trait.position")
  if(add.position & !is.null(trait.position)) {
    tmp <- attr(trait.position, "prefix")
    trait.position <- trait.position[dimnames(lod)[[2]], ]
    if(tmp != "")
      names(trait.position) <- substring(names(trait.position),
                                         nchar(tmp) + 2)
    
    n.pos <- nrow(trait.position)
    cexs <- max(0.5, min(2, 50 / n.pos))
    lwds <- max(1, min(2, 50 / n.pos))
  }
  else
    n.pos <- 0

  if(support.lod > 0) {
    col.breaks <- c("gray","orange","red")
    breaks <- rev(seq(max.lod + 0.00001, 0, by = -support.lod))
    breaks[1] <- 0
  }
  if (length(chr) > 1) {
    if(rescale == "support" & support.lod > 0) {
      image(1:nrow(lod), 1:n.traits , plod,
            ylab = ylab, 
            xlab = "",
            col = col.breaks,
            breaks = breaks,
            xaxt = "n", yaxt = "n")
      if(n.clust)
        abline(h = clust, col = "darkgray")
    }
    else {
      image(1:nrow(lod), 1:n.traits , plod,
            ylab = ylab, 
            xlab = "", col = cols, 
            xaxt = "n", yaxt = "n")
    }
    mtext("Chromosome", 1, 2)
    n.mar <- rep(0, length(chr))
    for (i in 1:length(chr))
      n.mar[i] <- sum(x$chr == chr[i])
    wh <- c(0.5, cumsum(n.mar) + 0.5)
    abline(v = wh, xpd = FALSE)
    a <- par("usr")
    abline(v = a[1:2], xpd = FALSE)
    abline(h = a[3:4], xpd = FALSE)
    for (i in 1:length(n.mar))
      axis(side = 1, at = mean(wh[i + c(0, 1)]), labels = chr[i])

    ## Add circle at gene if on this chromosome.
    if(n.pos) {
      names(wh) <- c(chr, "")
      for(chri in as.character(chr)) {
        tmp <- (trait.position$chr == chri) & !is.na(trait.position$cM)
        if(any(tmp)) {
          posi <- outer(trait.position$cM[tmp], x$pos[x$chr == chri],
                        function(x,y) abs(x-y))
          posi <- apply(posi, 1, which.min)
          points(wh[chri] + posi, seq(n.pos)[tmp],
                 cex = cex, lwd = lwd, col = ifelse(col.scheme == "gray", "blue", "black"))
        }
      }
    }
  }
  else {
    if(use.cM)
      ipos <- x$pos[x$chr %in% chr]
    else
      ipos <- qm.approx(maps, "cM", chr, x$pos[x$chr %in% chr])$y

    if(rescale == "support" & support.lod > 0) {
      image(ipos, 1:n.traits, plod,
            ylab = ylab, 
            xlab = "",
            col = col.breaks,
            breaks = breaks,
            xaxt = "n", yaxt = "n")
    }
    else {
      image(ipos, 1:n.traits, plod,
            ylab = ylab, 
            xlab = "",
            col = cols, yaxt = "n", xaxt = "n")
    }
    if(n.clust)
      abline(h = clust, col = "darkgray")
    add.rug(chr, "", maps, use.cM = use.cM, outer = TRUE,
            xlim = range(ipos), bottom.axis = TRUE)

    ## Add circle at gene if on this chromosome.
    units <- ifelse(use.cM,"cM","Mb")
    if(n.pos) {
      tmp <- (trait.position$chr == chr) & !is.na(trait.position[[units]])
      if(any(tmp))
        points(trait.position[[units]][tmp], seq(n.pos)[tmp],
               cex = cex, lwd = lwd, col = ifelse(col.scheme == "gray", "blue", "black"))
    }
  }

  ## Vertical axis labels.
  if(add.names)
    axis(side = 2, at = seq(n.traits),
       labels = mylabels(dimnames(lod)[[2]], max.names, ylab.choice),
       las = 1)
  else
    axis(side = 2, at = pretty(seq(n.traits)), labels = TRUE, las = 1)

  title(main, line = 1 + 2 * (length(chr) == 1))

  if(zscale %in% c("gray","value")) {
    ## Get max lod over displayed genome.
    lod.max <- apply(lod, 2, max, na.rm = TRUE)
    ## Rescale proportions to percents (hopefully).
    if(max(lod.max) <= 1)
      lod.max <- lod.max * 100
  }
  if(zscale == "value") {
    mtext(lod.name, 4, line = 0.5, las = 1, at = par("usr")[4])
    axis(round(lod.max, 1), side = 4, at = seq(n.traits), las = 1)
#    mtext(round(lod.max, 1), side = 4, at = seq(n.traits), las = 1, line = 0.5)
  }
  else if (zscale == "gray") {
    tmp <- par("mar")
    tmp[2] <- ifelse(rescale != "none", 0.1, 2.1)
    tmp[4] <- 0.1
    tmpar <- par(mar = {
      if ("mar2" %in% names(dots)) 
        dots[["mar2"]]
      else
        tmp
    })

    ## Add Scale for LOD using log(lod) on col scheme.
    if(rescale != "none") {
      image(1:1, 1:n.traits,
            t(as.matrix(rank(lod.max))),
            ylab = "", 
            xlab = "",
            col = rev(gray(seq(0, 1, len = 256))), yaxt = "n", xaxt = "n")
      tmp <- par("usr")
      abline(h = tmp[3:4], v = tmp[1:2])
      if(n.clust)
        abline(h = clust, col = "darkgray")
      mtext("LOD", 3, 0)
    }
    else {
      colorstep <- (max.lod - lo)/255
      image(x = 1:1, y = seq(lo, max.lod, colorstep),
            z = matrix(c(1:256), 1, 256),
            zlim = c(1, 256), ylab = "", xlab = "", 
            xaxt = "n", yaxt = "n", col = cols)
      u <- par("usr")
      abline(v = u[1:2], xpd = FALSE)
      abline(h = u[3:4], xpd = FALSE)
      yloc <- pretty(c(lo, max.lod), 4)
      yloc <- yloc[yloc >= u[3] & yloc <= u[4]]
      axis(side = 2, at = yloc, labels = yloc)
    }
    par(tmpar)
  }
  invisible(hc)
}


