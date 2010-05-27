################################################################
## My approximation routine. Use qm.approx, hide myapprox.
qm.approx <- function(maps, base = bases, chr, ..., non.seg = FALSE)
{
  bases <- c("cM","Mb")
  base <- pmatch(base, bases)
  if(is.na(base))
    stop("base must be cM or Mb")
  
  x <- bases[base]
  y <- bases[-base]
  non.seg <- ifelse(non.seg, "same", "map")
  x <- paste(x, non.seg, sep = ".")
  y <- paste(y, non.seg, sep = ".")
  myapprox(maps[[x]][[chr]], maps[[y]][[chr]], ...)
}
################################################################
## This is the only plot routine that refers to same and map.
add.rug <- function(chr, main, maps,
                    p = qm.approx(maps, off.base, chr),
                    use.cM,
                    outer = FALSE,
                    xlim = range(map),
                    bottom.axis = FALSE,
                    side = 1)
{
  ##*** Somehow the top axis is getting cM and Mb mixed up. Need to fix.
  bases <- c("cM","Mb")
  base <- bases[2 - use.cM]
  off.base <- bases[1 + use.cM]
  
  ## Add rugs, etc.
  if (length(chr) == 1) {
    ## Get map for chr using proper base.
    map <- maps[[paste(base, "map", sep = ".")]][[chr]]
    ticksize <- ifelse(outer, -0.02, 0.02)

    ## Get plot limits in plotting units.
    usr <- par("usr")

    ## Add grey ticks for non-segregating markers (if available).
    non.seg <- maps[[paste(base, "same", sep = ".")]]
    if(!is.null(non.seg)) {
      non.seg <- non.seg[[chr]]
      rug(non.seg, 0.75 * ticksize, quiet = TRUE, side = side, col = "gray")
      rug(non.seg, 0.75 * ticksize, quiet = TRUE, side = side + 2, col = "gray")
      if(side == 1)
        abline(h = usr[3:4])
      else
        abline(v = usr[1:2])
    }

    rug(map, ticksize, quiet = TRUE, side = side)
    if(bottom.axis) {
      axis(side, pretty(xlim, n = 30), line = ifelse(outer, 0.6, 0))
    }

    rug(map, ticksize, quiet = TRUE, side = side + 2)
    ## This is the culprit.
    axis(side + 2, p$y, p$x, line = ifelse(outer, 0.6, 0))

    usr <- usr[2 * side - c(1,0)]
    usr <- usr[1] - 0.01 * diff(usr[1:2])
    if(use.cM) {
      mtext("cM", side,     1.6, at = usr, adj = 1)
      mtext("Mb", side + 2, 1.6, at = usr, adj = 1)
    }
    else {
      
      mtext("cM", side + 2, 1.6, at = usr, adj = 1)
      mtext("Mb", side,     1.6, at = usr, adj = 1)
    }
    mtext(paste("Chromosome", chr), side, 1.35 + outer)
  }
  title(main, line = 0.5 + 2 * (length(chr) == 1))
}
################################################################
myapprox <- function(Mb, cM,
  pos = posn, n.pos = 30, ...)
{
  ## Translate Mb to cM within range.
  
  ## Some wierd bug because Mb is of class "A" or "X", but not "numeric".
  posn <- pretty(c(Mb), n.pos)

  ## Adjust pos to be within Mb range.
  tmp <- c(pos)
  tmp <- pmin(max(c(Mb)), pmax(min(c(Mb)), tmp))

  ## Linear interpolation between SNPs.
  require(stats)
  p <- approx(c(Mb), c(cM), tmp)

  ## Reset x to be pos.
  p$x <- pos
  p
}
################################################################
pull.pseudomarkers <- function(cross,
                               step = step.default,
                               off.end = off.end.default,
                               stepwidth = stepwidth.default,
                               traitnames = NULL, ...)
{
  if(is.null(cross$geno[[1]]$prob)) {
    step.default <- 2
    off.end.default <- 0
    stepwidth.default <- "variable"
  }
  else {
    step.default <- attr(cross$geno[[1]]$prob, "step")
    off.end.default <- attr(cross$geno[[1]]$prob, "off.end")
    stepwidth.default <- attr(cross$geno[[1]]$prob, "stepwidth")
  }
  ## Get loci = SNPs plus pseudomarkers between SNPs.
  ## Assume here that cross already run through calc.genoprob
  ##        or other arguments provided.
  if(is.null(step) | is.null(off.end) | is.null(stepwidth))
    cross <- calc.genoprob(cross)
  
  tmpfn <- function(x, step, off.end, stepwidth)
    create.map(x$map, step, off.end, stepwidth)

  loci <- lapply(cross$geno, tmpfn, step, off.end, stepwidth)
  loci.names <- unlist(lapply(loci, names))
  tmp <- grep("^loc", loci.names)
  loci.names[tmp] <- paste("c", names(unlist(loci))[tmp], sep = "")

  ## First create matrix. Want to preserve duplicate row names.
  n.traits <- length(traitnames)
  map <- matrix(NA, length(loci.names), 2 + n.traits)
  dimnames(map) <- list(loci.names, c("chr", "pos", traitnames))

  ## Now convert to data frame.
  map <- as.data.frame(map)
  
  ## And add chr, pos.
  map$chr <- ordered(rep(names(loci), sapply(loci, length)), names(loci))
  map$pos <- unlist(loci)

  map
}
########################################################################
read.maps <- function(cross, filename, chr.valid = names(cross$geno),
                   keep.nonseg = TRUE, drop.extra = TRUE,
                   reset.chr = TRUE, interp.loc = TRUE,
                   genotypes = c("A","H","B"),
                   verbose = TRUE, ...)
{
  ##*** cM.same and Mb.map do not plot properly.
  
  ## File should have snp, chr, loc, and orient columns (last is optional).
  geno <- read.table(filename, header = TRUE, fill = TRUE, comment.char = "")
  
  ## Add "rs" to SNP names.
  geno$snp <- paste("rs", geno$snp, sep = "")

  ## Pull map from cross object.
  cM.map <- pull.map(cross)
  cross.chr <- ordered(rep(names(cross$geno), nmar(cross)), names(cross$geno))

  ## Match SNP names to geno SNPs.
  cross.map <- unlist(sapply(cM.map, names))
  match.snp <- match(cross.map, geno$snp)

  ## Set up list of SNPs to keep (for non-segregating SNPs below).
  if(length(keep.nonseg) == 1) {
    if(!is.logical(keep.nonseg))
      keep.nonseg <- TRUE
    keep.nonseg <- array(keep.nonseg, nrow(geno))
    names(keep.nonseg) <- geno$snp
  }
  else {
    if(is.logical(keep.nonseg))
      tmp <- names(keep.nonseg)
    else {
      tmp <- as.character(keep.nonseg)
      keep.nonseg <- rep(TRUE, length(tmp))
      names(keep.nonseg) <- tmp
    }
    if(is.null(tmp))
      stop("keep option needs names to verify SNP order\n")
    tmp <- match(as.character(geno$snp), as.character(tmp))
    if(any(is.na(tmp)))
      stop("keep option names do not match file\n")
    keep.nonseg <- keep.nonseg[tmp]
  }

  ## Determine valid chrs. Damage control as needed for map.
  geno$chr <- ordered(geno$chr, chr.valid)
  tmp <- is.na(geno$chr)
  if(any(tmp)) {
    if(verbose) {
      cat(paste("Found", sum(tmp), "missing chr value for mapped marker(s):\n",
                paste(geno$snp[tmp], collapse = ","), "\n"))
    }
    tmp2 <- which(tmp[match.snp])
    if(reset.chr) {
      if(length(tmp2)) {
        tmp <- cross.chr
        geno$chr[match.snp[tmp2]] <- tmp[tmp2]
      }
      if(verbose)
        cat(" Value(s) reset to match cross chr:",
            paste(tmp[tmp2], collapse = ","), "\n")
    }
    else {
      if(verbose)
        cat(" marker(s) dropped from map.\n")

      ## Drop SNP from cM.map.
      for(i in unique(cross.chr[tmp2])) {
        ii <- tmp2[i == cross.chr[tmp2]]
        cM.map[[i]] <- cM.map[[i]][-match(cross.map[ii], names(cM.map[[i]]))]
      }

      ## Drop SNP from cross.chr and cross.map.
      cross.map <- cross.map[-tmp2]
      cross.chr <- cross.chr[-tmp2]
      ## Re-calibrate geno, keep.nonseg and match.snp.
      keep.nonseg <- keep.nonseg[-match.snp[tmp2]]
      geno <- geno[-match.snp[tmp2], ]
      match.snp <- match(cross.map, geno$snp)
    }
  }

  ## Get SNP names and match to geno.
  tmp <- is.na(match.snp)
  if(any(tmp)) {
    ## Find markers that don't match.
    if(verbose) {
      cat("Found", sum(tmp), "extra marker(s) not in file:\n",
          paste(cross.map[tmp], collapse = ","), "\n")
      cat(" chr:", paste(cross.chr[tmp], collapse = ","), "\n")
    }
    
    if(drop.extra) {
      ## Drop extra markers.
      if(verbose)
        cat(" Extra marker(s) dropped\n")
      match.snp <- match.snp[!tmp]
      tmp2 <- split(tmp, cross.chr)
      cross.map <- cross.map[!tmp]
      cross.chr <- cross.chr[!tmp]
      for(i in names(cM.map)) {
        if(any(tmp2[[i]])) {
          if(all(tmp2[[i]]))
            stop(paste("all SNPs missing for chr", i))
          cM.map[[i]] <- cM.map[[i]][!tmp2[[i]]]
        }
      }
    }
    else {
      if(verbose)
        cat(" Extra marker(s) interpolated in Mb.\n")
      
      ## Add extra markers to geno data frame.
      extra.markers <- which(tmp)
      geno <- rbind(geno,
                    data.frame(snp = cross.map[extra.markers],
                               chr = cross.chr[extra.markers],
                               loc = rep(NA, length(extra.markers)),
                               orient = rep(NA, length(extra.markers))))
      match.snp <- match(cross.map, geno$snp)
    }
  }

  ## Verify that chrs have not changed.
  tmp <- geno$chr[match.snp] == cross.chr
  if(!all(tmp)) {
    ## Some chr have changed--need to rethink map.
    d <- geno[match.snp,][!tmp,, drop = FALSE]
    d$cross.chr <- cross.chr[!tmp]
    d$cross.loc <- unlist(cM.map)[!tmp]
    tmp <- paste(nrow(d), "SNP moved to different chr")
    cat("\nWARNING:", tmp, "\n")
    print(d)
    cat("\n")
    stop(paste(tmp, "(see table above)\nDo you need to rebuild genetic map?"))
  }

  ## Need to take care of NA for geno$loc.
  ## Create Mb.map. Interpolate missing locations.
  tmp <- geno$loc[match.snp] / 1e6
  tmp2 <- is.na(tmp)
  if(any(tmp2)) {
    if(verbose) {
      cat(paste("Found", sum(tmp2), "missing loc value(s) for mapped SNP:\n",
                paste(geno$snp[match.snp][tmp2], collapse = ","), "\n"))
      cat(" chr:", paste(geno$chr[match.snp][tmp2], collapse = ","), "\n")
    }
  }

  Mb.map <- split(tmp, geno$chr[match.snp])
  class(Mb.map) <- "map"
  for(i in names(cM.map)) {
    class(Mb.map[[i]]) <- class(cM.map[[i]])
    names(Mb.map[[i]]) <- names(cM.map[[i]])
  }
  if(interp.loc) {
    if(verbose)
      cat(" Value(s) interpolated in Mb using cross map and other markers\n")
    for(i in names(cM.map)) {
      ## Reset missing values here.
      tmp <- is.na(Mb.map[[i]])
      if(any(tmp)) {
        if(all(tmp))
          stop(paste("All loc missing for chr", i))
        
        Mb.map[[i]][tmp] <- myapprox(cM.map[[i]][!tmp], Mb.map[[i]][!tmp],
                                     cM.map[[i]][tmp])$y
      }
    }
  }
  else {
    if(verbose)
      cat(" SNP dropped from map.\n")

    ## Drop SNP from cM.map.
    for(i in names(cM.map)) {
      ## Reset missing values here.
      tmp <- is.na(Mb.map[[i]])
      if(any(tmp)) {
        if(all(tmp))
          stop(paste("All loc missing for chr", i))
        
        Mb.map[[i]] <- Mb.map[[i]][!tmp]
        cM.map[[i]] <- cM.map[[i]][!tmp]
      }
    }

    ## Drop SNP from cross.map.
    cross.map <- unlist(sapply(cM.map, names))
    match.snp <- match(cross.map, geno$snp)

    ## Re-calibrate geno, keep.nonseg and match.snp.
    keep.nonseg <- keep.nonseg[-match.snp[tmp2]]
    geno <- geno[-match.snp[tmp2], ]
    match.snp <- match(cross.map, geno$snp)
  }

  ## Verify if any SNPs change order.
  is.amiss <- sapply(Mb.map, function(x) any(diff(x) < 0))
  if(any(is.amiss)) {
    is.amiss <- names(is.amiss[is.amiss])
    warning(paste("Marker order changed on chr:", paste(is.amiss, collapse = ",")))
  }
  
  ## NON-SEGREGATING SNPs (not on genetic map).
  
  ## Drop any non-segregating SNPs with missing chr or loc.
  keep.nonseg[is.na(geno$loc[keep.nonseg]) | is.na(geno$chr[keep.nonseg])] <- FALSE

  ## Now get SNPs that do not map (non-segregating).
  keep.nonseg[match.snp] <- FALSE
  if(any(keep.nonseg)) {
    Mb.same <- split(geno$loc[keep.nonseg] / 1e6, geno$chr[keep.nonseg])
    tmp <- split(geno$snp[keep.nonseg], geno$chr[keep.nonseg])
    for(i in names(Mb.same)) {
      tmp2 <- Mb.same[[i]]
      names(tmp2) <- tmp[[i]]
      Mb.same[[i]] <- sort(tmp2)
    }
    class(Mb.same) <- "map"
    for(i in names(Mb.same))
      class(Mb.same[[i]]) <- class(cM.map[[i]])
    
    ## Create cM.same.
    cM.same <- Mb.same
    for(i in names(cM.same))
      cM.same[[i]] <- myapprox(Mb.map[[i]], cM.map[[i]], cM.same[[i]])$y
  }
  else {
    cM.same <- Mb.same <- NULL
  }

  out <- list(cM.map = cM.map, Mb.map = Mb.map,
              cM.same = cM.same, Mb.same = Mb.same)
  attr(out, "genotypes") <- genotypes
  class(out) <- c("read.maps", "list")

  out
}
################################################################
print.read.maps <- function(x, ...) print(summary(x, ...))
################################################################
summary.read.maps <- function(object, ...)
{
  lapply(object, summary, ...)
}
################################################################
plot.read.maps <- function(x, main = "genetic vs. physical maps", ...)
{
  is.same <- !is.null(x$cM.same)
  if(is.same) {
    tmpar <- par(mfrow = c(2,1))
    on.exit(par(tmpar))
  }
  plot.map(x$cM.map, x$Mb.map, main = main, ...)
  if(is.same) {
    title("\n\nsegregating markers")
    plot.map(x$cM.same, x$Mb.same, main = main, ...)
    title("\n\nnon-segregating markers")
  }
}
