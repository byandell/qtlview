#######################################################################
calc.hc <- function(scans, cluster)
{
  if(cluster & ncol(scans) >= 4) {
    ## Sort by hierarchical cluster ordering.
    ## Can only cluster by LOD scores as traits may not be available.
    tmpfn <- function(x) if(all(x == 0)) x else (x / max(x, na.rm = TRUE))

    rlod <- apply(as.matrix(scans[, -(1:2)]), 2, tmpfn)
    hc <- hclust(dist(t(rlod)))
  }
  else
    hc <- NULL
  
}
#######################################################################
get.scans <- function(traitnames = NULL,
                      category="clinical",
                      filename=NULL,
                      cross,
                      maps,
                      chr = names(cross$geno),
                      threshold.lod = 0,
                      sex = "both",
                      method = "hk",
                      scan.name = paste(sex, method, sep = "."),
                      no.X = FALSE,
                      sep = ",",
                      trait.annotation = NULL,
                      drop.chr = TRUE,
                      ...)
{
  ## This assumes scans are saved in CSV with traits = rows and markers = columns.
  ## This format may not work for more than a few thousand markers...
  
  ## The category and tissue.name are synonomous.
  category <- tolower(category)
  tissue.name <- category

  ## Read CSV file with traits (rows) by markers (columns).
  lods <- read.table(filename, header = TRUE, sep = sep)
  
  ## Check if first (really any) column is called Tissue.
  ## If it is, extract it as tissue vector and remove it from data frame.
  tmp <- match("Tissue", names(lods))
  if(!is.na(tmp)) {
    tissues <- lods[[tmp]]
    lods <- lods[, -tmp]
  }
  else
    tissues <- NULL
  
  ## First column should now have trait names.
  lod.names <- as.character(lods[[1]])
  
  ## Transpose remaining columns as new lods.
  lods <- t(lods[,-1])
  
  ## Add tissue to lod names.
  if(!is.null(tissues))
    lod.names <- paste(tissues, lod.names, sep = ".")

  ## Put lod names as new column headers.
  ## Row names are marker names.
  dimnames(lods)[[2]] <- lod.names
  
  ## Get chr and pos and make room for traits.
  scans <- pull.pseudomarkers(cross, traitnames = lod.names, ...)
  
  ## Here expand lods using row.names to fill out missing markers.
  ## Put lods in to scans object.
  ## Have to do it this funny way as row.names may be duplicated.
  tmp <- match(make.names(row.names(lods)),
               make.names(row.names(scans)), nomatch = 0)
  if(any(tmp == 0)) {
    warning(paste(sum(tmp == 0), "markers do not match cross object"))
    if(all(tmp == 0))
      stop("no markers match cross object")
  }
  scans <- scans[tmp, ]
  scans[, lod.names] <- lods[tmp > 0, lod.names]
  tmp <- table(scans$chr)
  scans$chr <- ordered(scans$chr, names(tmp)[tmp > 0])
  rm(lods, lod.names)
  gc()
  
  attr(scans, "method") <- scan.name

  ## Reduce to desired trait names.
  if(!length(traitnames))
    traitnames <- make.names(names(scans)[-(1:2)])
  traitnames <-
    names(scans)[match(make.names(traitnames, allow_ = FALSE),
                       make.names(names(scans), allow_ = FALSE),
                       nomatch = 0)]

  ## Replace missing values with zeroes.
  scans[is.na(scans)] <- 0

  scan.attr <- attributes(scans)

  scans <- scans[, c(TRUE, TRUE,
                     make.names(names(scans)[-(1:2)], allow_ = FALSE)
                     %in% make.names(traitnames, allow_ = FALSE))]
  for(i in names(scan.attr))
    if(is.null(attr(scans, i)))
      attr(scans, i) <- scan.attr[[i]]

  attr(scans, "tissue.name") <- tissue.name
  attr(scans, "category") <- category
  
  tmp <- find.threshold(scans, chr, threshold.lod)

  if(drop.chr) {
    ## Reduce to chr with max(lod) > threshold.lod.
    scans <- scans[scans$chr %in% tmp$chr, ]
    if(!nrow(scans)) {
      tmp <- paste("\n\n*** No chromosomes pass threshold",
                   signif(threshold.lod, 5), "for any trait ***\n\n")
      cat(tmp)
      stop(tmp)
    }
    scans$chr <- ordered(scans$chr, tmp$chr)
  }
  
  ## Identify traits with max(lod) > threshold.lod.
  is.selected <- tmp$traits
  if(!any(is.selected)) {
    tmp <- paste("\n\n*** No trait passes threshold",
                 signif(threshold.lod, 5), "on any chromosome ***\n\n")
    cat(tmp)
    stop(tmp)
  }
  attr(scans, "is.selected") <- is.selected
  attr(scans, "maps") <- maps

  ## Add trait positions.
  trait.position <-
    find.trait.position(traitnames, maps, trait.annotation, ...)
  attr(scans, "trait.position") <- trait.position
    
  class(scans) <- c("aug.scanone", "scanone", "data.frame")
  
  scans
}
################################################################
composite.lod <- function(scans, is.selected, main, threshold.lod, support.lod,
                          scan.type = "LOD")
{
  ## Total over scans.
  scan.attr <- attributes(scans)
  scans <- scans[, c(TRUE, TRUE, is.selected)]
  
  tmpfn <- function(x) {
    x[x < 0] <- 0
    x[is.na(x)] <- 0
    sum(x)
  }
  if(scan.type == "BIM") {
    tmpfn2 <- function(x, support.lod)
      x / max(x, na.rm = TRUE)
###          x >= max(x, na.rm = TRUE) * (1 - min(1, support.lod))
  }
  else {
    tmpfn2 <- function(x, support.lod)
      x >= (max(x, na.rm = TRUE) - support.lod)
  }

  newdata <- scans[, 1:2]
#  names(newdata)[2] <- "pos.cM"
  newdata$"Composite LOD" <- rep(NA, nrow(newdata))
  newdata$"Peak LOD" <- rep(NA, nrow(newdata))
  for(i in levels(scans$chr)) {
    ii <- scans$chr == i
    if(any(ii)) {
      tmp <- apply(scans[ii, -(1:2), drop = FALSE], 2,
                   function(x) any(x > threshold.lod))
      if(any(tmp)) {
        tmp2 <- scans[ii, c(rep(FALSE, 2), tmp), drop = FALSE]
        
        tmp3 <- apply(tmp2, 2, tmpfn2, 0)
        ## Indicate support interval if using support.lod.
        if(!is.null(support.lod))
          tmp2 <- apply(tmp2, 2, tmpfn2, support.lod)
        ## Otherwise sum up LODs directly.
        
        newdata[[3]][ii] <- apply(tmp2, 1, tmpfn)
        newdata[[4]][ii] <- apply(tmp3, 1, tmpfn)
      }
      else {
        newdata[[3]][ii] <- newdata[[4]][ii] <- rep(0, sum(ii))
      }
    }
  }

  for(i in names(scan.attr))
    if(is.null(attr(newdata, i)))
       attr(newdata, i) <- scan.attr[[i]]

  newdata
}
################################################################
sumone.scanone <- function(object, chr = unique(sum$chr), ...)
{
  sum <- summary(object, ...)
  sum2 <- duplicated(sum$chr)
  if(any(sum2)) {
    ## Multiple ties to max lod.
    sum2 <- sum[!sum2, ]
    sum2$pos <- unlist(tapply(sum$pos, sum$chr, mean))[sum2$chr]
    sum <- sum2
  }
  sum[match(chr, sum$chr), ]
}
################################################################
best2pattern <- function(best)
{
  ## Infer best pattern for use in make.bestscan.
  ## This misses epistasis and sets all posteriors to 1.

  ## Pattern is comma-separated chromosome names.
  tmp <- ordered(best[[1]], unique(best[[1]]))
  mypat <- tapply(best$chrom, tmp,
                  function(x) paste(x, collapse = ","))
  mypat[mypat == "NA"] <- "NULL"

  ## Set up data frame.
  pattern <- data.frame(x = levels(tmp))
  names(pattern) <- names(best)[1]
  pattern$pattern <- ordered(mypat)
  pattern$posterior <- rep(1, nrow(pattern))
  
  pattern
}
################################################################
make.bestscan <- function(cross, best, pattern = best2pattern(best),
                          scan.type = "BIM",
                          cross.name = deparse(substitute(cross)),
                          maps = myget(cross.name, "maps"),
                          trait.annotation = NULL,
                          attenuate = 1, tolerance = 1e-5,
                          drop.null = TRUE, verbose = TRUE)
{
  ## Want to infer best if not provided.
  
  loci <- pull.pseudomarkers(cross)
  loci.attr <- attributes(loci)

  ## Assume chr split denoted by period (.) in chr name (discarded).
  chr <- sapply(strsplit(as.character(best$chrom), ".", fixed = TRUE),
                function(x) x[1])
  ncol.best <- ncol(best)

  traitnames <- as.character(pattern[[1]])
  if(drop.null)
    traitnames <- traitnames[pattern$pattern != "NULL"]

  n.traits <- length(traitnames)

  ## Create aug.scanone object of attenuated pattern scans.
  if(verbose)
    cat("Creating", n.traits, "x", nrow(loci), "aug.scanone object of pattern scans...\n")
  bestscan <- matrix(0, n.traits, nrow(loci))
  dimnames(bestscan) <- list(traitnames, row.names(loci))
      
  for(i in 1:19) {
    if(verbose) cat("chr", i, "\n")
    
    locus.i <- which(loci$chr == i)
    mcmc.i <- which((chr == i) & !is.na(chr))
    trait.match <- match(as.character(best[mcmc.i, 1]),
                         traitnames, nomatch = 0)
    
    ## Create length(locus.i) * length(mcmc.i) matrix of (1-2r)^2 = exp(-2*2*cM/100).
    atten.locus <- t(apply(as.matrix(mcmc.i), 1,
                           function(x, y) {
                             out <- best[x, ncol.best] *
                               exp(- 0.04 * attenuate * abs(best[x, "locus"] - y))
                             ## Next step saves space.
                             out[out <= tolerance] <- 0
                             out
                           },
                           loci[locus.i, "pos"]))
    bestscan[trait.match, locus.i] <- atten.locus
  }
  bestscan
}
