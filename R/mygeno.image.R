mygeno.image <- function(x = get(cross.name),
                         chr = names(x$geno),
                         cross.name = deparse(substitute(x)),
                         maps = myget(cross.name, "maps"),
                         reorder = TRUE, id = "MouseNum",
                         sex = c("both","male","female"),
                         genotypes = attr(maps, "genotypes"),
                         xlim = range(unlist(map)),
                         use.cM = FALSE,
                         reorder.by.genotype = any(reorder > 0),
                         recomb.only = FALSE,
                         keep.missing = TRUE,
                         ...)
{
  ## Problem if maps do not match genotypes in cross.
  ## Need to check?
  sex <- match.arg(sex)
  
  if(!inherits(x, "cross"))
    stop("Input should have class \"cross\".")
  
  cross <- subset(x, chr = chr, ...)
  if(sex != "both") {
    sexpgm <- getsex(cross)
    cross <- subset(cross, ind = (sexpgm$sex == (sex == "male")))
  }

#  if(nchr(cross) > 1)
#    stop("Genotype image for single chromosome only.")

  chr <- names(cross$geno)

  ## Get ids.
  ids <- cross$pheno[[id]]
  
  ## Reorder individuals based on phenotype?
  ## (borrowed from [qtl]{plot.missing})
  Geno <- pull.geno(cross)

  map <- if(use.cM) maps$cM.map else maps$Mb.map
  if(length(chr) == 1) {
    map <- map[[chr]]
    xlim <- xlim.range(xlim, map)
    in.range <- map >= xlim[1] & map <= xlim[2]
    mean.geno <- apply(Geno[, in.range, drop = FALSE], 1,
                       mean, na.rm = TRUE)
  }
  else {
    map <- map[chr]
    in.range <- rep(TRUE, ncol(Geno))
    mean.geno <- apply(Geno, 1, mean, na.rm = TRUE)
  }
  
  ordset <- NULL
  if(is.character(reorder))
    reorder <- find.pheno(cross, reorder)

  if(is.numeric(reorder)) {
    if (reorder < 1 || reorder > nphe(cross)) 
      stop("reorder should be an integer between 1 and", nphe(cross))
    ordset <- as.matrix(cross$pheno[, reorder, drop = FALSE])
  }
  else
    ordset <- NULL
  
  if(reorder.by.genotype) {
    ## Reorder by genotype.
    o <- order(mean.geno)
    
    Geno <- Geno[o, ]
    ids <- ids[o]
    mean.geno <- mean.geno[o]
    if(!is.null(ordset))
      ordset <- ordset[o,, drop = FALSE]
  }

  ## Drop all but recombinants.
  if(recomb.only) {
    ## This seems to be keeping those that have all missing values. Why?
    keep <- apply(Geno[, in.range, drop = FALSE], 1,
                  function(x) !all(x == mean(x, na.rm = TRUE)))
    keep[is.na(keep)] <- FALSE
    Geno <- Geno[keep,, drop = FALSE]
    ids <- ids[keep]
    mean.geno <- mean.geno[keep]
    if(!is.null(ordset))
      ordset <- ordset[keep,, drop = FALSE]
  }

  ## Drop missing?
  if(!keep.missing) {
    keep <- apply(Geno[, in.range, drop = FALSE], 1,
                  function(x) !any(is.na(x)))
    Geno <- Geno[keep,, drop = FALSE]
    ids <- ids[keep]
    mean.geno <- mean.geno[keep]
    if(!is.null(ordset))
      ordset <- ordset[keep,, drop = FALSE]
  }
  
  class(Geno) <- c("mygeno.image", "matrix")
  attr(Geno, "chr") <- chr
  attr(Geno, "id") <- ids
  attr(Geno, "geno") <- mean.geno
  attr(Geno, "ordset") <- ordset
  if(length(chr) == 1)
    attr(Geno, "xlim") <- xlim
  attr(Geno, "recomb.only") <- recomb.only
  attr(Geno, "keep.missing") <- keep.missing
  attr(Geno, "maps") <- maps

  ## Set genotype 4 to "same" for non-segregating markers.
  if(is.null(genotypes))
    genotypes <- c("A","H","B")
  genotypes <- genotypes[1:4]
  genotypes[4] <- "same"
  attr(Geno, "genotypes") <- genotypes

  col <- c("blue","green","red","lightgray")
  names(col) <- c(1:3,"NA")
  attr(Geno, "col") <- col
  
  Geno
}
##################################################################################
summary.mygeno.image <- function(object,
                                 ...)
{
      
  cat("Genotypes are coded 1, 2, 3; geno column is mean over markers\n\n")

  genotypes <- attr(object, "genotypes")
  col <- attr(object, "col")
  cat(paste(col, "=", names(col), "=", genotypes,
            collapse = ", "), "\n\n")

  cat(paste("recomb.only =", attr(object, "recomb.only"),
            "keep.missing =", attr(object, "keep.missing"),
            "\n\n"))

  ## Set up data frame with id, geno, any phenos.
  tmp <- data.frame(id = attr(object, "id"), geno = attr(object, "geno"))
  tmp2 <- seq(nrow(tmp))
  row.names(tmp) <- tmp2
  if(!is.null(attr(object, "ordset")))
    tmp <- cbind(tmp, attr(object, "ordset"))

  tmp[rev(tmp2), ]
}
print.mygeno.image <- function(x, ...) print(summary(x, ...))
##################################################################################
xlim.range <- function(xlim, map)
{
  ## Set up x limits if not provided or NULL.
  tmp <- range(map)
  if(is.null(xlim))
    xlim <- tmp
  if(length(xlim) != 2)
    xlim <- tmp
  else {
    xlim[1] <- max(xlim[1], tmp[1])
    xlim[2] <- min(xlim[2], tmp[2])
  }
  xlim
}
##################################################################################
plot.mygeno.image <- function(x,
                              xlim = attr(x, "xlim"),
                              chr = achr,
                              use.cM = FALSE,
                              equal.spacing = FALSE,
                              zscale = TRUE,
                              normal.score = TRUE,
                              main = "", ...)
{
  genotypes <- attr(x, "genotypes")
  genotypes <- genotypes[-length(genotypes)]
  geno <- attr(x, "geno")
  ordset <- attr(x, "ordset")

  maps <- attr(x, "maps")
  if(use.cM) {
    map <- maps$cM.map
    non.seg <- maps$cM.same
  }
  else {
    map <- maps$Mb.map
    non.seg <- maps$Mb.same
  }

  col <- attr(x, "col")
  breaks <- 0.5 + seq(0,length(col))

  ## Check on selected chr.
  achr <- attr(x, "chr")
  chr <- as.character(chr)
  if(!all(chr %in% achr))
    stop(paste("Chromosomes must be in selected set:",
               paste(achr, collapse = ",")))
  
  if(length(achr) > 1 & length(chr) < length(achr)) {
    ## Need to reduce x to selected chr.
    tmp <- sapply(map[achr], length)
    tmp <- apply(as.matrix(tmp), 2, function(x,y) rep(y, x),
                 achr %in% chr)
    x <- x[, tmp]
  }

  tmpar1 <- par(mar = c(4.1,4.1,3.1,0.1))
  on.exit(par(tmpar1))
  ## Make room for Z scale.
  dots <- list(...)
  if (zscale) {
    if ("layout" %in% names(dots)) 
      layout(dots[["layout"]][[1]], dots[["layout"]][[2]])
    else layout(cbind(1, 2), c(10, 1))
    on.exit(layout(1,1))
}

  if(length(chr) == 1) {
    ## Get map. Stitch in non-segregating markers if provided.
    map <- c(map[[chr]])
    non.seg <- c(non.seg[[chr]])

    xlim <- xlim.range(xlim, map)
  
    if(equal.spacing) {
      pos <- seq(length(map))
      xlim <- range(pos[map >= xlim[1] & map <= xlim[2]]) + c(-0.5, 0.5)

      imago <- t(x)
    }
    else {
      pos <- unique(c(map, non.seg))
      o <- order(pos)
      pos <- pos[o]
      
      ## Need to drop any non.seg that have identical positions.
      imago <- matrix(4, length(o), nrow(x))
      imago[o %in% seq(length(map)), ] <- t(x)
    }

    ## Plot image of chr.
    tmpar <- par(xaxt = "n")
    image(pos, seq(ncol(imago)), imago, col = col, breaks = breaks,
          xlab = "", ylab = "mouse", xlim = xlim, ...)
    par(tmpar)
    if(equal.spacing) {
      axis(1, pos)
      mtext(paste("Chromosome", chr), 1, 2)
    }
    else
      add.rug(chr, "", maps, use.cM = use.cM, outer = TRUE,
              xlim = xlim, bottom.axis = TRUE)
  }
  else { ## Multiple chromosomes.
    map <- map[chr]
    n.chr <- sapply(map, function(x) length(c(x)))
    imago <- matrix(4, ncol(x) + length(chr) + 1, nrow(x))
    imago[-cumsum(1 + c(0, n.chr)), ] <- t(x)
    pos <- seq(nrow(imago))

    ## Plot image of chr.
    tmpar <- par(xaxt = "n")
    image(pos, seq(ncol(imago)), imago, col = col, breaks = breaks,
          xlab = "", ylab = "mouse", ...)
    par(tmpar)

    ## Add chr name.
    tmp <- ordered(c("0", rep(chr, 1 + n.chr)), c("0", chr))
    tmp <- tapply(pos, tmp, mean)[-1] + 0.5
    axis(1, tmp, labels = chr)
  }
  if(main != "")
    title(main)

  if(!is.null(ordset)) {
    ord.names <- dimnames(ordset)[[2]]
    n.trait <- ncol(ordset)
  }
  
  if (zscale) {
    par(mar = c(4.1,0.6,3.1,0.1))

    imageno <- t(as.matrix(geno))
    dimnames(imageno)[[1]] <- "geno"
    imageno <- imageno - 1
    imageno <- imageno / max(imageno, na.rm = TRUE)
    
    if(!is.null(ordset)) {
      tmpfn <- function(x, geno) {
        fit <- loess(x~geno)
        x[!is.na(x) & !is.na(geno)] <- fitted(fit)
        (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
      }
      sordset <- apply(ordset, 2, tmpfn, geno)
      tmpfn <- function(x, geno) {
        x <- rank(x, na.last = "keep")
        (x - min(x, na.rm = TRUE)) / diff(range(x, na.rm = TRUE))
      }
      tmp <- apply(ordset, 2, tmpfn, geno)
      n.o <- ncol(tmp)
      sordset <- cbind(sordset, tmp)
      sordset[, c(2 * seq(n.o) - 1, 2 * seq(n.o))] <- sordset
      imageno <- rbind(imageno, t(sordset))
    }

    col.scheme <- "gray"
    cols <- switch(col.scheme,
                   gray = gray(seq(0, 0.9, len = 256)),
                   heat = heat.colors(256),
                   terrain = terrain.colors(256), 
                   topo = topo.colors(256), cm = cm.colors(256),
                   redblue = rev(rainbow(256, 
                     start = 0, end = 2/3)))
    
    ## Add Scale for mean geno on col scheme.
    image(seq(nrow(imageno)), seq(ncol(imageno)), imageno,
          ylab = "", xlab = "",
          col = cols, yaxt = "n", xaxt = "n")
    tmp <- par("usr")
    abline(h = tmp[3:4], v = tmp[1:2])
    if(!is.null(ordset)) {
      tmp <- abbreviate(c("geno", ord.names), 10)
      for(i in seq(1 + n.trait))
        mtext(tmp[i], 3, 0, at = 2 * i - 1, las = 3)
    }
  }
  if(!is.null(ordset)) {
    par(mfrow = c(1, n.trait), mar = c(3.1,3.1,2.1,0.1))
    o <- order(geno)
    geno <- geno[o]
    tmp <- geno - 1
    tmp <- 1 + round(255 * tmp / max(tmp, na.rm = TRUE))
    cols <- rev(rainbow(256, start = 0, end = 2/3))[tmp]
    col <- col[-length(col)]
    for(i in seq(n.trait)) {
      x <- ordset[o, i]
      if(normal.score)
        x <- normal.trans(x)
      plot(x ~ geno, xlab = "", ylab = "", col = cols)
      mtext(paste(seq(genotypes), col, genotypes,
                  sep = "=", collapse = ", "),
            1, 2)
      mtext(ord.names[i], 2, 2)
      x[!is.na(x) & !is.na(geno)] <- fitted(loess(x ~ geno))
      lines(x ~ geno, lwd = 4)
      if(main != "")
        title(main)
    }
  }
}
