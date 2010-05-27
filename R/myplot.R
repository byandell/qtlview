################################################################
myplot.scanone <- function(scans,
                           threshold.lod = 0, main = "",
                           col = c("black","blue","red","green",
                             "purple","magenta"),
                           ylim = ylims, add.position = FALSE,
                           lod.name = "LOD",
                           use.cM = TRUE,
                           chr = NULL, xaxs = "r",
                           ...)
{
  chr <- find.threshold(scans, chr, threshold.lod)$chr

  maps <- attr(scans, "maps")

  traitnames <- names(scans)[-(1:2)]

  n.traits <- length(traitnames)

  ## This is code duplication of myplot(). Fix later.
  class(scans) <- c("scanone", "data.frame")

  col <- array(col, n.traits)
  names(col) <- traitnames

  ## Set vertical limits to include all peaks.
  ylims <- range(unlist(scans[scans$chr %in% chr, -(1:2)]))
  if(ylims[1] < 0)
    ylims[1] <- 0
  
  ## Set up cM to Mb translation.  
  if(length(chr) == 1) {
    p <- qm.approx(maps, ifelse(use.cM, "Mb", "cM"), chr)
    if(!use.cM)
      scans$pos <- qm.approx(maps, "cM", chr, pos = scans$pos)$y
  }

  names(scans)[3] <- lod.name
  for(i in seq(traitnames)) {
    tmp <- tapply(scans$pos,scans$chr,function(x)diff(range(x)))
    tmp <- tmp[!is.na(tmp)]
    par(xaxs = xaxs)
    plot(scans, lodcolumn = i, add = (i > 1), chr = chr,
         col = col[i], ylim = ylim)

    if(i == 1)
      add.rug(chr, main, maps, p, use.cM)
  }
  ## Need to modify this for separate A and X thresholds.
  threshold.lines(scans, threshold.lod, ...)
  
  invisible(col)
}
################################################################
mysum.scanone <- function(scans, threshold.lod, lod.name = "LOD",
                          col = c("black","blue","red","green","purple","magenta"),
                          ...)
{
  traitnames <- names(scans)[-(1:2)]
  n.traits <- length(traitnames)

  col <- array(col, n.traits)
  names(col) <- traitnames

  sum <- summary.aug.scanone(scans, traitnames = traitnames,
                 threshold.lod = threshold.lod, long = TRUE, ...)

  if(nrow(sum$pos)) {
    if(n.traits == 1) {
      cat(traitnames, "\n")
      sum <- data.frame(pos = round(sum$pos[1,], 3), lod = sum$lod[1,])
      names(sum) <- c("pos", lod.name)
    }
    else {
      maps <- attr(scans, "maps")
      sum <- maxit(sum, maps, attr(scans, "trait.position"))
      sum$col = factor(col[row.names(sum)])
      sum <- sum[apply(sum, 1, function(x) !all(is.na(x))), ]
    }
  }
  sum
}
################################################################
maxit <- function(sum, maps, trait.position)
{
  wh <- apply(sum$lod, 1, which.max)
  traitnames <- dimnames(sum$lod)[[1]]
  n.trait <- length(traitnames)
  out <- data.frame(chr = factor(dimnames(sum$lod)[[2]][wh]),
                    pos = round(sum$pos[seq(n.trait) +
                      (wh - 1) * n.trait], 3),
                    lod = sum$lod[seq(n.trait) +
                      (wh - 1) * n.trait])
  row.names(out) <- traitnames

  ## Add trait positions in cM and Mb.
  add.trait.position(out, trait.position)
}
