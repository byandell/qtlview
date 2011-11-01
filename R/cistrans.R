###############################################################################
cistrans <- function(x = read.table(filename, header = TRUE, sep = sep),
                     filename,
                     cross.name,
                     maps = myget(cross.name, "maps"),
                     peak.chr = peak.chrs,
                     trait.chr = trait.chrs,
                     min.lod = 0,
                     main = paste(levels(x$Tissue), collapse = ", "),
                     sep = "\t",
                     use.density = TRUE,
                     trait.annotation = myget(cross.name, "annotation"),
                     use.annot = FALSE,
                     ...)
{
  ## Kludge until PHP code changes:
  dots <- list(...)
  if ("trans.chr" %in% names(dots))
    trait.chr <- trans.chr

  use.annot <- use.annot & !is.null(trait.annotation)

  ## Translate x names from maxlod to peak if needed.
  maxlods <- c("maxlod","maxlod.chr","maxlod.pos.Mb","maxlod.pos.cM",
               "trans.chr","trans.pos.Mb","trans.pos.cM")
  peaks <- c("peak.score","peak.chr","peak.pos.Mb","peak.pos.cM",
               "trait.chr","trait.pos.Mb","trait.pos.cM")
  tmp <- match(maxlods, names(x), nomatch = 0)
  names(x)[tmp] <- peaks[tmp > 0]

  is.trait <- "trait.chr" %in% names(x)
  if(is.trait)
    is.trait <- !all(is.na(x$trait.chr))
  
  getreal <- function(x, realx = c(1:19,"X")) {
    if(is.null(x))
      NULL
    else {
      x <- levels(as.ordered(x))
      realx[realx %in% x]
    }
  }

  ## Set up peak chromosomes.
  peak.chrs <- getreal(x$peak.chr)
  x$peak.chr <- ordered(x$peak.chr, peak.chrs)
  chr.names <- names(maps$cM.map)
  peak.chr <- chr.names[chr.names %in% peak.chr]
  if(!length(peak.chr))
    stop("must select at least on peak chr")

  ## Set up trait chromosomes and locations.
  if(!is.trait | use.annot) {
    ## Use trait.annotation if none provided.
    tmp <- match(x$a.gene.id, trait.annotation$a_gene_id, nomatch = 0)
    x$trait.chr <- x$trait.pos.Mb <- rep(NA, nrow(x))
    x$trait.chr[tmp > 0] <- as.character(trait.annotation$Chromosome[tmp])
    x$trait.pos.Mb[tmp > 0] <- trait.annotation$Chromosome_Position[tmp]
    ## Be sure to reset cM position with new Mb positions.
    x$trait.pos.cM <- NULL
    is.trait <- !all(is.na(x$trait.chr))
  }
  if(is.trait) {
    trait.chrs <- getreal(x$trait.chr)
    x$trait.chr <- ordered(x$trait.chr, trait.chrs)
    trait.chr <- chr.names[chr.names %in% trait.chr]
    if(!length(trait.chr))
      stop("must select at least one trait chr")
  }

  ## Reduce to entries within range of peak.chr and trait.chr.
  tmp <- (x$peak.chr %in% peak.chr)
  if(is.trait)
    tmp <- tmp & (x$trait.chr %in% trait.chr)
  if(!any(tmp))
    stop("no transcripts selected")
  x <- x[tmp, ]

  ## Reduce to entries above min.lod.
  tmp <- x$peak.score >= min.lod
  if(!any(tmp))
    return(mystop(paste("\n\n*** No transcripts above minimum score of", min.lod, ". ***\n\n")))
  x <- x[tmp, ]

  ## Drop cis traits?
  drop.cis <- FALSE
  if(drop.cis & is.trait) {
    dist.cis <- 20
    tmp <- as.character(x$peak.chr) == as.character(x$trait.chr)
    tmp <- tmp & (abs(x$peak.pos.Mb - x$trait.pos.Mb) <= dist.cis)
    if(any(tmp))
      x <- x[!tmp, ]
  }

  ## Distance in cM and Mb. Create if not present.
  if(is.null(x$peak.pos.cM))
    x$peak.pos.cM <- mypos(x, maps, "Mb", "peak")
  if(is.null(x$peak.pos.Mb))
    x$peak.pos.Mb <- mypos(x, maps, "cM", "peak")
  if(is.null(x$trait.pos.cM) & is.trait)
    x$trait.pos.cM <- mypos(x, maps, "Mb", "trait")
  if(is.null(x$trait.pos.Mb) & is.trait)
    x$trait.pos.Mb <- mypos(x, maps, "cM", "trait")

  x <- cumscore(x, ...)

  class(x) <- c("cistrans", "data.frame") 

  if(use.density)
    attr(x, "dens") <- cumscore.dens(x, ...)

  attr(x, "main") <- main
  attr(x, "peak.chr") <- peak.chr
  attr(x, "trait.chr") <- trait.chr
  attr(x, "is.trait") <- is.trait
  attr(x, "maps") <- maps
  attr(x, "use.annot") <- use.annot

  x
}
#####################################################################
mypos <- function(out, maps, orig = "Mb", prefix = "peak")
{
  ## Recast position from cM to Mb or vice versa.
  ## This is very specific to cistrans objects. (used in pattern.R)
  
  pos <- paste(prefix, "pos", orig, sep = ".")
  chr <- paste(prefix, "chr", sep = ".")
  
  res <- rep(NA, nrow(out))
  tmp <- !is.na(out[[pos]]) & !is.na(out[[chr]])
  for(i in levels(out[[chr]])) {
    ii <- (i == out[[chr]]) & tmp
    if(any(ii))
      res[ii] <- qm.approx(maps, orig, i, out[[pos]][ii])$y
  }
  res
}
#####################################################################
mycol <- function(x,
                  n.col = 256,
                  col.scheme = c("redblue","cis","gray","redblue2"),
                  attenuate = (col.scheme != "gray"),
                  ...)
{
  ## I think all this works well. However, it does not include a fixed window
  ## for cis vs. trans, as biologists seem to want in practice.
  ## This could be accomplished by (a) adding a window size and
  ## (b) making att 0/1 if the window approach is used.
  ## Need to think how best to accommodate this.
  
  ## Want to be able to modify the red so that it reflects attenuation.
  col.scheme <- match.arg(col.scheme)
  if(col.scheme == "redblue") 
    cols <- rev(rainbow(n.col, start = 0, end = 2/3))
  else {
    if(col.scheme == "redblue2")
      cols <- rainbow(n.col, start = 2/3, end = 1)
    else
      cols <- rev(gray(seq(0, 0.8, len = n.col)))
  }
  
  att <- rep(0, nrow(x))
  attenuate <- attenuate & !is.null(x$trait.chr)
  if(attenuate) {
    tmp <- (as.character(x$peak.chr) == as.character(x$trait.chr) &
            !is.na(x$peak.chr) & !is.na(x$trait.chr))
    if(any(tmp))
      att[tmp] <- exp(- 0.04 * abs(x$peak.pos.Mb[tmp] - x$trait.pos.Mb[tmp]))
    gb <- rev(seq(0, ifelse(col.scheme == "cis", 0.8, 1), len = n.col))
  }
    
  ## Set up colors along rank(peak.score).
  tmp <- rank(x$peak.score, na.last = "keep")
  index <- ceiling(n.col * tmp / max(tmp, na.rm = TRUE))

  ## This is not quite right.
  if(attenuate) {
    if(col.scheme == "cis")
      val <- rgb(att + (1 - att) * gb[index], gb[index], gb[index])
    else { ## black at cis, blue-red outside using bluered2
      n.col2 <- round(n.col / 2)
      r.col2 <- 1 / n.col2
      if(col.scheme == "redblue2") {
        ## redblue2 (via magenta)
        red <- c(seq(0, by = r.col2, len = n.col2), rep(1, n.col - n.col2))
        green <- rep(0, n.col)
        blue <- c(rep(1, n.col2), seq(1 - r.col2, by = -r.col2, len = n.col - n.col2))
      }
      else {
        ## redblue (via green--almost)
        n.col4 <- round(n.col / 4)
        n.col3 <- round(3 * n.col / 4)
        r.col4 <- 1 / n.col4
        red <- c(rep(0, n.col2),
                 seq(r.col2, by = r.col4, len = n.col3 - n.col2),
                 rep(1, n.col - n.col3))
        green <- c(seq(0, by = r.col4, len = n.col4),
                   rep(1, n.col3 - n.col4),
                   rev(seq(0, by = r.col4, len = n.col - n.col3)))
        blue <- c(rep(1, n.col4),
                  rev(seq(0, by = r.col4, len = n.col2 - n.col4)),
                  rep(0, n.col - n.col2))
      }
      
      val <- rgb(((1 - att) * red[index] + att * blue[index]) * (1 - 0.2 * att),
                 ((1 - att) * green[index] + att * blue[index]) * (1 - 0.2 * att),
                 blue[index] * (1 - 0.2 * att))
    }
  }
  else
    val <- cols[index]
  
  list(color = cols, index = index, value = val, att = att)
}
#####################################################################
plot.cistrans <- function(x, cex = ifelse(n.peak == 1 | n.trans == 1, 1, 0.2),
                          main = attr(x, "main"),
                          col.score = c("score","cis","trans"),
                          cis.only = (col.score == "cis" & n.peak > 1 & n.trans > 1),
                          xlim = xlims,
                          ylim = ylims,
                          use.cM = FALSE,
                          jitter = FALSE,
                          peak.chr = attr(x, "peak.chr"),
                          trait.chr = attr(x, "trait.chr"),
                          show.cumscore = TRUE,
                          pch = 1,
                          ...)
{
  ## Add equal.spacing here and for cumscore? Is it worth it?
  ## Need to aline x$peak.pos.Mb with Mb.map. I have some code somewhere...

  n.peak <- length(peak.chr)
  n.trans <- length(trait.chr)
  is.trait <- attr(x, "is.trait")
  if(!is.trait)
    n.trans <- 0

  maps <- attr(x, "maps")

  col.score <- match.arg(col.score)
  if(col.score %in% c("cis","trans") & is.trait) {
    ## Attenuate peak score.
    tmp <- as.character(x$peak.chr) == as.character(x$trait.chr)
    att <- exp(- 0.04 * abs(x$peak.pos.cM[tmp] - x$trait.pos.cM[tmp]))
    if(col.score == "cis") {
      x$peak.score[tmp] <- x$peak.score[tmp] * att
      x$peak.score[!tmp] <- 0
    }
    else
      x$peak.score[tmp] <- x$peak.score[tmp] * (1 - att)
  }
    
  ## Would be nice to incorporate the same chrom (diag) plots
  ## of plot.cis.dist in pattern.R.
  if(cis.only) {
    print(plot.cis.dist(x, col.score = col.score, main = main,
                        use.cM = use.cM, peak.chr = peak.chr, ...))
  }
  else {

    ## Default gap from R/qtl scanone.
    gap <- 25
    
    ## Select segregating and non-segregating maps.
    if(use.cM) {
      seg <- maps$cM.map
      non.seg <- maps$cM.same
    }
    else {
      seg <- maps$Mb.map
      non.seg <- maps$Mb.same
    }
    
    ## Set up peak and trans lengths
    catmap <- function(map.len, gap) {
      map.len <- cumsum(c(0, gap + map.len))
    names(map.len) <- c(names(map.len)[-1], "")
      ## Added one too many gaps.
      map.len[length(map.len)] <- map.len[length(map.len)] - gap
      map.len
    }
    map.len <- sapply(seg, max)
    trans.len <- catmap(map.len, gap)
    n.map <- length(trans.len)
    
    tmp <- names(map.len)[-n.map]
    if(is.trait)
      trans.len <- catmap(map.len[trait.chr], gap)
    peak.len <- catmap(map.len[peak.chr], gap)
    
    ## Cumulative position for easy plotting.
    tmp <- ifelse(use.cM, "cM", "Mb")
    pos <- paste("peak.pos", tmp, sep = ".")
    x$peak.cumpos <- x[[pos]] + peak.len[as.character(x$peak.chr)]
    if(is.trait) {
      pos <- paste("trait.pos", tmp, sep = ".")
      x$trans.cumpos <- x[[pos]] + trans.len[as.character(x$trait.chr)]
    }
    else
      x$trans.cumpos <- x$peak.score
    
    ## Force xlim, ylim if either includes more than one chr.
    xlims <- range(x$peak.cumpos, na.rm = TRUE)
    if(n.peak > 1 | is.null(xlim))
      xlim <- xlims
    ylims <- range(x$trans.cumpos, na.rm = TRUE)
    if(!is.trait)
      ylims <- range(0, ylims)
    if(n.trans > 1 | is.null(ylim))
      ylim <- ylims
    else
      par(mar = c(5.1,4.1,4.1,3.1))
    
    ## Set up colors.
    if(col.score == "trans")
      col.score <- "cis"
    colors <- mycol(x, ...)
    index <- order(colors$index)
    colors <- colors$value

    plot(trans.cumpos ~ peak.cumpos, x, cex = cex, type = "n",
         xlab = "chromosome of peak score",
         ylab = ifelse(is.trait, "chromosome of transcript", "peak.score"),
         xlim = xlim, ylim = ylim,
         xaxt = "n", yaxt = ifelse(is.trait, "n", "s"),
         xaxs = ifelse(n.peak == 1, "r", "i"),
         yaxs = ifelse(n.trans <= 1, "r", "i"))
    if(n.peak == n.trans)
      abline(0, 1,  col = "gray")
    tmp <- if(jitter & n.peak == 1)
      formula(trans.cumpos ~ jitter(peak.cumpos))
    else
      formula(trans.cumpos ~ peak.cumpos)
    points(tmp, x[index,], col = colors[index], cex = cex, pch = pch)
    
    if(n.peak == 1)
      add.rug(peak.chr, main, maps, use.cM = use.cM, outer = TRUE,
              bottom.axis = TRUE, xlim = xlim)
    else {
      axis(1, ((peak.len[-(1 + n.peak)] + peak.len[-1]) / 2) - gap / 2, peak.chr)
      abline(v = peak.len[2:n.peak] - gap / 2)
      title(main)
    }
    if(is.trait) {
      if(n.trans == 1)
        add.rug(trait.chr, "", maps, use.cM = use.cM, outer = TRUE,
                bottom.axis = TRUE, xlim = ylim, side = 2)
      else {
        axis(2, ((trans.len[-(1 + n.trans)] + trans.len[-1]) / 2) - gap / 2,
             trait.chr)
        abline(h = trans.len[2:n.trans] - gap / 2)
      }
    }
  }

  if(show.cumscore)
    plot.cumscore(x, peak.chr = peak.chr, main = main, xlim = xlim,
                  use.cM = use.cM, maps = maps, ...)
  
  invisible()
}
#####################################################################
cumscore <- function(x, ...)
{
  x <- x[order(x$peak.chr, x$peak.pos.Mb), ]

  ## Cumulate score.
  x$peak.cumscore <- unlist(tapply(x$peak.score,
                                   x$peak.chr,
                                   function(x) {
                                     out <- cumsum(x)
                                     out / max(out)
                                   }))
  ## Cumulate count.
  x$peak.cumcount <- unlist(tapply(x$peak.score > 0,
                                   x$peak.chr,
                                   function(x) {
                                     out <- cumsum(x)
                                     out / max(out)
                                   }))

  is.trait <- !is.null(x$trait.chr)
  if(is.trait)
    is.trait <- !all(is.na(x$trait.chr))
  if(is.trait) {
    tmp <- (as.character(x$peak.chr) == as.character(x$trait.chr) &
            !is.na(x$peak.chr) & !is.na(x$trait.chr))
    
    ## Cumulate attenuated score.
    attscore <- rep(0, nrow(x))
    if(any(tmp))
      attscore[tmp] <- x$peak.score[tmp] *
        exp(- 0.04 * abs(x$peak.pos.cM[tmp] - x$trait.pos.cM[tmp]))
    x$peak.attscore <- unlist(tapply(attscore, x$peak.chr,
                                     function(x) {
                                       out <- cumsum(x)
                                       out / max(out)
                                     }))
    
    ## Cumulate attenuated score.
    attscore <- x$peak.score
    if(any(tmp))
      attscore[tmp] <- x$peak.score[tmp] *
        (1 - exp(- 0.04 * abs(x$peak.pos.cM[tmp] - x$trait.pos.cM[tmp])))
    x$peak.enhscore <- unlist(tapply(attscore, x$peak.chr,
                                     function(x) {
                                       out <- cumsum(x)
                                       out / max(out)
                                     }))
  }

  x
}
#####################################################################
print.cistrans <- function(x, ...) print(summary(x, ...))
#####################################################################
print.summary.cistrans <- function(x, digits = 3, ...)
{
  use.cM <- attr(x, "use.cM")

  for(i in dimnames(x)[[2]]) {
    cat(paste("\n", i, ":\n", sep = ""))
    print(round(x[,i,], digits))
  }
  
  legend <- c(
      "peak.score (black) = score used to assess peak for trait\n",
      "peak.attscore (blue)  = peak.score reduced (exp decay with cM) away from cis position\n",
      "peak.enhscore  (red)   = peak.score reduced near cis position by (1 - attenuation)\n",
      "unweighted (green)  = unweighted estimate of most common position\n")
  cat("\n")
  for(i in dimnames(x)[[3]])
    cat(legend[pmatch(i, legend, nomatch = 0)])
  cat("reference (gray)   = reference line of no pattern over chromosome\n")
  cat("\nPositions are in", ifelse(use.cM, "cM", "Mb"), "\n")
  invisible()
}
#####################################################################
cumscore.dens <- function(object,
                          chr = chrs,
                          use.cM = FALSE,
                          window = 5,
                          ...)
{
  chrs <- table(object$peak.chr)
  chrs <- names(chrs[chrs > 0])
  
  ## Moving sum with window.
  pos <- "peak.pos.cM"
  scores <- c("unweighted","peak.enhscore","peak.attscore","peak.score")
  chr <- as.character(chr)
  window <- window / 2

  index <- paste(object$peak.chr, object$peak.pos.cM, sep = ":")
  index <- ordered(index, unique(index))
  dens <- object[!duplicated(index), c("peak.chr","peak.pos.cM","peak.pos.Mb")]
  dens.scores <- matrix(NA, nrow(dens), length(scores))
  dimnames(dens.scores) <- list(NULL, scores)

  drop <- rep(FALSE, length(scores))
  names(drop) <- scores
  for(score in scores[-1]) {
    drop[score] <- is.null(object[[score]])
    if(!drop[score])
      dens[[score]] <- tapply(object[[score]], index, sum, na.rm = TRUE)
  }
  scores <- scores[!drop]
      
  dens$unweighted <- table(index)

  tmpfn <- function(x, scores, pos, window) {
    scores <- scores[match(names(x), scores, nomatch = 0)]
    
    dens.scores <- matrix(NA, nrow(x), length(scores))
    dimnames(dens.scores) <- list(NULL, scores)

    index <- outer(x[[pos]], x[[pos]], function(x,y,w) abs(x-y)<= w, window)
    index <- apply(index, 1, function(x) range(seq(x)[x]))

    n.values <- sum(x$unweighted)
    for(score in scores) {
      weights <- cumsum(x[[score]])
      if(!is.null(weights)) if(nrow(x) > 1 & !any(is.na(weights))) {
        weights <- n.values * weights / weights[length(weights)]
        dens.scores[, score] <- weights[index[2,]] - weights[index[1,]]
      }
    }
    dens.scores
  }

  for(chri in chr) {
    ii <- dens$peak.chr == chri

    if(any(ii))
      dens[ii, scores] <- tmpfn(dens[ii, ], scores, pos, window)
  }
  
  dens
}
#####################################################################
summary.cistrans <- function(object, chr = levels(ordered(object$peak.chr)),
                             use.cM = FALSE,
                             dens = attr(object, "dens"), ...)
{
  if(is.null(dens))
    dens <- cumscore.dens(object, chr, use.cM)
  
  dens.peak <- function(dens, scores, pos) {
    if(is.null(dens))
      return(NULL)
    
    stats <- c("pos","left","right")
    poscnt <- matrix(NA, 2, length(scores))
    dimnames(poscnt) <- list(c("position","count"), scores)

    for(score in scores) {
      if(nrow(object) > 1) {
        wh <- which.max(dens[[score]])[1]
        poscnt["position", score] <- dens[[pos]][wh]
        poscnt["count", score] <- dens[[score]][wh]
      }
    }
    poscnt
  }

  pos <- paste("peak.pos", ifelse(use.cM, "cM", "Mb"), sep = ".")
  scores <- c("peak.enhscore","peak.attscore","peak.score","unweighted")
  scores <- scores[match(names(dens), scores, nomatch = 0)]
  chr <- as.character(chr)

  out <- array(NA, c(length(chr), 2, length(scores)))
  dimnames(out) <- list(chr, c(pos,"count"), scores)
  
  for(chri in chr) {
    tmp <- dens.peak(dens[chri == dens$peak.chr,], scores, pos)
    out[chri,, ] <- tmp
  }
  class(out) <- c("summary.cistrans","matrix")
  attr(out, "use.cM") <- use.cM
  out
}
#####################################################################
plot.cumscore <- function(x, peak.chr = NULL,
                          main = mains, xlim = xlims, use.cM = FALSE,
                          use.density = TRUE, window = 5,
                          maps = maps, ...)
{
  mains <- attr(x, "main")

  if(use.density) {
    dens <- attr(x, "dens")
    if(is.null(dens))
      dens <- cumscore.dens(x, window = window, ...)
  }
  
  if(!is.null(peak.chr))
    x <- x[x$peak.chr %in% peak.chr, ]
  x$peak.chr <- ordered(x$peak.chr)
  peak.chr <- levels(ordered(x$peak.chr))

  n.peak <- unique(x$peak.chr)
  n.peak <- length(peak.chr)
  if(nrow(x)) {
    if(n.peak > 1) {
      if(use.density) {
        form <- if(use.cM)
          formula("unweighted + peak.score ~ peak.pos.cM | peak.chr")
        else
          formula("unweighted + peak.score  ~ peak.pos.Mb | peak.chr")
        print(xyplot(form, dens,
                     type = "l", col = c("green","black"),
                     main = main, ...), ...)
      }
      else {
        print(xyplot(peak.cumcount + peak.cumscore ~ peak.pos.cM | peak.chr, x,
                     type = "l", col = c("green","black"),
                     panel = function(x,y,...) {
                       panel.lines(range(x, na.rm = TRUE), range(y, na.rm = TRUE), col = "gray")
                       panel.xyplot(x,y,...)
                     },
                     main = main, ...), ...)
      }
    }
    else {
      pos <- paste("peak.pos", ifelse(use.cM, "cM", "Mb"), sep = ".")
      xlims <- range(x[[pos]], na.rm = TRUE)
      if(use.density) {
        ylim <- range(dens$peak.score, dens$peak.enhscore, dens$peak.attscore, dens$unweighted, na.rm = TRUE)
        plot(dens[[pos]], dens$peak.score, type = "n", xlim = xlim, xaxt = "n",
             ylim = ylim,
             xlab = "", ylab = paste("traits with peaks in ", window, "cM window", sep = ""))

        abline(h = 0, col = "gray", lwd = 2)

        lines(dens[[pos]], dens$peak.score, lwd = 2)
        lines(dens[[pos]], dens$unweighted, lwd = 2, col = "green")
        if(!is.null(dens$peak.attscore))
           lines(dens[[pos]], dens$peak.attscore, lwd = 2, col = "blue")
        if(!is.null(dens$peak.enhscore))
           lines(dens[[pos]], dens$peak.enhscore, lwd = 2, col = "red")
      }
      else {
        ylims <- range(x$peak.cumscore)
        plot(x[[pos]], x$peak.cumscore, type = "n", xlim = xlim, xaxt = "n",
             xlab = "", ylab = "cumulative score")

        ## Add gray line for null hypothesis of no QTL pattern across chromosome.
        if(use.cM)
          lines(range(x$peak.pos.cM, na.rm = TRUE), c(0,1), col = "gray", lwd = 2)
        else
          lines(Mb.map[[peak.chr]], cM.map[[peak.chr]] / max(cM.map[[peak.chr]]), col = "gray", lwd = 2)

        lines(x[[pos]], x$peak.cumscore, lwd = 2)
        lines(x[[pos]], x$peak.cumcount, lwd = 2, col = "green")
        if(!is.null(x$peak.attscore))
           lines(x[[pos]], x$peak.attscore, lwd = 2, col = "blue")
        if(!is.null(x$peak.enhscore))
           lines(x[[pos]], x$peak.enhscore, lwd = 2, col = "red")
      }
      
      add.rug(peak.chr, main, maps, use.cM = use.cM, outer = TRUE,
              bottom.axis = TRUE, xlim = xlim)
    }
  }
}
######################################################################
plot.cis.dist <- function(x, peak.chr = NULL, n.qtl = NULL, cutoff = 2,
                          main = "", ...)
{
  require(lattice)
  
  if(!is.null(n.qtl))
    x <- x[x$n.qtl %in% n.qtl, ]

  if(is.null(peak.chr))
    peak.chr <- levels(x$trait.chr)

  tmp <- as.character(x$trait.chr) == as.character(x$peak.chr)
  if(length(peak.chr) == 1)
    tmp <- tmp & x$trait.chr == peak.chr

  x <- x[tmp,, drop = FALSE]

  ## Set up colors.
  colors <- mycol(x, ...)
  x$color <- colors$index
  trellis.par.set(superpose.symbol = list(col = colors$color))

  xyplot(trait.pos.Mb ~ peak.pos.Mb | trait.chr, x,
         group = color,
         panel = function(x,y,...) {
           panel.abline(0,1, col = "gray")
           panel.xyplot(x,y,...)
         }, main = main, ...)
}
