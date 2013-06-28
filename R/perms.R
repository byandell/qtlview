#############################################################
read.perms <- function(cross, filenames, row.names = NULL,
                    perms = read.permfile(filenames, row.names, ...),
                    method = "hk", model = "normal", ...)
{
  ## The acquisition is not documented.
  ## Need to check email with Jee Young Moon.
  ## Need to allow for X chr as well as both, male, female.
  ## Could need to read several files?
  
  read.permfile <- function(filenames, row.names, ...) {
    ## Create a list of data frames. Each data frame should have 
    sexes <- c("both","female","male")
    if(!all(names(filenames) %in% sexes))
      stop("read.perms filenames need to be given as sex-named character vector")
  
    perms <- list()
    for(sex in names(filenames)) {
      filename <- filenames[sex]
      if(!missing(filename)) {
        res <- read.csv(filename, header = TRUE, row.names = row.names, ...)

        ## Drop index column if present.
        if(ncol(res) > 2 | all(res[[1]] == seq(nrow(res))))
          res <- res[, -1, drop = FALSE]
        ## Make sure X is not X.1.
        names(res)[names(res) == "X.1"] <- "X"
      }
      perms[[sex]] <- res
    }
    perms
  }

  for(sex in names(perms)) {
    res <- perms[[sex]]
    
    ## Recast perm as scanoneperm object. See [qtl]scanone.perm.
    xchr <- names(res) == "X"
    is.X <- any(xchr) && any(names(res) == "A") & ncol(res) > 1
    
    df <- rep(1, ncol(res))
    names(df) <- names(res)
    
    if(is.X) {
      ## Looks like X-chr specific permutations performed.
      res <- lapply(res, function(x) matrix(x, ncol = 1,
                                            dimnames = list(NULL, "lod")))
      
      attr(res, "xchr") <- xchr
      
      L <- sapply(pull.map(cross), function(x) diff(range(x)))
      tmp <- names(L) == "X"
      La <- sum(L[!tmp])
      Lx <- sum(L[tmp])
      attr(res, "L") <- c(A = La, X = Lx)
      attr(res, "df") <- rbind(df, df)
    }
    else {
      res <- as.matrix(res)
      attr(res, "df") <- df
    }
    attr(res, "method") <- method
    attr(res, "model") <- model
    attr(res, "type") <- class(cross)[1]
    if(is.X)
      class(res) <- c("scanoneperm", "list")
    else
      class(res) <- c("scanoneperm", "matrix")
    
    perms[[sex]] <- res
  }
  class(perms) <- c("read.perms", "list")
  perms
}
#############################################################
plot.read.perms <- function(x, ...)
{
  ## Caution: plot.scanoneperm forces mfrow=c(2,1).
  for(sex in names(x))
    plot(x[[sex]], ...)
  invisible()
}
#############################################################
summary.read.perms <- function(object, ...)
{
  ## Assumption is that each permutation run is for one trait.
  
  sum <- lapply(object, summary, ...)
  out <- sum[[1]]
  if(is.matrix(out)) {
    ## No X chromosome; summary.scanoneperm object is matrix.
    for(i in names(sum)[-1])
      out <- cbind(out, sum[[i]])
    dimnames(out)[[2]] <- names(sum)
    attr(out, "n.perm") <- attr(sum[[1]], "n.perm")
    class(out) <- "summary.scanoneperm"
  }
  else {
    ## Autosome and X chromosome; summary.scanoneperm object is data frame.
    types <- names(out)
    for(type in types) {
      tbl <- as.data.frame(lapply(sum, function(x) x[[type]]))
      names(tbl) <- names(sum)
      out[[type]] <- tbl
    }
  }
  out
}
#############################################################
threshold.perm <- function(perms, sex, threshold.level, type = "A")
{
  threshold.lod <- 0
  names(threshold.lod) <- "A"
  if(threshold.level > 0 & threshold.level < 1) {
    tmp <- summary(perms[[sex]], threshold.level)
    if("A" %in% type)
      threshold.lod["A"] <- ifelse(is.list(tmp), tmp$A, tmp[1])
    if("X" %in% type) {
      if(length(type) == 1)
        threshold.lod <- NULL
      threshold.lod["X"] <- tmp$X
    }
  }
  threshold.lod
}
################################################################
find.threshold <- function(scans, chr = "", threshold.lod)
{
  chr <- find.chr(chr, levels(scans$chr))

  if(any(threshold.lod > 0)) {
    ## Allow for separate thresholds for autosomes and X chromosome.
    thresholds <- rep(threshold.lod[1], length(chr))
    names(thresholds) <- chr
    if(!is.na(threshold.lod["X"]) & "X" %in% chr)
      thresholds["X"] <- threshold.lod["X"]

    ## Find what chr cross threshold.
    keep.chr <- character()
    for(i in chr) {
      ii <- scans$chr == i
      if(any(ii)) {
        tmp <- apply(scans[ii, -(1:2), drop = FALSE], 2,
                     function(x) any(x > thresholds[i]))
        if(any(tmp)) {
          keep.chr <- c(keep.chr, i)
        }
      }
    }
  }
  else
    keep.chr <- chr
  
  ## Now find what traits cross the threshold.
  if(length(keep.chr)) {
    if(threshold.lod > 0) {
      traits <- apply(scans[scans$chr %in% keep.chr, -(1:2),
                            drop = FALSE], 2,
                      function(x) any(x > threshold.lod))
      traits[is.na(traits)] <- FALSE
    }
    else
      traits <- rep(TRUE, ncol(scans) - 2)
  }
  else
    traits <- NULL

  list(chr = keep.chr, traits = traits)
}
################################################################
threshold.lines <- function(scans, threshold.lod, gap = 25,
                           lwd = 2, lty = 2, col.scheme, cluster, ...)
{
  ## Assume scans has been reduced to those chrs being plotted.
  chrs <- levels(scans$chr)
  
  if(any(threshold.lod > 0)) {
    if(!("X" %in% chrs) | all(chrs %in% "X") | (length(threshold.lod) == 1)) {
      if(!("X" %in% chrs)) {
        ## Single threshold for all autosomes.
        if(length(threshold.lod) > 1)
          threshold.lod <- threshold.lod["A"]
      }
      else {
        ## Only X chromosome.
        if(!sum("X" != chrs)) {
          if(length(threshold.lod) > 1)
            threshold.lod <- threshold.lod["X"]
        }
      }
      ## One threshold line.
      abline(h = threshold.lod, lwd = lwd, lty = lty, ...)
    }
    else {
      ## Need to separate A and X thresholds.
      chr.len <- tapply(scans$pos,scans$chr, function(x) diff(range(x)))
      ## Assume X is at the end.
      is.X <- names(chr.len) == "X"
      n.chr <- length(chrs)
      left <- c(0, cumsum(gap + chr.len[-n.chr]))
      names(left) <- names(chr.len)
      right <- (cumsum(gap + chr.len)) - gap
      lines(range(left[!is.X], right[!is.X]), rep(threshold.lod["A"], 2),
            lwd = lwd, lty = lty, ...)
      lines(range(left[is.X], right[is.X]), rep(threshold.lod["X"], 2),
            lwd = lwd, lty = lty, ...)
    }
  }
}
################################################################
threshold.pass.test <- function(x, lod)
{
  ifelse(length(x) > 0,
         ifelse(all(is.na(x)), FALSE,
                max(x, na.rm = TRUE) >= lod),
         FALSE)
}
################################################################
threshold.pass <- function(scans, threshold.lod, chr, dim = 2)
{
  ## Find traits that pass threshold.
  if(length(threshold.lod) == 1 | !("X" %in% levels(chr)))
    apply(scans, dim, threshold.pass.test, threshold.lod[1])
  else {
    ## Check both Autosomes and X chromosome.
    apply(scans, dim,
          function(x, chr, lod) {
            threshold.pass.test(x[chr != "X"], lod["A"]) |
            threshold.pass.test(x[chr == "X"], lod["X"])
          },
          chr, threshold.lod)
  }
}

