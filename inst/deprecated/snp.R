## These routines are specific to SNP records.
######################################################################
plot.snp.record <- function(x = snp.record[[diagname]],
                            diagname = deparse(substitute(x)),
                            chr = c(1:19,"X"),
                            title = {
                              if(length(chr) == 1)
                                paste("chr ", chr, ":", diagname)
                              else
                                diagname
                            },
                            map,
                            ...)
{
  if(is.character(x) & length(x) == 1) {
    diagname <- x
    x <- snp.record[[diagname]]
  }
  if(is.null(x))
    stop("invalid x or diagname")
  
  require(qtl)
  
  markers <- unlist(lapply(map, names))
  
  scanone.obj <- data.frame(chr = ordered(rep(names(map), sapply(map, length)),
                              names(map)),
                            pos = unlist(map))
  tmp <- rep(NA, nrow(scanone.obj))
  tmp2 <- match(snp.record$Build36.SNP, markers, nomatch = 0)

  tmp[tmp2] <- x[tmp2 > 0]
  scanone.obj[[diagname]] <- tmp
  class(scanone.obj) <- c("scanone","data.frame")

  plot(scanone.obj, chr=chr, ...)

  if(length(chr) == 1) {
    rug(cM.map[[chr]],0.02,side=3)
    p <- myapprox(Mb.map[[chr]], cM.map[[chr]])
    axis(3, p$y, p$x)
    mtext("Mb", 3, 1, at = -5, adj = 1)
  }

  title(title, line = 2)
  invisible()
}

######################################################################
## Examine segregation disorder around single SNPs.
find.segdisorder <- function(segdis, position, maxdist = 5)
{
  ## First have to sort by position.
  o <- order(position)
  segdis <- segdis[o]

  segdis[is.na(segdis)] <- 0
  
  n.segdis <- length(segdis)
  left <- c(segdis[-1], segdis[n.segdis])
  right <- c(segdis[1], segdis[-n.segdis])

  position <- position[o]
  position <- diff(position)
  position <- pmax(c(0, position), c(position, 0))

  res <- rep(0, length(segdis))
  res[o] <- pmax(0, segdis - left, segdis - right) * (position < maxdist)
  res
}
find.segdischr <- function()
{
  double.snp <- list()
  chroms <- levels(snp.record$Build36.Chrom)
  for(i in chroms) {
    tmp2 <- snp.record$Build36.Chrom == i & snp.record$Segregating & snp.record$Keep
    tmp3 <- find.segdisorder(snp.record$SegDisorder.logP[tmp2],
                             snp.record$Build36.Position[tmp2])
    
    ## Double crossovers by SNP and mouse.
    double.snp[[i]] <- list(snp = tmp3, seq = seq(tmp2)[tmp2])
  }
  ## Reorder snp values into res$snp.
  res <- rep(0, nrow(snp.record))
  res[unlist(lapply(double.snp, function(x) x$seq))] <-
    unlist(lapply(double.snp, function(x) x$snp))
  res
}

######################################################################
## Examine crossover pattern around single SNPs.
find.crosspat <- function(geno, position, maxdist = 5)
{
  ## First have to sort by position.
  o <- order(position)
  geno <- geno[o]
  
  geno[geno > 3] <- NA

  n.geno <- length(geno)
  left <- c(geno[-1], NA)
  right <- c(NA, geno[-n.geno])

  posleft <- position[o]
  posleft <- diff(posleft) < maxdist
  posright <- c(posleft, FALSE)
  posleft <- c(FALSE, posleft)

  ## Add up crossovers to left and right of SNP.
  ## Count  meaning
  ## 0      no change or flanking too far away.
  ## 1      single crossover from left.
  ## 2      double crossover from left.
  ## 4      single crossover from right.
  ## 8      double crossover from right.
  ## 5      single crossover from left and right (double crossover around SNP).
  ## 6,9,10 single and double or two doubles.
  tmp1 <- abs(geno - left) * posleft
  tmp1[is.na(tmp1)] <- 0
  tmp2 <- abs(geno - right) * posright
  tmp2[is.na(tmp2)] <- 0
  count <- tmp1 + 4 * tmp2
  ## 41     single CO
  ## 82     double between two SNPs (DC2)
  ## 451    single on both sides (DC1): 451 = 410 + 041

  ## 861    single right, double left: 861 = 820 + 041
  ## 492    double right, single left: 492 = 410 + 082
  ## 4551   three singles: 4551 = 4100 + 0410 + 0041

  ## Flanking SNPs.
  left <- c(count[-1], 0)
  right <- c(0, count[-n.geno])

  ## Pattern ..4[59].. and ..[56]1.. indicate double crossover around SNP.
  ## Don't count single CO: change to ..0[59].. and ..[56]0..
  redo <- left == 4 & (count %in% c(5,9))
  count[c(redo[-1], FALSE)] <- 0
  redo <- right == 1 & (count %in% c(5,6))
  count[c(FALSE, redo[-n.geno])] <- 0

  count[count %in% c(6,9,10)] <- 9

  res <- rep(0, length(geno))
  res[o] <- count

  ## Summaries by mouse across chromosome.
  ## Pattern ..41..  indicates single crossover between two SNPs.
  ##                   count for both SNPs, but once per mouse.
  tmp <- ceiling(sum(count %in% c(1,4)) / 2)
  ## Pattern ..5..   indicates two t
  tmp[2] <- sum(count == 5)
  ## Pattern ..82..  indicates double crossover between two SNPs.
  ##                   count for both SNPs, but once per mouse.
  tmp[3] <- ceiling(sum(count %in% c(2,8)) / 2)
  ## Other patterns for mouse.
  tmp[4] <- sum(count == 9)
              
  c(tmp, res)
}

######################################################################
find.patternprob <- function(snp.record, mouse.record, fname, type,
                           build = 32, ...)
{
  if(build == 32) {
    chrom.name <- "Chrom.Name"
    chrom.pos <- "Chrom.Position"
  }
  else { ## build 36
    chrom.name <- "Build36.Chrom"
    chrom.pos <- "Build36.Position"
  }
  double.snp <- list()
  chroms <- levels(snp.record[[chrom.name]])
  snp.scores <- c(1,4,5,2,8,9)
  names(snp.scores) <- paste("snp", c("COL","COR","DC1","DCL","DCR","DD"), sep = ".")
  for(i in chroms) {
    tmp2 <- (snp.record[[chrom.name]] == i &
             snp.record$Segregating & snp.record$Keep)
    tmp <- apply(snp.f2[tmp2, ], 2, fname, snp.record[[chrom.pos]][tmp2], ...)

    ## Double crossovers by SNP and mouse.
    ## Pull out single and double count per mouse.
    double.snp[[i]] <- list(mouse.CO = tmp[1, ],
                            mouse.DC1 = tmp[2, ],
                            mouse.DC2 = tmp[3, ],
                            mouse.DD = tmp[4, ],
                            seq = seq(tmp2)[tmp2])
    
    tmp <- tmp[-(1:4), ]
    tmp3 <- mouse.record$Keep[mouse.record$Strain == "F2"]

    for(j in names(snp.scores))
      double.snp[[i]][[j]] <- apply(tmp[, tmp3], 1,
                                    function(x) sum(x == snp.scores[j]))
  }
  ## Reorder snp values into res$snp.
  res <- list()
  tmp2 <- unlist(lapply(double.snp, function(x) x$seq))
  for(type in names(snp.scores)) {
    res[[type]] = rep(0, nrow(snp.record))
    res[[type]][tmp2] <- unlist(lapply(double.snp, function(x) x[[type]]))
  }

  mouse.scores <- paste("mouse", c("CO","DC1","DC2","DD"), sep = ".")
  ## Keep mouse values separately by chr.
  for(type in mouse.scores) {
    res[[type]] <- as.data.frame(matrix(0, nrow(mouse.record), length(chroms)))
    names(res[[type]]) <- paste(type, chroms, sep = ".")
    res[[type]][mouse.record$Strain == "F2",] <-
      as.data.frame(lapply(double.snp, function(x) x[[type]]))
  }
  res
}

######################################################################
## Examine double crossovers around single SNPs.
find.doublecross <- function(geno, position, maxdist = 5)
{
  ## First have to sort by position.
  o <- order(position)
  geno <- geno[o]
  
  geno[is.na(geno)] <- 9

  n.geno <- length(geno)
  left <- c(geno[-1], 9)
  right <- c(9, geno[-n.geno])

  position <- position[o]
  position <- diff(position)
  position <- pmax(c(0, position), c(position, 0))

  res <- rep(0, length(geno))
  res[o] <- (geno - left != 0 & geno - right != 0 &
             geno < 4 & left < 4 & right < 4 &
             position < maxdist)
  res
}

######################################################################
## Want to also look at crossovers 0->2 and 2->0 at adjacent SNPs.
find.crossover <- function(geno, position, maxdist = 5, change = 2)
{
  ## First have to sort by position.
  o <- order(position)
  geno <- geno[o]
  
  position <- position[o]
  position <- c(diff(position), 0)

  geno[is.na(geno)] <- 9
  
  n.geno <- length(geno)
  right <- c(9, geno[-n.geno])

  res <- rep(0, length(geno))
  res[o] <- (abs(geno - right) == change &
             geno < 4 & right < 4 &
             position < maxdist)
  res
}

######################################################################
find.crossprob <- function(snp.record, mouse.record, fname, type,
                           build = 32, pattern = FALSE, ...)
{
  if(build == 32) {
    chrom.name <- "Chrom.Name"
    chrom.pos <- "Chrom.Position"
  }
  else { ## build 36
    chrom.name <- "Build36.Chrom"
    chrom.pos <- "Build36.Position"
  }
  double.snp <- list()
  chroms <- levels(snp.record[[chrom.name]])
  for(i in chroms) {
    tmp2 <- snp.record[[chrom.name]] == i & snp.record$Segregating & snp.record$Keep
    tmp <- apply(snp.f2[tmp2, ], 2, fname, snp.record[[chrom.pos]][tmp2], ...)

    if(pattern) {
      single <- tmp[1, ]
      double <- tmp[2, ]
      tmp <- tmp[-(1:2), ]
    }
    ## Double crossovers by SNP and mouse.
    double.snp[[i]] <- list(snp = apply(tmp, 1, sum),
                            seq = seq(tmp2)[tmp2],
                            mouse = apply(tmp, 2, sum))
  }
  ## Reorder snp values into res$snp.
  res <- list(snp = rep(0, nrow(snp.record)))
  res$snp[unlist(lapply(double.snp, function(x) x$seq))] <-
    unlist(lapply(double.snp, function(x) x$snp))

  ## Keep mouse values separately by chr.
  res$mouse <- as.data.frame(matrix(0, nrow(mouse.record), length(chroms)))
  names(res$mouse) <- paste(type, chroms, sep = ".")
  res$mouse[mouse.record$Strain == "F2",] <-
    as.data.frame(lapply(double.snp, function(x) x$mouse))
  res
}

#################################################################
## Make cross object.
snp2cross <- function(snp.geno, mouse.pheno, snp.id, include.x = TRUE,
                      build36 = TRUE)
{
  ## Recode Sex to remove extra level.
  mouse.pheno$Sex <- factor(mouse.pheno$Sex)
  
  ## Recode any errors as missing values.
  snp.geno[snp.geno > 3] <- NA
  
  geno <- list()
  if(build36) {
    chrom <- snp.id$Build36.Chrom
    pos <- snp.id$Build36.Position
    snp <- snp.id$Build36.SNP
  }
  else {
    chrom <- snp.id$Chrom.Name
    pos <- snp.id$Chrom.Position
    snp <- snp.id$External.Id
  }
  o <- order(chrom, pos)
  snp.geno <- snp.geno[o,]
  chrom <- chrom[o]
  pos <- pos[o]
  snp <- as.character(snp[o])

  for(i in levels(chrom)) {
    if(i != "X" | include.x) {
      ii <- chrom == i
      data <- as.matrix(t(snp.geno[ii, ]))
      dimnames(data) <- list(seq(nrow(data)), snp[ii])
      map <- pos[ii]
      names(map) <- snp[ii]
      geno[[i]] <- list(data = data, map = map)
      class(geno[[i]]) <- ifelse(i == "X", "X", "A")
    }
  }
  cross <- list(geno = geno, pheno = mouse.pheno)
  class(cross) <- c("f2", "cross")
  cross
}
