#################################################################
## Also see snp.R and qb.R.
#################################################################
## My chr-specific interpolation.
cM2Mb <- function(chr, pos)
{
  myapprox(cM.map[[chr]], Mb.map[[chr]], pos)$y
}
Mb2cM <- function(chr, pos)
{
  myapprox(Mb.map[[chr]], cM.map[[chr]], pos)$y
}
#################################################################
transcript.plot <- function(cross.name = "B6BTBR00",
                            tissue.name = "liver",
                            traitnames,
                            heatmap = length(traitnames) < 8,
                            chr = c(1:19,"X"), ...)
{
  if(cross.name == "B6BTBR00") {
    if(tissue.name != "liver")
      stop("only have liver for B6BTBR00 intercross")
    
    ## Be sure to use method = "hk" for this cross!
    require(B6BTBR00)
    data(list = paste(cross.name, tissue.name, sep = "."))
    newdata <- data.frame(MouseNum = dimnames(B6BTBR00.liver)[[2]])

    ## Add traits.
    tmp <- as.matrix(B6BTBR00.liver[traitnames,])
    if(length(traitnames) > 1)
      tmp <- t(tmp)
    newdata <- cbind(newdata, tmp)
    traitnames <- make.names(traitnames)
    names(newdata) <- c("MouseNum", traitnames)

    tmp <- match(6, chr)
    if(!is.na(tmp))
      chr <- chr[-tmp]
  }
  else {
    stop("no transcript data for other crosses (yet)")
  }

  tmp <- multraw(traitnames = traitnames, chr = chr,
    newdata = newdata, cross.name = cross.name, ...)

  plot(tmp, heatmap = heatmap, ...)
  summary(tmp)
}
#################################################################
## Umbrella routine. No longer used.
multtrait.plot <- function(traitnames = NULL,
                           main = "",
                           threshold.level = 0.05,
                           filename = NULL,
                           cross.name = "B6BTBR07",
                           ...)
{
    out <- multtrait(traitnames = traitnames,
                     main = main,
                     filename = filename,
                     threshold.level = threshold.level,
                     cross.name = cross.name, ...)
  print(summary(out))
  plot(out, ...)
}
################################################################
## My plot routine.
myplot <- function(traitnames = NULL, cross.name = "B6BTBR07",
                   ...)
{
  out <- multraw(traitnames = traitnames,
                 cross.name = cross.name, ...)
  print(summary(out))
  plot(out, ...)
}
###########################################################  
raw2mult <- function(x)
{
  ## Translate multraw to multtrait.
  ## Useful to get heatmaps.
}
#################################################################################
add.batch<-function(cross, tissue.name)
{
  ## Set up batch (put in RData later).
  ## Hyb.Date2 has proper levels.
  batch <- read.csv(file.path("Rosetta",
                              paste("F2", tissue.name, "hybsex.csv", sep = ".")),
                              header=TRUE, row.names=1)
  Hyb.Date <- levels(batch$Hyb.Date2)

  newdata <- matrix(NA, nind(cross), nlevels(batch$Hyb.Date2)-1)
  newdata2 <- matrix(0, dim(batch)[1], nlevels(batch$Hyb.Date2)-1)
  tmp <- match(as.character(batch$MouseNum),
               as.character(cross$pheno$MouseNum), nomatch=0)
  if (any(tmp==0))
    stop(paste("Mouse ID(s) not in list:",
               paste(dimnames(tissue)[[2]][tmp == 0], collapse = ",")))
  
  for (i in 1:(nlevels(batch$Hyb.Date2)-1))
    newdata2[as.numeric(batch$Hyb.Date2)==i,i] <- 1
  newdata[tmp,] <- newdata2
  dimnames(newdata) <- list(NULL,
                            paste('batch', seq(nlevels(batch$Hyb.Date2)-1),
                                  sep = ''))
  tmp <- match(dimnames(newdata)[[2]], names(cross$pheno))
  
  if (!all(is.na(tmp))) 
    stop(paste("trait", paste(dimnames(newdata)[[2]][!is.na(tmp)], collapse = ","), "already in database"))

  cross$pheno <- cbind(cross$pheno, newdata)
  cross   
}
#################################################################################
## My xyplot of cM to Mb.
show.maps <- function(chr = c(1:19,"X"), ...)
{
  data <- data.frame(chr = ordered(rep(c(1:19,"X"), sapply(cM.map, length)),
                                   c(1:19,"X")),
    cM = unlist(cM.map), Mb = unlist(Mb.map))
  row.names(data) <- unlist(lapply(cM.map, names))
  data$assembly <-
    snp.record$Build36.Assembly[match(row.names(data), snp.record$Build36.SNP)]

  require(lattice)
  xyplot(Mb ~ cM | chr, data[data$chr %in% chr, ], group = assembly,
    col = c("blue","red"),
    panel = function(x,y,...) {
      panel.xyplot(x,y,...)
      p <- myapprox(y, x)
      panel.lines(p$y, p$x, col = "gray")
    }, ...)
}
######################################################################
simplify.pattern <- function(pattern)
{
  ## Split on comma separating chr and epistatic pairs of chr.
  x <- strsplit(as.character(pattern), ",", fixed = TRUE)

  tmpfn <- function(x) {
    ## Drop epistatic interactions.
    g <- grep(".*:.*", x)
    if(length(g))
      x <- x[-g]
    ## Drop split on chromosomes.
    x <- sapply(x, function(x) strsplit(x, ".", fixed = TRUE)[[1]][1])
    paste(x, collapse = ",")
  }
  x <- sapply(x, tmpfn)
  x[x == "NULL"] <- ""
  x
}
######################################################################
make.cispeaks <- function(cross, islet.dist, filename = "islet.bim.cis.csv", ...)
{
  ## Get peak.locus. Need to find (pseudo)marker closest to peak.
  tmp <- mypull.loci(cross, ...)
  ## Match where we can.
  tmp2 <- match(paste(islet.dist$peak.chr, islet.dist$peak.pos.cM, sep = "."),
                paste(tmp$chr, tmp$pos, sep = "."))
  tmp3 <- is.na(tmp2)
  for(i in levels(islet.dist$peak.chr)) {
    ii <- (islet.dist$peak.chr == i) & tmp3
    jj <- (tmp$chr == i)
    if(any(ii) & any(jj)) {
      print(i)
      ## Now we try to match.
      for(j in seq(sum(ii)))
        tmp2[ii][j] <- which.min(abs(islet.dist$peak.pos.cM[ii][j] - tmp$pos[jj]))
    }
  }
  islet.dist$peak.locus <- row.names(tmp)[tmp2]
  islet.peaks <- islet.dist[, c("gene.id","peak.locus","peak.chr","peak.pos.cM","peak.pos.Mb","peak.score")]
  names(islet.peaks)[1] <- "a_gene_id"
  row.names(islet.peaks) <- seq(nrow(islet.peaks))
  write.csv(islet.peaks, file = filename)
}
######################################################################
cis.dist <- function(tissue.best, gene.annot, maps, geno.names = 1:19,
                     tissue = "islet", cis.call = 5)
{
  ## Get chr (drop "_random" stuff).
  ## NB: Some out.trans can be missing (unknown position of SNP?).
  chr.trans <- sapply(strsplit(as.character(gene.annot$Chromosome), "_",
                         fixed = TRUE),
                function(x) x[1])

  ## Drop null models.
  tissue.best <- tissue.best[tissue.best$n.qtl > 0, ]
  ## Match a_gene_id between tissue.best and gene.annot.
  aGeneId <- match(as.character(tissue.best$a_gene_id),
                   as.character(gene.annot$a_gene_id),
                   nomatch = 0)
  chr.trans <- ordered(chr.trans[aGeneId], geno.names)
  
  chr.best <- as.character(tissue.best$chrom[aGeneId > 0])
  chr.best <- sapply(strsplit(chr.best, ".", fixed = TRUE),
                     function(x) x[1])
  chr.best <- ordered(chr.best, geno.names)
  n.qtl <- c(table(tissue.best$a_gene_id))
  n.qtl <- rep(n.qtl, n.qtl)
    
  ## NB: pos is gene position in Mb; locus is mapped locus in cM.
  ## Make output compatible with cistrans.
  if(length(tissue) == 1)
    tissue <- rep(tissue, sum(aGeneId > 0))
  
  out <- data.frame(Tissue = tissue,
                    peak.chr = chr.best,
                    peak.pos.cM = tissue.best$locus[aGeneId > 0],
                    peak.score = tissue.best$variance[aGeneId > 0],
                    symbol = gene.annot$Symbol[aGeneId],
                    gene.id = gene.annot$a_gene_id[aGeneId],
                    trans.chr = chr.trans,
                    trans.pos.Mb = gene.annot$Chromosome_Position[aGeneId])
  o <- order(out$trans.chr, out$trans.pos.Mb)
  out <- out[o, ]

  ## Recast peak.pos as Mb and trans.pos as cM for later use.
  out$peak.pos.Mb <- mypos(out, maps, "cM", "peak")
  out$trans.pos.cM <- mypos(out, maps, "Mb", "trans")

  out$cis <- abs(out$peak.pos.Mb - out$trans.pos.Mb) < cis.call
  out$gene.name..link.to.MGI. <- gene.annot$Description[aGeneId][o]
  out$n.qtl <- n.qtl[o]

  class(out) <- c("cis.dist", "data.frame")
  out
}
