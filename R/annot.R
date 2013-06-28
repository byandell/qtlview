###############################################################################
## Key annotation routine used by package qtlmult.
## Object trait.annotation must be provided (see read.annotation below)
## or the annotation list is returned as NULL.

##################################################################
myget <- function(cross.name, suffix)
{
  tmp <- paste(cross.name, suffix, sep = ".")
  if(exists(tmp))
    get(tmp)
  else
    NULL
}
##################################################################
add.trait.position <- function(out, trait.position)
{
  if(is.null(trait.position))
    out
  
  traitnames <- row.names(out)
  tmp <- trait.position[traitnames, ]
  tmp2 <- sapply(tmp, function(x) all(is.na(x)))
  if(all(tmp2))
    out
  else {
    tmp <- as.data.frame(lapply(tmp, function(x)
                                {
                                  x <- as.character(x)
                                  x[is.na(x)] <- "."
                                  factor(x)
                                }))
    cbind(out, tmp)
  }
}
##################################################################
find.trait.position <- function(traitnames, maps,
                                trait.annotation = NULL,
                                pos = c("cM", "Mb"),
                                prefix = "trait", digits = 3, ...)
{
  ## Find match of traitname and trait.annotation.
  is.selected <- find.trait.annot(unique(traitnames), trait.annotation, ...)
  n.pos <- sum(is.selected > 0)
  if(n.pos) {
    trait.position <- data.frame(
      chr = as.character(trait.annotation$Chromosome[is.selected]))
    row.names(trait.position) <- unique(traitnames)[is.selected > 0]

    ## Set up physical positions for correspondence.
    Mb.pos <- trait.annotation$Chromosome_Position[is.selected]
    chrs <- names(maps$Mb.map)

    is.selected <- is.na(match(trait.position$chr, chrs))
    if("Mb" %in% pos) {
      trait.position$Mb <- round(Mb.pos, digits)
    }
    if("cM" %in% pos) {
      trait.position$cM <- rep(NA, n.pos)

      for(chri in chrs) {
        is.selected <- (trait.position$chr == chri &
                        !is.na(trait.position$chr))
        if(any(is.selected))
          trait.position$cM[is.selected] <-
            round(qm.approx(maps, "Mb", chri, Mb.pos[is.selected])$y, digits)
      }
    }
    if(prefix != "")
      names(trait.position) <- paste(prefix, names(trait.position), sep = ".")
    attr(trait.position, "prefix") <- prefix
    trait.position
  }
  else
    NULL
}
##########################################################
## Create trait.annotation object.
## Used in find.trait.position.
read.annotation <- function(filename, update.names = NULL,
                          drop.extra = TRUE, ...)
{
  trait.annotation <- read.csv(filename, header = TRUE, ...)
  annot.names <- names(trait.annotation)

  ## Need a_gene_id, Symbol, Chromosome, Start_Coordinate, End_Coordinate, Strand
  ## Chromosome_Position  = Start_Coordinate / 10^6 (if missing)
  ## Start_Coordinate, End_Coordinate are then dropped
  ## Only the following are curently used in package qtlview:
  ## a_gene_id,Symbol,Chromosome,Chromosome_Position
  
  final.names <- c("a_gene_id", "Symbol", "Chromosome",
                   "Start_Coordinate", "End_Coordinate", "Strand")

  if(length(update.names)) {
    ## Check that update names are in final names.
    m <- pmatch(tolower(names(update.names)), tolower(final.names))
    if(any(is.na(m)))
      stop(paste("annotation internal names do not match:",
                 paste(names(update.names)[is.na(m)], collapse = ", ")))
    names(update.names) <- final.names[m]

    ## Check that values of update.names match column names of file.
    m <- pmatch(update.names, annot.names)
    if(any(is.na(m)))
      stop(paste("annotation column names do not match:",
                 paste(update.names[is.na(m)], collapse = ", ")))
    
    annot.names[m] <- names(update.names)
    
    names(trait.annotation) <- annot.names
  }

  ## Make sure all final names are on trait.annotation.
  tmp <- is.na(match(final.names, annot.names))
  if(any(tmp)) {
    warning(paste("missing trait annotation names:",
                  paste(final.names[tmp], collapse = ", ")))
    final.names <- final.names[!tmp]
  }

  if(drop.extra)
    trait.annotation <- trait.annotation[, final.names]

  if(is.null(trait.annotation$Chromosome_Position))
    trait.annotation$Chromosome_Position <- 10^-6 * trait.annotation$Start_Coordinate
#    trait.annotation$Chromosome_Position <- 0.5 * 10^-6 *
#      (trait.annotation$Start_Coordinate + trait.annotation$End_Coordinate)
  if(drop.extra)
    trait.annotation$Start_Coordinate <- trait.annotation$End_Coordinate <- NULL

  trait.annotation
}

################################################################
## Used here only.
get.geneid <- function(traitnames, trait.annotation, blank = TRUE, ...)
{
  tissues <- find.trait.annot(traitnames, trait.annotation, out.type = "tissue", ...)

  geneid <- sapply(strsplit(traitnames, ".", fixed = TRUE),
                   function(x) x[length(x)])
  ## Find which are actually a_gene_id's.
  tmp <- geneid %in% c(trait.annotation$a_gene_id, make.names(trait.annotation$a_gene_id))
  if(blank)
    geneid[!tmp] <- ""
  else {
    if(any(!tmp))
      geneid[!tmp] <- traitnames[!tmp]
    if(any(tmp)) {
      ## Add tissue name if more than one used.
      if(length(table(tissues[tmp])) > 1)
        geneid[tmp] <- paste(tissues[tmp], geneid[tmp], sep = ".")
    }
  }
      
  geneid
}
################################################################
## Used here only.
get.symid <- function(geneid, trait.annotation, catenate = TRUE)
{
  tmp <- geneid %in% trait.annotation$a_gene_id
  if(any(tmp)) {
    tmp2 <- match(geneid[tmp], trait.annotation$a_gene_id)
    geneid[tmp] <- if(catenate)
      paste(trait.annotation$Symbol[tmp2],
            trait.annotation$a_gene_id[tmp2], sep = ".")
    else
      as.character(trait.annotation$Symbol[tmp2])
  }
  tmp <- is.na(geneid)
  if(any(tmp))
    geneid[tmp] <- ""
  geneid
}
################################################################
mylabels <- function(traitnames, max.names = 200, ylab){
  if(ylab == "none" | length(traitnames) > max.names)
    return(TRUE)
  
  if(ylab == "symbol.a_gene_id")
    return(traitnames)
  if(ylab == "symbol" | ylab == "a_gene_id"){
    traitspl <- strsplit(traitnames,"\\.")
    if(ylab == "symbol")
      return(sapply(traitspl,function(x) paste(x[1],x[2],sep=".")))
    else
      return(sapply(traitspl, function(x) paste(x[1],x[length(x)],sep=".")))
  }
}
###############################################################################
find.trait.annot <- function(traitnames, trait.annotation, out.type = "selected",
                             summarize = FALSE,
                             tissue.list = NULL,
                             ...)
{
  if(is.null(trait.annotation))
    return(NULL)
  
  ## This routine has multiple uses:
  ## out.type="selected": Find a_gene_id as last part of traitnames (sep = ".").
  ## out.type="tissues":  Find tissues as first part of traitnames (set to unknown if missing).
  ## summarize=TRUE: Summarise tissues for use in main plot title.
  if(out.type == "selected") {
    ## a_gene_id is last part of traitnames for transcripts.
    tmp <- sapply(strsplit(traitnames, ".", fixed = TRUE),
                  function(x) x[length(x)])
    tissues <- match(tmp, trait.annotation$a_gene_id, nomatch = 0)
  }
  else{
    ## Tissue is first part of traitnames.
    tissues <- tolower(sapply(strsplit(traitnames, ".", fixed = TRUE),
                   function(x) x[1]))
    if(!is.null(tissue.list)) {
      tmp <- !(tissues %in% tissue.list)
      if(any(tmp))
        tissues[tmp] <- "unknown"
    }
  }
  if(summarize) {
    tissues <- table(tissues)
    tissues <- paste(tissues, names(tissues), collapse = ", ")
  }
  tissues
}
