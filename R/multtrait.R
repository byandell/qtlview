###########################################################
## Routine multtrait assumes cross, perms and maps are available.
## Default names are cross.name, cross.name.perms, cross.name.maps.

multtrait <- function(traitnames = NULL,
                      cross = get(cross.name),
                      cross.name = deparse(substitute(cross)),
                      perms = myget(cross.name, "perms"),
                      maps = myget(cross.name, "maps"),
                      category="clinical",
                      filename=NULL,
                      chr="",
                      main = "",
                      summary.plot = c("total","chr","none","phenotype"),
                      threshold.level = 0.05,
                      threshold.lod = calc.lod,
                      cluster = TRUE,
                      support.lod = 1.5,
                      n.clust = ifelse(is.null(hc), 0, ceiling(log2(length(hc$order)))),
                      heatmap = FALSE,
                      scan.type = c("LOD","LPD","BF","RAW","MOM","PAT","BIM"), 
                      sex = c("both","male","female"),
                      ylab = c("symbol","a_gene_id","symbol.a_gene_id","none","on"),
                      trait.annotation = myget(cross.name, "annotation"),
                      ...)
{
  ## Match arguments.
  sex <- match.arg(sex)
  scan.type <- match.arg(scan.type)
  if(scan.type == "PAT")
    scan.type <- "BIM"
  ylab <- match.arg(ylab)
  if(ylab == "on")
    ylab <- "none"

  ## Construct scanone from raw data (much slower).
  if(is.raws(filename)) {
    ## Really should have called multraw (formerly myplot).
    if(scan.type == "RAW")
      scan.type <- "LOD"
    return(multraw(traitnames = traitnames,
                   cross = cross, perms = perms, maps = maps,
                   main = main,
                   filename = filename,
                   threshold.level = threshold.level,
                   chr = chr,
                   scan.type = scan.type, sex = sex, ylab = ylab,
                   trait.annotation = trait.annotation,
                   ...))
  }

  ## Threshold for LOD: calc.lod is default for threshold.lod.
  if(scan.type == "LOD") {
    calc.lod <- threshold.perm(perms, sex, threshold.level)
  }
  else {
    ## No meaningful threshold for non-LOD approaches.
    support.lod <- 0.5
    calc.lod <- 1e-6
  }
  
  ## Set up pre-computed scans.
  scans <- get.scans(traitnames, category, filename, cross,
                     maps, chr, threshold.lod, sex,
                     trait.annotation = trait.annotation, ...)
  is.selected <- attr(scans, "is.selected")
  
  ## Extract category, tissue.name, traitnames.
  category <- attr(scans, "category")
  tissue.name <- attr(scans, "tissue.name")
  traitnames <- names(scans)[-(1:2)]
  tmp <- length(traitnames)
  if(tmp == 0) {
    tmp <- paste("\n\n*** No traits found in file", filename, "***\n\n")
    cat(tmp)
    stop(tmp)
  }

  ## Set up default title if needed.
  if(main == "")
    main <- paste(cross.name, ": ", length(traitnames), " phenotypes", sep = "")

  ## Set up hierarchical clustering to order traits if cluster=TRUE.
  hc <- calc.hc(scans[, c(TRUE, TRUE, is.selected)], cluster)

  summary.plot <- match.arg(summary.plot)
  if(summary.plot == "total" & any(is.selected))
    comp.lod <- composite.lod(scans, is.selected, main, threshold.lod, support.lod,
                              scan.type = scan.type)
  else
    comp.lod <- NULL

  out <- list(scans = scans, hc = hc, n.clust = n.clust,
              threshold.level = threshold.level,
              threshold.lod = threshold.lod,
              comp.lod = comp.lod,
              support.lod = support.lod,
              summary.plot = summary.plot,
              main = main, heatmap = heatmap,
              ylab = match.arg(ylab),
              sex = match.arg(sex),
              scan.type = match.arg(scan.type))
  class(out) <- "multtrait"
  out
}
###########################################################  
summary.multtrait <- function(object, n.clust = object$n.clust, ...)
{
  ## Pass through to summary.aug.scanone.
  summary(object$scans, hc = object$hc,
          n.clust = n.clust, threshold.level = object$threshold.level,
          threshold.lod = object$threshold.lod,
          comp.lod = object$comp.lod, heatmap = object$heatmap)
}
###########################################################  
print.multtrait <- function(x, ...) print(summary(x, ...))
###########################################################  
plot.multtrait <- function(x,
                           main = x$main,
                           heatmap = x$heatmap,
                           n.clust = x$n.clust,
                           rescale = TRUE,
                           use.cM = FALSE,
                           summary.plot = x$summary.plot,
                           ylab = x$ylab,
                           ...)
{
  ## Find traits that pass threshold.
  if(sum(threshold.pass(x$scans[, -(1:2)],
                        x$threshold.lod,
                        x$scans$chr)) <= 1)
    heatmap <- FALSE
  
  if(heatmap) {
    if(x$scan.type %in% c("MOM","PAT","BIM")) {
      rescale <- "peaks"
      lod.name <- paste("Count of", x$scan.type, "peaks")
    }
    else {
      lod.name <- ifelse(is.null(x$support.lod),
                         "Composite LOD",
                         paste("Count of", x$support.lod,
                               "LOD Support Intervals (blue = peaks)"))
    }
    ## Pass through to plot.aug.scanone.
    plot(x$scans, threshold.lod = x$threshold.lod,
         main = main, hc = x$hc, n.clust = n.clust,
         support.lod = x$support.lod, comp.lod = x$comp.lod,
         rescale = rescale, use.cM = use.cM, ylab = ylab,
         lod.name = x$scan.type, ...)

    ## Summary plot if requested.
    if(summary.plot != "none") {
      if(is.null(x$comp.lod)) {
        ## Pass through to summary.aug.scanone and plot.summary.aug.scanone.
        plot(summary(x$scans, threshold.lod = 0,
                     long = TRUE, ...),
             main = main, by = summary.plot)
      }
      else {
        myplot.scanone(x$comp.lod, threshold.lod = 0, lod.name = lod.name,
                       main = main, use.cM = use.cM, xaxs = "i", ...)
      }
    }
  }
  else { ## Individual scans.
    myplot.scanone(x$scans, threshold.lod = x$threshold.lod,
                   main = main, use.cM = use.cM, lod.name = x$scan.type, ...)
  }
  invisible()
}
