listof.complod <- function(file.groups, dir = NULL,
                        scan.type = "LOD", ...)
{
  group.names <- names(file.groups)
  file.groups <- file.path(dir, file.groups)
  n.groups <- length(file.groups)

  scan.type <- array(scan.type, n.groups)
  names(scan.type) <- group.names
  
  scan.groups <- list()
  for(i in seq(n.groups)) {
    cat("group", group.names[i], "\n")
    tmp <- multtrait(file=file.groups[i],
                     scan.type = scan.type[i],
                     cross.name = "B6BTBR07",
                     drop.chr = FALSE, ...)
    ## Turn count into percent.
    for(j in 3:4)
      tmp$comp.lod[[j]] <- 100 * tmp$comp.lod[[j]] / (ncol(tmp$scans) - 2)

    if(i == 1)
      support.lod <- tmp$support.lod

    scan.groups[[i]] <- tmp$comp.lod
  }
  names(scan.groups) <- group.names

  if(!any(scan.type %in% c("LPD","LOD")))
    support.lod <- NULL
  attr(scan.groups, "support.lod") <- support.lod
  attr(scan.groups, "scan.type") <- scan.type
  
  class(scan.groups) <- c("listof.complod", "list")
  scan.groups
}
plot.listof.complod <- function(x,
                                lodcolumn = 1,
                                ylim = ylims,
                                ylab = paste("Percent of",
                                  attr(x, "support.lod"),
                                  attr(x, "scan.type")[1],
                                  "Support Intervals"),
                                col = seq(n.x),
                                lty = seq(n.x),
                                pre.main = "",
                                main = paste(pre.main, mains),
                                ...)
{
  tmp <- ncol(x[[1]]) - 2
  if(lodcolumn > tmp | lodcolumn < 1)
    stop(paste("lodcolumn must be between 1 and", tmp))
     
  n.x <- length(x)

  ## Colors.
  cols <- palette()
  if(is.numeric(col))
    col <- cols[1 + ((col - 1) %% length(cols))]
  col <- array(col, n.x)

  ## Line types.
  ltys <- c("solid","dashed","dotted","dotdash","longdash","twodash")
  if(is.numeric(lty))
    lty <- ltys[1 + ((lty - 1) %% length(ltys))]
  
  mains <- paste(col, names(x), sep = "=", collapse = ", ")
  
  ylims <- 0
  for(i in seq(n.x)) 
    ylims <- range(ylims, x[[i]][[lodcolumn + 2]])
  
  for(i in seq(n.x))
    plot.scanone(x[[i]], lodcolumn = lodcolumn, ylim = ylim, ylab = ylab,
                 col = col[i], lty = lty[i], add = (i > 1), main = main, ...)
}

