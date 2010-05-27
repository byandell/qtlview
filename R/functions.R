mystop <- function(x, ...)
{
  class(x) <- c("mystop", class(x))
  x
}
summary.mystop <- function(object, ...)
  cat(object)
plot.mystop <- function(x, ...)
  cat(x)
###########################################################  
normal.trans <- function(x) {
  x <- rank(x, na.last = "keep")
  qnorm(x / (1 + sum(!is.na(x))))
}
################################################################
find.chr <- function(chr, chr.names)
{
  ## Double-check chr.
  if(length(chr)) {
    if(chr[1] != "")
      chr <- as.character(chr[!is.na(match(chr, chr.names))])
    else
      chr <- chr.names
  }
  else
    chr <- chr.names
  chr
}
################################################################
mytrait <- function(traitnames, pheno, warn = TRUE, exact = FALSE)
{
  pheno.names <- dimnames(pheno)[[2]]
  if(missing(traitnames))
    return(pheno.names)
  
  ## Expand trait names if not complete.
  tmp <- pheno.names[if(length(traitnames) > 1 | exact)
                      {
                        pmatch(tolower(make.names(traitnames)),
                               tolower(make.names(pheno.names)),
                               duplicates.ok = TRUE)
                      }
                      else
                      {
                        grep(tolower(make.names(traitnames)),
                             tolower(make.names(pheno.names)))
                      }
                      ]

  tmp1 <- is.na(tmp)
  if(any(tmp1)) {
    ## Try prepending "wk" or "X".
    ## First remove "wk" from name.
    tmp3 <- sapply(strsplit(traitnames[tmp1], "wk"), paste, collapse = "")
    tmp2 <- pmatch(tolower(make.names(paste("wk", tmp3, sep = ""))),
                   tolower(pheno.names), nomatch = 0)
    tmp2 <- pmax(pmatch(tolower(make.names(paste("X", tmp3, sep = ""))),
                        tolower(pheno.names), nomatch = 0),
                 tmp2)
    tmp2[tmp2 == 0] <- NA
    tmp2 <- pheno.names[tmp2]
    tmp[tmp1] <- tmp2
    
    if(any(is.na(tmp)) & warn) {
      print(traitnames)
      print(tmp)
      warning("traitnames do not match names(pheno)")
    }
  }
  tmp
}
