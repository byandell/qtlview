###########################################################################3
my.length.X <- function(cross)
  dim(cross$geno$X$prob)[2]
###########################################################################3
my.add.X <- function(cross, newdata, transpose = FALSE)
{
  ## QTLBIM calcs do not (yet) use X chr. Need to add it.
  X.names <- dimnames(cross$geno$X$prob)[[2]]
  if(transpose) {
    tmp <- cbind(newdata, matrix(0, nrow(newdata), length(X.names)))
    dimnames(tmp) <- list(NULL,
                          c(dimnames(newdata)[[2]], X.names))
  }
  else {
    tmp <- rbind(newdata, matrix(0, length(X.names), ncol(newdata)))
    dimnames(tmp) <- list(NULL, dimnames(newdata)[[2]])
  }
  tmp
}
                  
