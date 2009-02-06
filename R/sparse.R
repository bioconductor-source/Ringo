## functions for retrieving the indices of non-zero elements in sparse
##  matrices; originally used the matrix.csr class from SparseM, but now
##  switched to dgCMatrix from package Matrix

setMethod("nonzero", signature(x="matrix"), function(x){
    return(which(x != 0, arr.ind=TRUE))
})

setMethod("nonzero", signature(x="matrix.csr"), function(x){
  ## function to get a two-column matrix containing the indices of the
  ### non-zero elements in a "matrix.csr" class matrix
  if (all(x@ra==0))
      return(matrix(0, nrow=0, ncol=2, dimnames=
                    list(character(0), c("row","col"))))
  res <- cbind(rep(seq(dim(x)[1]),diff(x@ia)), # row indices
               x@ja )# column indices directly saved in matrix.csr
  colnames(res) <- c("row","col")
  ## remove zero elements
  res <- res[x@ra != 0,,drop=FALSE]
  return(res)
})

setMethod("nonzero", signature(x="dgCMatrix"), function(x){
    ## function to get a two-column matrix containing the indices of the
    ### non-zero elements in a "dgCMatrix" class matrix
    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0))
        return(matrix(0, nrow=0, ncol=2,
               dimnames=list(character(0), c("row","col"))))
    res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
})
