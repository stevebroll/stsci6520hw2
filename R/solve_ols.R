solve_ols <- function(A, b, method = 'gs', iter = 5000, ncores = NULL){
  n = nrow(A)
  x = rep(0,n)
  D = diag(diag(A))
  L = U = A
  U[upper.tri(U, diag = T)] = 0
  L[lower.tri(L, diag = T)] = 0

  if(method = 'gs'){
    LDinv = solve(L + D)
    for(i in 1:iter){
      x = LDinv %*% (b - U %*% x)
    }
  }
  else if(method = 'jacobi'){
    LU = L + U
    Dinv = solve(D)
    for(i in 1:iter){
      x = Dinv%*%(b- (LU)%*%x)
    }
  }
  else if(method = 'parallel'){
    require(doParallel)
    if(is.null(ncores)){
      ncores = as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
    }

    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    update = function(i,x){
      x[i] = (b[i] - A[i,-i]%*% x[-i])/(A[i,i])
      return(x)
    }
    for(j in 1:iter){
      x = as.numeric(foreach(i = 1:n) %dopar% {update(i,x)[i]})
      stopCluster(cl)
    }
  }
  else{
    stop("Please choose method = 'gs', 'jacobi', or 'parallel'")
  }

  return(x)
}
