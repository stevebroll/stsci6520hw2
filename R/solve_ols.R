#' Iterative solver for Linear Systems of Equations
#'
#' solve_ols can use Gauss-Seidel , Jacobi(sequential), or Jacobi
#' (parallel) to solve a system of equations Ax = b given A and b.
#'
#' @param A Matrix with n rows and p columns
#' @param b Vector of length n
#' @param method Iterative method for solving Ax = b, either 'gs' for
#' Gauss-Seidel (sequential), 'jacobi' for Jacobi (sequential), or
#' 'parallel' for parallel Jacobi.
#' @param iter Number of iterations for the solver.
#' @param ncores Optionally set number of cores. If left empty, solve_ols will
#' get the number of cores on the user's system and subtract 1.
#'
#' @return Numeric vector x, solution to Ax = b.
#' @export
#'
#' @examples
#' A = diag(1, nrow = 100)
#' A[abs(row(A) -col(A)) == 1] = -1
#' b = A%*%(rep(c(1,0), 50))
#' solve_ols(A,b)
#' solve_ols(A,b, method = 'jacobi', iter = 1000, ncores = 6)
#'

solve_ols <- function(A, b, method = 'gs', iter = 5000, ncores = NULL){
  n = ncol(A)
  x = rep(0,n)
  D = diag(diag(A))
  L = U = A
  U[upper.tri(U, diag = T)] = 0
  L[lower.tri(L, diag = T)] = 0

  if(method == 'gs'){
    LDinv = solve(L + D)
    for(i in 1:iter){
      x = LDinv %*% (b - U %*% x)
    }
  }
  else if(method == 'jacobi'){
    LU = L + U
    Dinv = solve(D)
    for(i in 1:iter){
      x = Dinv%*%(b- (LU)%*%x)
    }
  }
  else if(method == 'parallel'){
    if(is.null(ncores)){
      ncores = as.numeric(Sys.getenv("NUMBER_OF_PROCESSORS"))
    }

    cl = makeCluster(ncores)
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
