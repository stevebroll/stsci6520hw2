#' Algorithmic Leveraging
#'
#' algo_leverage fits linear regression models on subsets of data sampled by
#' both uniform and leverage sampling
#'
#' @param x Predictor matrix x, with n rows and p columns
#' @param y Response vector y, with n elements
#' @param subset_size Size of subset(s) on which regression model(s) will be fit
#' @param num_sample Number of subsets to take when method = 'both'
#' @param method One of 'both', 'uniform', or 'leverage'. 'both' will produce
#' summary outputs for both uniform and leverage sampling over num_sample
#' subsets, while 'uniform' and 'leverage' output one draw of estimated Betas
#' for uniform and leverage sampling, respectively.
#'
#' @return When method = 'both' and X has 1 column, boxplots will be produced
#' showing the distribution of betas from models fit on samples drawn using both
#' uniform and leverage sampling. When method = 'both' and X has more than 1
#' column, a line graph will compare the average values of each Beta for both
#' uniform and leverage sampling. When method = 'uniform' or 'leverage', the
#' model will be fit for one sample using uniform or leverage sampling,
#' respectively, and the Beta values will be returned.
#' @export
#'
#' @examples
algo_leverage <- function(x, y, subset_size,
                          num_sample = 500,
                          method = 'both'
                          ){
  x = as.matrix(x)
  n = nrow(x); p = ncol(x)
  r = subset_size
  if(r >= n){
    stop("r must be smaller than number of rows in X")
  }

  if(method == 'both'){

    hmat = x %*% MASS::ginv((t(x)%*% x))%*%t(x)
    hvec = diag(hmat)/sum(diag(hmat))
    unif = MASS::ginv(diag(1/n, nrow = n))
    levg = MASS::ginv(diag(hvec, nrow = n))

    uniout = levout = matrix(NA, num_sample, p)

  for(i in 1:num_sample){
    suni = sample(1:n, size = r, replace = T, prob = rep(1/n, n))
    slvg = sample(1:n, size = r, replace = T, prob = hvec)
    xuni = x[suni,]; yuni = y[suni]
    xlev = x[slvg,]; ylev = y[slvg]

    umat = unif[suni, suni]
    lmat = levg[slvg, slvg]

    uniout[i,] = MASS::ginv(t(xuni)%*%umat%*%xuni) %*% t(xuni) %*% umat %*% yuni

    levout[i,] = MASS::ginv(t(xlev)%*%lmat%*%xlev) %*% t(xlev) %*% lmat %*% ylev
  }
  cat(paste('Mean Betas for ', num_sample,
              ' Subsamples: \n'))
  print(data.frame(Uniform = colMeans(uniout), Leverage = colMeans(levout)))
  cat('\n')


   if(p == 1){
      boxplot(uniout[,1], levout[,1], main = paste(num_sample, "Draws of Beta
      for Uniform and Levarage Subsamples"),names = c("Uniform", "Leverage") )
    }
   else if(p>1){
    plot(colMeans(uniout), col = 'red', main = "Uniform vs Leverage Subsampling",
         xlab = 'Beta', ylab = 'Coefficient', type = 'l')
    lines(colMeans(levout), col = 'blue')
    legend("topleft", c('Uniform', 'Leverage'), inset = .04,
         col = c('red', 'blue'), lty = 1)
    }
  }
  else if(method == 'uniform'){

    unif = MASS::ginv(diag(1/n, nrow = n))
    suni = sample(1:n, size = r, replace = T, prob = rep(1/n, n))
    xuni = x[suni,]; yuni = y[suni]
    umat = unif[suni, suni]
    return(data.frame(Beta = MASS::ginv(t(xuni)%*%umat%*%xuni) %*% t(xuni) %*% umat %*% yuni))

  }
  else if(method == 'leverage'){
    hmat = x %*% MASS::ginv((t(x)%*% x))%*%t(x)
    hvec = diag(hmat)/sum(diag(hmat))
    levg = MASS::ginv(diag(hvec, nrow = n))
    slvg = sample(1:n, size = r, replace = T, prob = hvec)
    xlev = x[slvg,]; ylev = y[slvg]
    lmat = levg[slvg, slvg]
    cat(paste('Single Draw of Betas: \n'))
    return(data.frame(Beta = MASS::ginv(t(xlev)%*%lmat%*%xlev) %*% t(xlev) %*% lmat %*% ylev))
  }
  else{
    stop("Please choose method = 'both', 'leverage', or 'uniform'")
  }

}

