algo_leverage <- function(x, y, subset_size,
                          method = 'both',
                          num_sample = 500){
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

