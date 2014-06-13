SGL <- function(data, index, weights=NULL, type="linear", maxit=1000, thresh=0.001, min.frac=0.1, nlam=20, delta=2, gamma=0.8, standardize=TRUE, verbose=FALSE, step=1, reset=10, alpha=0.95, lambdas=NULL, adaptive=TRUE, unpenalized=NULL){
  X.transform <- NULL
  if (is.null(weights)) {weights = rep(1,nrow(data$x))}
  
  if (adaptive) {
    X = data$x
    Y = data$y
    
    #Center the X matrix:
    meanx = colMeans(X)
    X.centered = sweep(X, 2, meanx, '-')
    
    #Scale the X matrix to unit norm:
    normx = apply(X.centered, 2, function(x) {sqrt(sum(x**2))})
    X.normalized = sweep(X.centered, 2, 1/normx, '*')

    #Center the Y matrix:
    meany = mean(Y)
    Y.centered = Y - meany

    #Scale the X matrix adaptively for the group lasso:
    adamodel = lsfit(y=as.matrix(Y.centered), x=as.matrix(X.normalized), intercept=FALSE, wt=weights)
    s2 = sum(weights * adamodel$residuals**2)/sum(weights)
    adapt = adamodel$coef
    
    groups = unique(index)
    n.g = length(groups)
    adaweights = rep(1, n.g)
    for (i in 1:n.g) {
      g = groups[i]
      indx = which(groups == g)
      adaweights[indx] = 1 / sqrt(sum(adapt[indx]**2))**delta
    }
    
    #Indicate the groups whose coefficients are unpenalized:
    for (g in unpenalized) {
      indx = which(groups == g)
      adaweights[indx] = 0
    }
    
    data = list(x=X.normalized, y=Y.centered)
  }

  if (standardize) {
    X <- data$x
    means <- apply(X,2,mean)
    X <- t(t(X) - means)
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    X <- t(t(X) / var)
    data$x <- X
    X.transform <- list(X.scale=var, X.means=means)
  }

  if (type == "linear") {
    if (standardize) {
      intercept <- mean(data$y)
      data$y <- data$y - intercept
    }
    
    Sol <- oneDim(data, index, weights, adaweights=adaweights, thresh, inner.iter=maxit, outer.iter=maxit, outer.thresh=thresh, min.frac=min.frac, nlam=nlam, lambdas=lambdas, gamma=gamma, verbose=verbose, step=step, reset=reset, alpha=alpha)
    
    res = list()
    if (adaptive) {
      beta = sweep(Sol$beta, 1, 1/normx, '*')
      #beta = sweep(Sol$beta, 1, adaweights/normx, '*')
      intercept = as.vector(meany - t(as.matrix(meanx)) %*% beta)
      Sol$beta = beta
  
      res[['fitted']] = fitted = sweep(as.matrix(X) %*% beta, 2, intercept, '+')
      res[['residuals']] = resid = sweep(fitted, 1, Y, '-')
      
      #Calculate the degrees of freedom used in estimating the coefficients.
      #See Wang and Leng, 2008 (Computational Statistics and Data Analysis (52) pg5279), for details 
      not.zip = matrix(0, nrow=0, ncol=length(lambdas))
      group.df = matrix(0, nrow=0, ncol=ncol(beta))
      
      groups = unique(index)
      for (i in 1:length(groups)) {
        indx = which(index == groups[i])
        adaweight = adaweights[i]
        
        #group.df = rbind(group.df, apply(beta, 2, function(b) ifelse(!all(b[indx]==0), 1 + (length(indx)-1) * sqrt(sum(b[indx]**2)) / adaweight, 0)))
        group.df = rbind(group.df, apply(beta, 2, function(b) ifelse(!all(b[indx]==0), 1 + (length(indx)-1) * sqrt(sum(b[indx]**2)) * adaweight, 0)))
      }
      
      #res[['df']] = df = apply(group.df, 2, sum)
      #Naive df (add one for the intercept, which is estimated but does not appear in b)
      res[['df']] = df = drop(apply(beta, 2, function(b) sum(b!=0) + 1))
      res[['BIC']] = apply(resid, 2, function(x) sum(weights * x**2)) / s2 + log(sum(weights))*df
      res[['AIC']] = apply(resid, 2, function(x) sum(weights * x**2)) / s2 + 2*df
      res[['AICc']] = apply(resid, 2, function(x) sum(weights * x**2)) / s2 + 2*df + 2*df*(df+1)/(sum(weights)-df-1)
    }
    
    if (adaptive) {
      Sol <- list(beta=Sol$beta, lambdas=Sol$lambdas, type=type, intercept=intercept, X.transform=X.transform, LS.coefs=adapt, adaweights=adaweights, weights=weights, results=res)
    } else if (standardize) {
      Sol <- list(beta=Sol$beta, lambdas=Sol$lambdas, type=type, intercept=intercept, X.transform=X.transform)
    } else {
      Sol <- list(beta=Sol$beta, lambdas=Sol$lambdas, type=type, X.transform=X.transform)
    }
  }

  if (type == "logit") {
    Sol <- oneDimLogit(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)
    Sol <- list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, intercept = Sol$intercept, X.transform = X.transform)
  }

  if (type == "cox") {
    Sol <- oneDimCox(data, index, thresh = thresh, inner.iter = maxit, outer.iter = maxit, outer.thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)
    Sol = list(beta = Sol$beta, lambdas = Sol$lambdas, type = type, X.transform = X.transform)
  }
  
  class(Sol) = "SGL"
  return(Sol)
}


cvSGL <- function(data, index = rep(1,ncol(data$x)), type = "linear", maxit = 1000, thresh = 0.001, min.frac = 0.05, nlam = 20, gamma = 0.8, nfold = 10, standardize = TRUE, verbose = FALSE, step = 1, reset = 10, alpha = 0.95, lambdas = NULL){

  if(standardize == TRUE){
    X <- data$x
    means <- apply(X,2,mean)
    X <- t(t(X) - means)
    var <- apply(X,2,function(x)(sqrt(sum(x^2))))
    X <- t(t(X) / var)
    data$x <- X
  }

  if(type == "linear"){
   if(standardize == TRUE){
      intercept <- mean(data$y)
      data$y <- data$y - intercept
    }   
    Sol <- linCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, reset = reset, alpha = alpha)

    if(standardize == TRUE){
       Sol$fit = list(beta = Sol$fit$beta, lambdas = Sol$fit$lambdas, intercept = intercept, step = step)
    }
  }
  if(type == "logit"){
    Sol <- logitCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)

  }
  if(type == "cox"){
    Sol <- coxCrossVal(data, index, nfold = nfold, maxit = maxit, thresh = thresh, min.frac = min.frac, nlam = nlam, lambdas = lambdas, gamma = gamma, verbose = verbose, step = step, alpha = alpha, reset = reset)

  }

  Sol = list(fit = Sol$fit, lldiff = Sol$lldiff, lambdas = Sol$lambdas, type = type, llSD = Sol$llSD)

class(Sol) = "cv.SGL"

    return(Sol)
}
