test.logit = function() {
  data(longley)
  
  yr = 1955
  bw = 10
  wt = 1 - ((longley$Year - yr)/bw)**2
  w2 = rep(1,16)
  
  X = longley[,c(3,4)]
  Y = rbinom(16, prob=longley$GNP/(1.5*max(longley$GNP)), size=1)
  
  n = nrow(X)
  p = ncol(X)
  
  X.aug = as.matrix(cbind(X, X * (longley$Year-yr), longley$Year-yr))
  group = c(1,2,1,2,0)
  
  data = list(x=X.aug, y=Y)
  model = SGL(data, group, wt, alpha=0, min.frac=0.0001, nlam=100, standardize=FALSE, adaptive=TRUE, type='logit', unpenalized=c(0))

  m2 = SGL(data, group, w2, alpha=0, min.frac=0.0001, nlam=100, standardize=FALSE, adaptive=TRUE, type='logit')

  model.glm = glm(Y~X.aug, family='binomial', weights=wt)
  m2.glm = glm(Y~X.aug, family='binomial', weights=w2)
}