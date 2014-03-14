test.logit = function() {
  data(longley)
  
  yr = 1955
  bw = 10
  wt = (1 - ((longley$Year - yr)/bw)**2)**2
  
  Y = rbinom(16, prob=0.5, size=1)
  X = longley[,c(5)]
  
  n = nrow(X)
  p = ncol(X)
  
  X.aug = as.matrix(cbind(X, X * (longley$Year-yr)))
  group = c(1,1)
  
  data = list(x=X.aug, y=Y)
  model = SGL(data, group, wt, alpha=0, min.frac=0.0001, nlam=100, standardize=FALSE, adaptive=TRUE, type='logit')
}