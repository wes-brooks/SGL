test = function() {
    data(longley)
    
    yr = 1955
    bw = 10
    wt = (1 - ((longley$Year - yr)/bw)**2)**2
    
    Y = longley$GNP
    X = longley[,c(3,4,5,7)]
    
    n = nrow(X)
    p = ncol(X)
    
    X.aug = as.matrix(cbind(matrix(longley$Year-yr), X, X * (longley$Year-yr)))
    group = c(0,1,2,3,4,1,2,3,4)
    
    data = list(x=X.aug, y=Y)
    model = SGL(data, group, wt, alpha=0, min.frac=0.0001, nlam=100, standardize=FALSE, adaptive=TRUE, unpenalized=c(0))
}
