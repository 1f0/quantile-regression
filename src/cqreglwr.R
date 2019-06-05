mycqr <-
function (x, y, taus = c(.1,.3,.5), weights = NULL,  
    R= NULL, r = NULL, beta = 0.99995, eps = 1e-06) 
{
    if (is.null(weights))
        stop("cannot handle NULL weights")
    kern_weights = weights
    ntau = length(taus)
    weights = rep(1 / ntau, ntau)
    if (length(kern_weights) != length(y))
        stop("kern weights and y differ in length")

    n <- length(y)
    n2 <- NROW(R)
    m <- length(taus)
    p <- ncol(x)+m
    if (n != nrow(x)) 
        stop("x and y don't match n")
    if (m != length(weights)) 
        stop("taus and weights differ in length")
    if (any(taus < eps) || any(taus > 1 - eps)) 
        stop("taus outside (0,1)")
    W <- diag(weights)
    if(m == 1) W <- weights
    x <- as.matrix(x)
    X <- cbind(kronecker(W,rep(1,n)),kronecker(weights,x))
    y <- kronecker(weights,y)
    rhs <- c(weights*(1 - taus)*n, sum(weights*(1-taus)) * apply(x, 2, sum))
    if(n2!=length(r))
    stop("R and r of incompatible dimension")
    if(!is.null(R))
    if(ncol(R)!=p)
        stop("R and X of incompatible dimension")
    d <- rep(1, m*n)
    u <- rep(1, m*n)
    if(length(r)){
       wn1 <- rep(0, 10 * m*n)
       wn1[1:(m*n)] <- .5
       wn2 <- rep(0,6*n2)
       wn2[1:n2] <- 1 
       z <- .Fortran("rqfnc", as.integer(m*n), as.integer(n2), as.integer(p), 
           a1 = as.double(t(as.matrix(X))), c1 = as.double(-y), 
           a2 = as.double(t(as.matrix(R))), c2 = as.double(-r), 
           rhs = as.double(rhs), d1 = double(m*n), d2 = double(n2), 
           as.double(u), beta = as.double(beta), eps = as.double(eps), 
           wn1 = as.double(wn1), wn2 = as.double(wn2), wp = double((p + 3) * p), 
       it.count = integer(3), info = integer(1))
    }
    else{
        wn <- rep(0, 10 * m*n)
        wn[1:(m*n)] <- .5
        z <- .Fortran("rqfnb", as.integer(m*n), as.integer(p), a = as.double(t(as.matrix(X))), 
        c = as.double(-y), rhs = as.double(rhs), d = as.double(d), as.double(u), 
        beta = as.double(beta), eps = as.double(eps), wn = as.double(wn), 
        wp = double((p + 3) * p), it.count = integer(2), info = integer(1))
    }
    # if (z$info != 0) 
    #     warning(paste("Info = ", z$info, "in stepy: singular design: iterations ", z$it.count[1]))

    coefficients <- -z$wp[1:p]
    if(any(is.na(coefficients)))stop("NA coefs:  infeasible problem?")
    list(coefficients = coefficients, nit = z$it.count, flag = z$info)
}

# copy from quantreg
cqreglwr <- function(form,taumat=c(.10,.25,.50,.75,.90),window=.25,bandwidth=0,kern="tcub",distance="Mahal",target=NULL,data=NULL) {
  ntau = length(taumat)

  mat <- model.frame(form,data=data)
  y <- mat[,1]
  n = length(y)
  xmat <- as.matrix(model.matrix(form,data=data)[,-1])
  nk = ncol(xmat)
  if (nk==1) {vxmat <- var(xmat) }
  if (nk==2) {
    vxmat <- cov(xmat) 
    if (distance=="Euclid"|distance=="E") {vxmat <- diag(diag(vxmat)) }
  }

  if (kern=="rect")  { wgt <- function(psi) {ifelse(abs(psi)>=0,1,0) } }
  if (kern=="tria")  { wgt <- function(psi) {1 - abs(psi) } }
  if (kern=="epan")  { wgt <- function(psi) { 1-psi^2 } }
  if (kern=="bisq")  { wgt <- function(psi) { (1-psi^2)^2 } }
  if (kern=="tcub")  { wgt <- function(psi) { (1 - abs(psi)^3)^3 } }
  if (kern=="trwt")  { wgt <- function(psi) { (1 - psi^2)^3 } }
  if (kern=="gauss") { wgt <- function(psi) { exp(-((2.5*psi)^2)/2) } }

  if (bandwidth>0) {window = 0}

  if (identical(target,NULL)){
    target <- maketarget(form,window=window,bandwidth=bandwidth,kern="tcub",data=data)$target
  }
  alldata = FALSE
  if (identical(target,"alldata")){
    target <- xmat
    alldata = TRUE
  }
  if (bandwidth>0){window = 0}
  target <- as.matrix(target)
  nt = nrow(target)

  
  if (distance=="Latlong"|distance=="L") {
    tvect <- attr(terms(form),"term.labels")
    if (substr(tvect[1],1,2)=="la"|substr(tvect[1],1,2)=="La"|substr(tvect[1],1,2)=="LA") {
      la  <- 2*pi*xmat[,1]/360 
      la1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="la"|substr(tvect[2],1,2)=="La"|substr(tvect[2],1,2)=="LA") {
      la  <- 2*pi*xmat[,2]/360 
      la1 <- 2*pi*target[,2]/360 
    }
    if (substr(tvect[1],1,2)=="lo"|substr(tvect[1],1,2)=="Lo"|substr(tvect[1],1,2)=="LO") {
      lo  <- 2*pi*xmat[,1]/360 
      lo1 <- 2*pi*target[,1]/360
    }
    if (substr(tvect[2],1,2)=="lo"|substr(tvect[2],1,2)=="Lo"|substr(tvect[2],1,2)=="LO") {
      lo  <- 2*pi*xmat[,2]/360 
      lo1 <- 2*pi*target[,2]/360
    }
  }

  ytarget     <- array(0,dim=c(nt,ntau))
  ytarget.se  <- array(0,dim=c(nt,ntau))
  dtarget1    <- array(0,dim=c(nt,ntau))
  dtarget1.se <- array(0,dim=c(nt,ntau))
  dtarget2    <- array(0,dim=c(nt,ntau))
  dtarget2.se <- array(0,dim=c(nt,ntau))
  yhat     <- array(0,dim=c(n,ntau))
  yhat.se  <- array(0,dim=c(n,ntau))
  dhat1    <- array(0,dim=c(n,ntau))
  dhat1.se <- array(0,dim=c(n,ntau))
  dhat2    <- array(0,dim=c(n,ntau))
  dhat2.se <- array(0,dim=c(n,ntau))

  for (i in seq(1:nt)) {
    if (distance!="Latlong"&distance!="L")  dist <- sqrt(mahalanobis(xmat, target[i,], vxmat))
    if (distance=="Latlong"|distance=="L") {
      dist <- pmin(sin(la)*sin(la1[i]) + cos(la)*cos(la1[i])*cos(lo1[i]-lo),  1)
      dist <- acos(dist)*3958
    }
    if (window>0) {h = quantile(dist,window) }
    if (bandwidth>0) {h = bandwidth}
    samp <- dist<=h
    if (kern=="gauss") {samp <- dist<=max(dist)}
    
    x1 <- xmat[samp,1]-target[i,1]
    if (nk==2) {x2 <- xmat[samp,2]-target[i,2] }
    k <- wgt(dist[samp]/h)

    # hacking replace
    rq.fit.hogg <- mycqr
    print(summary(rq.fit.hogg(as.matrix(x1), y[samp], weights=k, taus=taumat), covariance=TRUE))
    stop("debug")

    if (nk==1) {fit <- summary(rq.fit.hogg(y[samp]~x1, weights=k, taus=taumat), covariance=TRUE) }
    if (nk==2) {fit <- summary(rq.fit.hogg(y[samp]~x1+x2, weights=k, taus=taumat), covariance=TRUE) }
    ytarget[i] = fit$coef[1]
    ytarget.se[i] = sqrt(fit$cov[1,1])
    dtarget1[i] = fit$coef[2]
    dtarget1.se[i] = sqrt(fit$cov[2,2])
    if (nk==2) {
      dtarget2[i] = fit$coef[3]
      dtarget2.se[i] = sqrt(fit$cov[3,3])
    }
  }

  if (alldata==FALSE) {
    yhat <- smooth12(target,ytarget,xmat)
    dhat1 <- smooth12(target,dtarget1,xmat)
    dhat2 <- 0
    if (nk==2){dhat2 <- smooth12(target,dtarget2,xmat)}
    yhat.se <- smooth12(target,ytarget.se,xmat)
    dhat1.se <- smooth12(target,ytarget.se,xmat)
    dhat2.se <- 0
    if (nk==2){dhat2.se <- smooth12(target,dtarget2.se,xmat)}
  }

  if (alldata==TRUE) {
    yhat <- ytarget
    dhat1 <- dtarget1
    dhat2 <- dtarget2
    yhat.se  <- ytarget.se
    dhat1.se <- dtarget1.se
    dhat2.se <- dtarget2.se
  }

  out <- list(target,ytarget,dtarget1,dtarget2,ytarget.se,dtarget1.se,dtarget2.se,
                    yhat,dhat1,dhat2,yhat.se,dhat1.se,dhat2.se)
  names(out) <- c("target","ytarget","dtarget1","dtarget2","ytarget.se","dtarget1.se","dtarget2.se",
                  "yhat","dhat1","dhat2","yhat.se","dhat1.se","dhat2.se")
  return(out)    
}
# copy from McSpatial

