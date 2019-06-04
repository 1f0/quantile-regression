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


    print(NROW(x1))
    print(nrow(as.matrix(x1)))
    print(class(x1))
    print(typeof(x1))
    print(length(taumat))
    print(length(k))


    print(rq.fit.hogg(as.matrix(x1), y[samp], weights=k, taus=taumat))
    stop()

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

