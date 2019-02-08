
library(optimx)
library(SpatialExtremes)
library(numDeriv)
library(MASS)
dyn.load(paste("spt-gamma",.Platform$dynlib.ext,sep=""))




data2frechet.semipar<-function(data,prob=0.95) {
	
	nsites<-dim(data)[2]
	ntimes<-dim(data)[1]
	newdata<-matrix(NA,ntimes,nsites)
	for (i in 1: nsites) {
			# only positive value are considered in the threshold calculation
      newdata[,i]<-gpd2frechet.semipar(data[,i],prob = prob)
	}
	return(newdata)
}

data2gpdcensored.single<-function(y,prob=0.95, kappa,method ="Nelder-Mead",type = 8) {
  sel.na<-is.na(y)
  y.tmp<-y[!sel.na]
  threshold<-quantile(y.tmp,prob = prob,type=type)
  marginal.param <-as.numeric(gpdmle(y.tmp,threshold=threshold,method = method))
  sel<- y.tmp > threshold

  y.tmp[!sel]<-0
 # we consider only the excedancees
  z <- (kappa + 1.0)*
    ((1.0+marginal.param[2]*(y.tmp[sel]-threshold)
      /marginal.param[1])^(1/marginal.param[2])-1) 
  y.tmp[sel]<-z 
  
  y[!sel.na]<-y.tmp
  return(y)
}



data2gpdcensored<-function(data,prob=0.95, kappa,method ="BFGS",type = 8) {
  
  newdata<-apply(data,2,data2gpdcensored.single, prob=prob,kappa=kappa,type=type)
  
  return(newdata)
}


data2gpdcensored.integral<-function(data,prob=0.95, kappa) {
  
  newdata<-apply(data,2,data2gpdcensored.integral.single, prob=prob,kappa=kappa)
  return(newdata)
}

data2gpdcensored.integral.single<-function(y,prob=0.95, kappa) {
  
    n<-length(y)  
    x<-rep(0,length(y))
    
    rank.data<-rank(y,ties.method = "random")
    sel<-rank.data > prob * n
    z<-(rank.data[sel] - prob * n)/(sum(sel)+1)
    x[sel]<-(kappa + 1.0)*((1-z)^(-1)-1)
  
  return(x)
}
gpd2frechet.semipar<-function(x,prob=0.95) {
	
  sel<-!is.na(x)
	z<-x[sel]
	threshold<-quantile(z,prob = prob)
	marginal.param <-as.numeric(gpdmle(z,threshold))
	n<-length(z)
	sel.1<-(z > threshold)
	p<-1-mean(sel.1,na.rm=TRUE)
	cdf<-rank(z,ties.method = "random")/(n + 1)
	cdf[sel.1]<-1-(1-p)*(1+marginal.param[2]*(z[sel.1]-threshold)/marginal.param[1])^(-1/marginal.param[2])
	z<- -1/log(cdf)
    	x[sel]<-z
        return(x)
}

data2frechet.nonpar<-function(data) {
  
  nsites<-dim(data)[2]
  ntimes<-dim(data)[1]
  newdata<-matrix(NA,ntimes,nsites)
  
  for (i in 1: nsites) {
    # only positive value are considered in the threshold calculation
    sel<-!is.na(data[,i])
    z<-data[sel,i]
    n<-length(z)
    cdf<-rank(z,ties.method = "random")/(n + 1)
    z<- -1/log(cdf)
    newdata[,i]<-data[,i]
    newdata[sel,i]<-z  
  }
  return(newdata)
}

# this function evaluate the theoretical chi
chi.spt.gamma<-function(param,v, xcoords, ycoords, tcoords) {
  val<-numeric(1)
  tmp<-.C("chi",v=as.double(v),xcoords=as.double(xcoords), 
          ycoords=as.double(ycoords), 
          tcoords=as.double(tcoords), param = as.double(param),
          val = as.double(val), DUP =FALSE)
  return(tmp$val)	
}

chi.spt.gamma.wrap<-function(theta,v,xcoords, ycoords, tcoords,param, mask)
{
  param[mask]<-theta
  val<-chi.spt.gamma(param=param,v=v,xcoords=xcoords, ycoords=ycoords, tcoords=tcoords)
  return(val)
}




dchi.spt.gamma<-function(param,v,xcoords, ycoords, tcoords,mask)
{
  theta<-param[mask]
  val<-grad(x=theta,func = chi.spt.gamma.wrap,v=v,xcoords=xcoords, ycoords=ycoords, tcoords=tcoords,param=param, mask= mask)
  return(val)
}

# this function calculates the empirical estimates of chi(u) and stacks the results in a vector
# useful for calculating LS estimates of the parameters
emp.f<-function(alldata,coords,prob=0.95,maxlag=4,maxdist=NULL,nbreaks=13, angles =c(0,pi/4,pi/2,3*pi/4),tmaxlag=10)
{
  emp.all<-NULL
  d.all<-NULL
  n<-NULL
  lags<-0:maxlag
  emp.chi<-chi.temporal.emp(alldata,u=prob,maxlag =tmaxlag)
  med<-apply(emp.chi[,1:tmaxlag],2,median)
  emp.all<-c(emp.all,med)
  n<-c(n,length(med))
  d.all<-c(d.all,rep(0,length(med)))
  for (lag in lags) {
    for (ang.rad in angles) {
      
      emp.chi<-chi.spatial.emp(data=alldata,coords=coords,u=prob,maxdist=maxdist, lag = lag,
                               omnidirectional=FALSE, ang.rad=ang.rad,tol.rad =pi/8,nbreaks = nbreaks)
      n<-c(n,length(emp.chi[,2]))
      emp.all<-c(emp.all,emp.chi[,2])
      d.all<-c(d.all,emp.chi[,1])
    }
  }
return(list(v=as.numeric(emp.all),d=as.numeric(d.all),angles =angles,lags=lags,n=n,
            maxlag=maxlag,tmaxlag=tmaxlag,
            nsites=ncol(alldata),ntimes=nrow(alldata)))  
}

# this function calculates the theoretical values of chi(u) and stacks the results in a vector
# useful for calculating LS estimates of the parameters
wls.f.wrap<-function(param,chi)
{
  
  if (sum(param[3:6] <0) > 0) {
    val<-Inf
    return(val)
  }
  if (param[5] >= pi/2) {
    val<-Inf
    return(val)
  }
  
  theo.all<-NULL
  n<-c(1,cumsum(chi$n))
  theo<-numeric(chi$n[1])
  
  for (i in 1:chi$n[1]) {
    theo[i]<-chi.spt.gamma(param = param,v=0,
                           xcoords = c(0,chi$d[i]),
                           ycoords = c(0,chi$d[i]),
                           tcoords = c(0,i))  
  }
  theo.all<-c(theo.all,theo)
  k<-2
  for (lag in chi$lags) {
    for (ang.rad in chi$angles) {
      
      theo<-numeric(chi$n[k])
      for (i in 1:chi$n[k]) {
        theo[i]<-chi.spt.gamma(param = param,v=0,
                  xcoords = c(0,chi$d[n[k]+i]*cos(ang.rad)),
                  ycoords = c(0,chi$d[n[k]+i]*sin(ang.rad)),
                  tcoords = c(0,0+lag))  
      }
      theo.all<-c(theo.all,theo)
      k<-k+1
    }
    
  } 
  
  return(theo.all)
}
# this function provide the criterion to be minimized in the LS estimates
wls.f<-function(theta,chi, param,mask)
{
  param[mask]<-theta	
  theo.all<-wls.f.wrap(param,chi)
  val<-sum((theo.all-chi$v)^2)
  return(val)
}



PLneg.spt<-function(theta,y,z, threshold, xcoords, ycoords, tcoords, 
                    delta,param, mask) {
  n<-length(y)
  loglik <- 0
  param[mask]<-theta
  
  v1<-y-threshold # exceedances
  v1[z == 0]<-0
	tmp<-.C("pwl_spatiotemporal",y=as.double(v1),z=as.integer(z),
          n=as.integer(n),theta=as.double(param), loglik=as.double(loglik),u=as.double(threshold), 
          xcoords=as.double(xcoords), ycoords=as.double(ycoords), 
          tcoords=as.double(tcoords), delta = as.double(delta), DUP =FALSE)
  	
  return(-tmp$loglik)	
}

PLneg.spt.parallel<-function(theta,y,z, threshold, xcoords, ycoords, tcoords,
                             delta, param, mask, ncores,nsites,ntimes) {
  
  tmp<-PLneg.spt.parallel.wrap(theta,y,z, threshold, xcoords, ycoords, tcoords,
                             delta, param, mask, ncores,nsites,ntimes)
  return(tmp$val)                           
  }

PLneg.spt.parallel.wrap<-function(theta,y,z, threshold, xcoords, ycoords, tcoords,
                             delta, param, mask, ncores,nsites,ntimes) {
  
  # ncores number of cores
  param[mask]<-theta 
  m<-ntimes %/%ncores
  lag<-delta[2]
  if (lag > m) stop("lag too high")
  ntimes.for.cores<-rep(m,ncores)
  ntimes.for.cores[ncores]<-ntimes.for.cores[ncores]+(ntimes-m*ncores)
  start<-c(0,cumsum(ntimes.for.cores)[-ncores])
  end<-cumsum(ntimes.for.cores)
  end[-ncores]<-end[-ncores]+lag
  results<-numeric(ncores)
  loglik<-rep(0,ncores)
  npairs<-rep(0,ncores)
  results<-foreach(j=1:ncores,.combine=rbind) %dopar%{ 
    
  ind<-(start[j]*nsites+1):(end[j]*nsites)
  v1<-y[ind]-threshold # exceedances
  nobs<-length(v1)
  z1<-z[ind]
  v1[z1 == 0]<-0
    
  val<-.C("pwl_spatiotemporal_parallel",y=as.double(v1),z=as.integer(z1), nobs=as.integer(nobs),
            nsites=as.integer(nsites), nblocks=as.integer(ntimes.for.cores[j]), 
            lag=as.integer(delta[2]), theta=as.double(param), loglik=as.double(loglik[j]),u=as.double(threshold), 
            xcoords=as.double(xcoords[ind]), 
            ycoords=as.double(ycoords[ind]), tcoords=as.double(tcoords[ind]), 
            delta = as.double(delta[1]),count=as.double(npairs[j]),NAOK = TRUE,  DUP =FALSE)	
  
  c(-val$loglik,val$count)
  }  
  
  return(list(val=sum(results[,1]),ntot=sum(results[,2])))	
}




# this function estimates the model

sptgamma.fit<-function(y, xcoords, ycoords, tcoords,  init, delta=c(Inf,Inf), 
                       threshold,
                       maxit = 20000, method ="BFGS", mask=NULL, 
                       marg.param.fixed=TRUE,nsites =NULL, ntimes =NULL, ncores = 1) {

    parnames<-c("xi","sigma","semiax1","semiax2","phi","duration","velocity1","velocity2","alpha","beta","kappa")	
    param<-init
    z <- as.numeric((y > threshold))
    y[z == 0]<-threshold
    res<-list()
    res$threshold<-threshold
    theta<-param[mask]
    res$start<-theta
    names(res$start)<-parnames[mask]  
    ptm <- proc.time()  
    
      parscale<-10^floor(log10(abs(theta+1e-4)))
      a <-optimx(par=theta, fn=PLneg.spt.parallel, method = "Nelder-Mead", control=list(maxit=100,parscale=parscale),y = y, z=z, xcoords = xcoords,
                ycoords= ycoords,tcoords =tcoords,delta = delta, threshold =threshold,
                param=param,mask=mask,ncores = ncores, nsites = nsites, ntimes =ntimes)
      theta<-as.numeric(coef(a))
      a <-optimx(par=theta, fn=PLneg.spt.parallel, method = method, control=list(maxit=1000,parscale=parscale),y = y, z=z, xcoords = xcoords,
                ycoords= ycoords,tcoords =tcoords,delta = delta, threshold =threshold,
                param=param,mask=mask,ncores = ncores, nsites = nsites, ntimes =ntimes)
      
      
      res$negplik <- a$value	
      res$hessian<-attr(a,"details")[1,]$nhatend
      res$grad<-attr(a,"details")[1,]$ngatend
      res$delta<-delta
      param[mask]<-as.numeric(coef(a))
	    res$thetahat<-as.numeric(coef(a))
	    names(res$thetahat)<-parnames[mask]
	    res$param<-param
	    res$mask<-mask
	 	res$convergence<- a$convcode 
	 	res$ncores<-ncores
	 	res$method<-method
	 	res$elapsed<-proc.time() - ptm
	    return(res)	
}



PLneg.spt.parallel.vec.grad<-function(theta,y,z, threshold, xcoords, ycoords, tcoords,
                                 delta, param, mask, ncores,nsites,ntimes, npairs,eps = 1e-4,filename=NULL) {
  
  # ncores number of cores
  param[mask]<-theta
  m<-ntimes %/%ncores
  lag <- delta[2]
  if (lag > m) stop("lag too high")
  ntimes.for.cores<-rep(m,ncores)
  ntimes.for.cores[ncores]<-ntimes.for.cores[ncores]+(ntimes-m*ncores)
  start<-c(0,cumsum(ntimes.for.cores)[-ncores])
  end<-cumsum(ntimes.for.cores)
  end[-ncores]<-end[-ncores]+lag
  
  results.lik<-foreach(j=1:ncores,.combine=c) %dopar%{ 
    ntot<-(npairs+(2*npairs+nsites)*lag)*(((end[j]-start[j])-lag)+lag*(lag-1)/2)  
    loglik<-rep(0,ntot)
    count <- 0
    ind<-(start[j]*nsites+1):(end[j]*nsites)
    v1<-y[ind]-threshold # exceedances
    nobs<-length(v1)
    z1<-z[ind]
    v1[z1 == 0]<-0
    val<-.C("pwl_spatiotemporal_parallel_vec",y=as.double(v1),z=as.integer(z1), nobs=as.integer(nobs),
            nsites=as.integer(nsites), nblocks=as.integer(ntimes.for.cores[j]), 
            lag=as.integer(delta[2]), theta=as.double(param), loglik=as.double(loglik),u=as.double(threshold), 
            xcoords=as.double(xcoords[ind]), 
            ycoords=as.double(ycoords[ind]), tcoords=as.double(tcoords[ind]), 
            delta = as.double(delta[1]),count=as.integer(count),NAOK = TRUE,  DUP =FALSE)	
    
    loglik<--val$loglik[1:(val$count)]
   save(loglik,file = paste(filename,"loglik-tmp-",j,".tmpRdata",sep=""))
    loglik
  }  
  
  
  nparam<-length(theta)

  for (k in 1:nparam) {
    dtheta <- theta
    dtheta[k] <- dtheta[k] + eps
    param[mask]<-dtheta
    results<-foreach(j=1:ncores,.combine=c) %dopar%{ 
      ntot<-(npairs+(2*npairs+nsites)*lag)*(((end[j]-start[j])-lag)+lag*(lag-1)/2)  
      loglik<-rep(0,ntot)
      count<-0
      ind<-(start[j]*nsites+1):(end[j]*nsites)
      v1<-y[ind]-threshold # exceedances
      nobs<-length(v1)
      z1<-z[ind]
      v1[z1 == 0]<-0
      val<-.C("pwl_spatiotemporal_parallel_vec",y=as.double(v1),z=as.integer(z1), nobs=as.integer(nobs),
            nsites=as.integer(nsites), nblocks=as.integer(ntimes.for.cores[j]), 
            lag=as.integer(delta[2]), theta=as.double(param), loglik=as.double(loglik),u=as.double(threshold), 
            xcoords=as.double(xcoords[ind]), 
            ycoords=as.double(ycoords[ind]), tcoords=as.double(tcoords[ind]), 
            delta = as.double(delta[1]),count=as.integer(count),NAOK = TRUE,  DUP =FALSE)	
      loglik<--val$loglik[1:(val$count)]
      save(loglik,file = paste(filename,"loglik-tmp-",k,"-",j,".tmpRdata",sep=""))
      val$count
    } 
  }  
  return(list(negloglik=sum(results.lik),ntot=sum(results)))	
}

PLneg.spt.parallel.vec.grad.tmp<-function(theta,y,z, threshold, xcoords, ycoords, tcoords,
                                 delta, param, mask, ncores,nsites,ntimes, npairs,eps = 1e-4,filename=NULL) {
  
  # ncores number of cores
  param[mask]<-theta
  m<-ntimes %/%ncores
  lag <- delta[2]
  if (lag > m) stop("lag too high")
  ntimes.for.cores<-rep(m,ncores)
  ntimes.for.cores[ncores]<-ntimes.for.cores[ncores]+(ntimes-m*ncores)
  start<-c(0,cumsum(ntimes.for.cores)[-ncores])
  end<-cumsum(ntimes.for.cores)
  end[-ncores]<-end[-ncores]+lag
  
  results.lik<-foreach(j=1:ncores,.combine=c) %dopar%{ 
    ntot<-(npairs+(2*npairs+nsites)*lag)*(((end[j]-start[j])-lag)+lag*(lag-1)/2)  
    loglik<-rep(0,ntot)
    count <- 0
    ind<-(start[j]*nsites+1):(end[j]*nsites)
    v1<-y[ind]-threshold # exceedances
    nobs<-length(v1)
    z1<-z[ind]
    v1[z1 == 0]<-0
    val<-.C("pwl_spatiotemporal_parallel_vec",y=as.double(v1),z=as.integer(z1), nobs=as.integer(nobs),
            nsites=as.integer(nsites), nblocks=as.integer(ntimes.for.cores[j]), 
            lag=as.integer(delta[2]), theta=as.double(param), loglik=as.double(loglik),u=as.double(threshold), 
            xcoords=as.double(xcoords[ind]), 
            ycoords=as.double(ycoords[ind]), tcoords=as.double(tcoords[ind]), 
            delta = as.double(delta[1]),count=as.integer(count),NAOK = TRUE,  DUP =FALSE)	
    
    loglik<--val$loglik[1:(val$count)]
    
    val$count
  }  
  
  
  return(list(ntot=sum(results.lik)))	
}


PLneg.spt.grad.recombine<-function(theta,ncores,ntot,eps = 1e-4,filename=NULL) {
  
  nparam<-length(theta)
  mat<-matrix(0,ntot,nparam)
  v<-numeric(ntot)
  count<-0
  for (j in 1:ncores) {
    load(paste(filename,"loglik-tmp-",j,".tmpRdata",sep=""))
    n<-length(loglik)
    sel<-seq(count+1,count+n)
    v[sel]<-loglik
    count<-count+n
  }
  for (k in 1:nparam) {
    count<-0
    for (j in 1:ncores) {
      load(file = paste(filename,"loglik-tmp-",k,"-",j,".tmpRdata",sep=""))
      n<-length(loglik)
      sel<-seq(count+1,count+n)
      mat[sel,k]<-loglik
      count<-count+n
     
    }  
    mat[,k]<-(mat[,k]-v)/eps
  }
  system(paste("rm ",filename,"loglik-tmp-*.tmpRdata",sep=""))
  return(mat)  
} 


j.hat.boot<-function(gradient.mat,ntot, size=1000,nboot =100) {
  
  nparam<-ncol(gradient.mat)
  
  start<-sample(0:(ntot-size),size = nboot,replace = FALSE)
  
  results<-matrix(0,nboot,nparam)
  for(i in 1:nboot) {
    sel<-(start[i]+1):(start[i]+size)
    g<-apply(gradient.mat[sel,],2,sum)
    results[i,]<-g
  }
  j.hat<-cov(results)/size
  
  return(j.hat)
}

j.hat.boot.parallel<-function(gradient.mat,ntot, size=1000,nboot =100,ncores=1) {
  
  nparam<-ncol(gradient.mat)
  
  start<-sample(0:(ntot-size),size = nboot,replace = FALSE)
  
  results<-matrix(0,nboot,nparam)
  results.mat<-foreach(i=1:nboot,.combine=rbind) %dopar%{ 
    sel<-(start[i]+1):(start[i]+size)
    g<-apply(gradient.mat[sel,],2,sum)
    g
  }
  #j.hat<-cov(results)/size
  j.hat<-cov(as.matrix(results.mat))/size
  return(j.hat)
}


j.hat.spt.gamma<-function(param, mask, y,z, threshold, xcoords, ycoords, tcoords,
                         delta,  ncores, nsites,ntimes, 
                         block.number=100, block.size =1000,
                         eps=1.0e-4) {
  
  if (block.size > ntimes)
    stop("block.size is too big")
  # y=(data[s_1,t_1],data[s_2,t_1,..., data[s_nsites,t_1], ..., data[s_nsites,t_ntimes])'
  # k number 
  # m number of times
  # k m-k overlapping observations  
  # g gradient
  lag<-delta[2]
  if (block.size < lag)
	  stop("block.size is too small")
  theta<-param[mask]
  ntheta<-length(theta)
  loglik<-rep(0,block.number)
  npairs<-rep(0,block.number)
  start<-sample(lag:(ntimes-block.size),replace = FALSE,size = block.number)
  mat<-foreach(i=1:block.number,.combine=rbind) %dopar%{ 
	sel <- (start[i]*nsites+1):((start[i]+block.size)*nsites)
	
	y1<-y[sel]
	z1<-z[sel]
	nobs<-length(y1)
  loglik[i]<-npairs[i]<-0
    val<-.C("pwl_spatiotemporal_parallel",y=as.double(y1),
            z=as.integer(z1), nobs=as.integer(nobs),
            nsites=as.integer(nsites), nblocks=as.integer(block.size), 
            lag=as.integer(lag), theta=as.double(param), loglik=as.double(loglik[i]),u=as.double(threshold), 
            xcoords=as.double(xcoords[sel]), 
            ycoords=as.double(ycoords[sel]), tcoords=as.double(tcoords[sel]), 
            delta = as.double(delta[1]),count=as.double(npairs[i]),NAOK = TRUE,  DUP =FALSE)	
	g<-rep(-val$loglik,ntheta)
	ntot<-val$count
  for (k in 1:ntheta) {
		dtheta <- theta
		dtheta[k] <- dtheta[k] + eps
		param[mask]<-dtheta
		loglik[i]<-npairs[i]<-0
		
		val<-.C("pwl_spatiotemporal_parallel",y=as.double(y1),z=as.integer(z1), 
		        nobs=as.integer(nobs),
            nsites=as.integer(nsites), nblocks=as.integer(block.size), 
            lag=as.integer(lag), theta=as.double(param), 
		        loglik=as.double(loglik[i]),u=as.double(threshold), 
            xcoords=as.double(xcoords[sel]), 
            ycoords=as.double(ycoords[sel]), 
		        tcoords=as.double(tcoords[sel]), 
            delta = as.double(delta[1]),
		        count=as.double(npairs[i]),NAOK = TRUE,  DUP =FALSE)	
			g[k]<-(-val$loglik-g[k])/eps # gradient      
		}
     c(g,ntot)
	}
  ntot<-mean(mat[,ntheta+1])
  j.hat<-cov(as.matrix(mat[,1:ntheta]))/ntot
  return(j.hat)
}
