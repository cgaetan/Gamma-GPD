# Simulation of a hierarchical space-time model
# For the simulation of the space-time gamma random field
# we adapt an algorithm in
# Wolpert, R. L., and Ickstadt, K. (1998b), “Simulation of 
# Lévy Random Fields,” in Practical Nonparametric
# and Semiparametric Bayesian Statistics, eds. D. Dey, 
# P. Müller, and D. Sinha, New York, NY: Springer New York,
# pp. 227–242.
dyn.load(paste("gpd-latent",.Platform$dynlib.ext,sep=""))
require("RandomFields")

e1<- function(x, d = 1e-9) {(2/d) * (1-pchisq(2*x, d))}
e1inv<- function(y, d = 1e-9) {qchisq(1-(d/2)*y, d)/2}



gammarf3D<- function(shape , rate, xlim=c(0,1),
                     ylim=c(0,1), zlim=c(0,1), nx ,ny, nz) {
  # jump locations 	
  x <- seq(xlim[1], xlim[2],length.out = nx)
  y <- seq(ylim[1], ylim[2],length.out = ny)
  z <- seq(zlim[1], zlim[2],length.out = nz)
  volume<-(x[2]-x[1])*(y[2]-y[1])*(z[2]-z[1])
  xyz<-expand.grid(x,y,z)
  n<-dim(xyz)[1]
  # a= shape b=rate 
  
  u <-  rgamma(n,shape=shape*volume,rate=rate)
  invisible(list(x=xyz[,1],y=xyz[,2],z=xyz[,3],u=u)) 
}


check.C3D<-function(x1,y1,z1,X,semiax1,semiax2,phi,duration,velocity)
{
  n<-length(x1) # number of observations
  m<-length(X$x) # number of elements of the grid
  f<-rep(0,n)
  val<-.C("integration3d",x = as.double(x1), y = as.double(y1), 
          z=as.double(z1), f = as.double(f),n = as.integer(n), 
          a = as.double(semiax1),b = as.double(semiax2), 
          phi =as.double(phi), duration =as.double(duration), 
          X = as.double(X$x), Y = as.double(X$y), Z = as.double(X$z), 
          U = as.double(X$u), V = as.double(velocity),  
          m = as.integer(m), DUP=FALSE)
  
  return(val$f)	
}

check.C3D.parallel<-function(x1,y1,z1,X,semiax1,semiax2,phi,duration,velocity,ncores)
{
  n<-length(x1) # number of observations
  nblocks<-n%/%ncores
  start<-seq(1,n,by=nblocks)
  end<-c(start[-1]-1,n)
  m<-length(X$x) # number of elements of the grid
  results<-foreach(j=1:length(start),.combine=c) %dopar%{ 
    n1<-end[j]-start[j]+1
    f<-rep(0,n1)
    val<-.C("integration3d",x = as.double(x1[start[j]:end[j]]), y = as.double(y1[start[j]:end[j]]), 
          z=as.double(z1[start[j]:end[j]]), f = as.double(f),n = as.integer(n1), 
          a = as.double(semiax1),b = as.double(semiax2), 
          phi =as.double(phi), duration =as.double(duration), 
          X = as.double(X$x), Y = as.double(X$y), Z = as.double(X$z), 
          U = as.double(X$u), V = as.double(velocity),  
          m = as.integer(m), DUP=FALSE)
  val$f
  }
  return(results)	
}





# this function simulates a Gamma random field
# we use an algorithm from Wolpert


gammarf3D.WI<- function(a, b, xlim=c(0,1),ylim=c(0,1), zlim=c(0,1), M=1000) {
  # jump locations 	
  x <- runif(M,min = xlim[1], xlim[2])
  y <- runif(M,min = ylim[1], ylim[2])
  z <- runif(M,min = zlim[1], zlim[2])
  # a = shape; b = rate  
  u <-  e1inv(cumsum(rexp(M)/a))/b
  u[(is.nan(u) | is.na(u) | is.infinite(u))]<-0
  invisible(list(x=x,y=y,z=z,u=u)) 
}

gpdrf.gamma3D.WI<-function(x1, y1,z1, model, M = 1000,parallel =FALSE, ncores =1) {
  # M number of points for the Poisson process
  alpha <- model$alpha
  beta <- model$beta
  # alpha= shape beta = rate 
  kappa <- model$kappa
  semiax1<-model$semiax1
  semiax2<-model$semiax2
  phi<-model$phi
  duration<- model$duration
  velocity<-model$velocity
  threshold <- model$threshold
  d<-max(c(semiax1,semiax2))
  xlim<-range(x1)+c(-(d+duration*abs(velocity[1])),(d+duration*abs(velocity[1])))
  ylim<-range(y1)+c(-(d+duration*abs(velocity[2])),(d+duration*abs(velocity[2])))
  zlim<-range(z1)+c(-duration,0)
  volume.ellipsoid<-pi*semiax1*semiax2*duration 
  volume<-diff(xlim)*diff(ylim)*diff(zlim)
  alpha.gamma<-alpha*volume/volume.ellipsoid
  
  # so the marginal distribution of \Gamma(dx)\sim  Gamma(alpha/volume.ellipsoid,beta)
  X<-gammarf3D.WI(a=alpha.gamma, beta,xlim=xlim,ylim=ylim, zlim=zlim, M=M)
  # we calculate the rgamma surface
  if (parallel) {
    f<-check.C3D.parallel(x1,y1,z1,X,semiax1,semiax2,phi,duration,velocity,ncores=ncores)   
  } else {
    f<-check.C3D(x1,y1,z1,X,semiax1,semiax2,phi,duration,velocity)   
  }
  
  
  rft<-rexp(length(f),f)+threshold
  zt<-rbinom(length(f),1,prob=exp(-kappa*f))
  rf<-rep(NA,length(f))
  rf[zt==1]<-rft[zt == 1]
  return(list(z=rf,f=f,coords=cbind(x1,y1,z1)	))
}

