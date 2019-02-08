dyn.load(paste("gpd-latent",.Platform$dynlib.ext,sep=""))
require("RandomFields")
ellipse<-function(x0,y0,a,b,alpha=0,l=100)
{
  # the semi-axis length along the x-axis is a 
  # the semi-axis length along the y-axis is b for alpha=0
  # counterclockwise angle alpha:
  theta <- seq(0, 2 * pi, length=l)
  x <- x0 + a * cos(theta) * cos(alpha) - b * sin(theta) * sin(alpha)
  y <- y0 + a * cos(theta) * sin(alpha) + b * sin(theta) * cos(alpha)
  return(cbind(x,y))
}


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

# simulation of a space-time random field with GDP marginal distribution
gpdrf.gamma3D<-function(x1, y1, z1, model,  nx=50, ny=50, nz=50) {
  alpha <- model$alpha
  beta <- model$beta
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
  X<-gammarf3D(shape=alpha.gamma, rate=beta,xlim=xlim,ylim=ylim, zlim=zlim, nx=nx,ny=ny, nz=nz)
  # we calculate the rgamma surface
  # X contains the simulated values
  f<-check.C3D(x1,y1,z1,X,semiax1,semiax2,phi,duration,velocity)
  rft<-rexp(length(f),f)+threshold
  zt<-rbinom(length(f),1,prob=exp(-kappa*f))
  rf<-rep(NA,length(f))
  rf[zt==1]<-rft[zt == 1]
  return(list(z=rf,f=f,X=X,coords=cbind(x1,y1,z1)))
}


gpdrf.chi2<-function(model,coords,ntimes) {
    nsites<-nrow(coords)	
	kappa <- model$kappa
	psi<-model$psi
    threshold <- model$threshold
    RFoptions(seed=NA) ## *ANY* simulation will have the random seed 0; set
                       ##                   RFoptions(seed=NA) to make them all random again
 #model <- RMexp(scale=0.7,Aniso=RMangle(angle=pi/4, ratio=3))
 #model <- RMexp(scale=0.5)
 #simu<-RFsimulate(model, x=xy[,1], y=xy[,2],n = ntimes)
 # C(h,u)= (ψ(u)+1)^{-δ/2} φ(h /√(ψ(u) +1))
 # phi(h)=exp(-h^2) psi(u)=r^\alpha 
 model <- RMnsst(phi=RMgauss(scale = psi[1]), psi=RMfbm(alpha=1,scale=psi[2]),
                 delta=2)
 # Gneiting model 
 # /* \frac{1}{(dt^\alpha+1)^\beta}\exp\{\-frac{ds^\gamma}{(dt^\alpha+1)^{\beta\gamma/2}\}
 # alpha=1, beta=0.5, gamma=2
 model <- RMexp(scale = psi[1],proj = "space") * RMexp(scale = psi[2],proj = "time")
 simu1<-RFsimulate(model, x=coords[,1], y=coords[,2],T = 1:ntimes)
 f<-as.matrix(simu1)^2
 simu2<-RFsimulate(model, x=coords[,1], y=coords[,2],T = 1:ntimes)
 f<-f+as.matrix(simu2)^2
 f<-0.5*as.vector(t(f))
 ftmp<-as.vector(t(as.matrix(simu2)))
 xy<- coords %x% rep(1,ntimes)
 z1<-rep(1,nsites)%x% (1:ntimes)
 
 rft<-rexp(length(f),f)+threshold
 zt<-rbinom(length(f),1,prob=exp(-kappa*f))
 rf<-rep(NA,length(f))
 rf[zt==1]<-rft[zt == 1]
 
 
 return(list(z=rf,f=f,coords=cbind(x1=xy[,1],y1=xy[,2],z1=z1),ftmp=ftmp))
}	
	
gpdrf.chi2.2<-function(model,coords,ntimes) {
  nsites<-nrow(coords)	
  kappa <- model$kappa
  psi<-model$psi
  threshold <- model$threshold
  simu1 <- RFsim(coordx = coords[,1], coordy = coords[,2], 
                 coordt = 1:ntimes, corrmodel="exp_exp", grid=FALSE, 
               param=list(nugget=0,mean=0,scale_s=psi[1], 
                          scale_t=psi[2],sill=1))$data
  
  f<-simu1^2
  simu2 <- RFsim(coordx = coords[,1], coordy = coords[,2], 
                 coordt = 1:ntimes, corrmodel="exp_exp", grid=FALSE, 
                 param=list(nugget=0,mean=0,scale_s=psi[1], 
                            scale_t=psi[2],sill=1))$data
  f<-f+simu2^2
  f<-0.5*f
  
  xy<- coords %x% rep(1,ntimes)
  z1<-rep(1,nsites)%x% (1:ntimes)
  
  rft<-rexp(length(f),f)+threshold
  zt<-rbinom(length(f),1,prob=exp(-kappa*f))
  rf<-rep(NA,length(f))
  rf[zt==1]<-rft[zt == 1]
  
  
  return(list(z=rf,f=f,coords=cbind(x1=xy[,1],y1=xy[,2],z1=z1)))
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

rf.gaussian<-function(param,coords,ntimes, corrmodel, D=0.00000001*diag(2)) {
  nsites<-nrow(coords)	
  a<-param["alpha"]
  lambda<-param["lambda"]
  A <-diag(c(1, 1/sqrt(lambda))) %*% matrix(ncol=2, c(cos(a), sin(a), -sin(a), cos(a)))
  coordsnew<-t(A%*%t(coords))
  if (corrmodel == 2) {
    model <- RMnsst(phi=RMgauss(scale = param["phi.s"]), psi=RMfbm(alpha=1,scale=param["phi.t"]),
                    delta=2)  
  }
  if (corrmodel == 4) {
    model <- RMexp(scale = param["phi.s"], proj = "space") * RMexp(scale = param["phi.t"],proj = "time")
  }
  
  if (corrmodel == 7) {
	
	mu<-c(param["velocity1"],param["velocity2"])
    model <- RMcoxisham(RMexp(scale = param["phi1"]),mu,D=D)
  }
  if (corrmodel == 8) {
	mu<-c(param["velocity1"],param["velocity2"])
	
    model <- RMcoxisham(RMspheric(scale = param["phi1"]),mu,D=D)
  }
  simu<-RFsimulate(model, x=coordsnew[,1], y=coordsnew[,2],T = 1:ntimes)
  rf<-t(as.matrix(simu))
  
  
  return(rf)
}	


