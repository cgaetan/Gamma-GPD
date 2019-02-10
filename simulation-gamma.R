rm(list=ls())

source('functions-gamma-simulation-WI.R')
source('functions-gamma.R')

# we use 6 cores
library("parallel")
library("doParallel")
ncores<-6 # number of cores 
registerDoParallel(cores=ncores)

nsites<-30 # number of sites
ntimes<-2000 # number of repetitions 
prob<-0.90 # Pr(Z <= th)=p where th is the threshold 
          #  at which the model is estimated
prob.sim<-0.8  # Pr(Z <= th)=p where th is the threshold 
          #  at which the data are simulated 
# we choose  locations in a square [0,1]x[0,1] at random
xlim<-c(0,1)
ylim<-c(0,1)
# times start=1, stop=ntimes
tlim<-c(1,ntimes) 
nx<-20
ny<-20
x<-seq(xlim[1],xlim[2],length.out = nx)
y<-seq(ylim[1],ylim[2],length.out = ny)
coords<-expand.grid(x,y) # we build a grid
set.seed(1)
sel<-sample(1:(nx*ny),nsites) 
coords<-as.matrix(coords[sel,]) # we randomly select nsites points of the grid
d<-dist(coords)


####################################
# initialization of the parameters #
####################################
phi<-pi/4
semiax1<-0.26
semiax2<-0.2
duration<-10
velocity<-c(0.1,0.1)
area.ellipsoid<-pi*semiax1*semiax2*duration # so the marginal distribution is Gamma(alpha*pi*radius[1]*radius[2]*duration,beta)
beta<-1
alpha<-1
kappa<-(1-prob)^(-1/(alpha))*beta-beta
kappa.sim<-(1-prob.sim)^(-1/(alpha))*beta-beta
threshold<-0 
xi<-1 # GPD shape parameter
sigma<-kappa+1 # GPD scale parameter
semiaxes<-c(semiax1,semiax2)
param<-c(xi,sigma,semiaxes,phi,duration,velocity,alpha,beta,kappa)
# list of parameters in the simulation stage
model.sim<-list(alpha=alpha,beta=beta,kappa=kappa.sim,semiax1=semiax1,semiax2=semiax2,phi=phi,
                duration=duration, threshold=threshold,velocity=velocity, threshold=threshold)

# list of parameters in the estimation stage
model<-list(alpha=alpha,beta=beta,kappa=kappa,semiax1=semiax1,semiax2=semiax2,phi=phi,
            duration=duration, threshold=threshold,velocity=velocity, threshold=threshold)

# we simulate from the space-time model
z<-seq(tlim[1],tlim[2],length.out = ntimes)
xyz<-cbind(coords%x%rep(1,ntimes),rep(z,nsites))		
x1<-xyz[,1]
y1<-xyz[,2]
z1<-xyz[,3]  

dlim<-max(c(semiax1,semiax2))
xlim<-range(x1)+c(-(dlim+duration*abs(velocity[1])),(dlim+duration*abs(velocity[1])))
ylim<-range(y1)+c(-(dlim+duration*abs(velocity[2])),(dlim+duration*abs(velocity[2])))
zlim<-range(z1)+c(-duration,0)
volume.ellipsoid<-pi*semiax1*semiax2*duration 
volume<-diff(xlim)*diff(ylim)*diff(zlim)
# number of points to be simulated in the whole volume
M<-round(20*volume/volume.ellipsoid)
maxdist<-max(dist(coords))
xy<-rep(1,ntimes)%x%coords
tcoords<-(1:ntimes)%x%rep(1,nsites)


prob.delta<-1 # a simple way for selecting \Delta_S in formula (11), i.e. all distances
delta.coords<-as.numeric(quantile(d,prob=prob.delta))
delta.times<-10
delta<-c(delta.coords,delta.times)
################################
# simulation of the exeedances #
################################
X3D<-gpdrf.gamma3D.WI(x1=x1,y1=y1,z1=z1, model=model.sim, M = M,parallel = TRUE,ncores = ncores) 
X3D$z[is.na(X3D$z)]<-0
alldata<-matrix(X3D$z,ntimes,nsites)	
# parameters to be estimated using maximum CL
mask<-rep(FALSE,length(param))
mask[c(3:8)]<-TRUE
# we stock all results in the list res
res<-list()
res$ncores<-ncores
res$param<-param

res$param[mask]<-c(0.2,0.1,pi/5,8,0.04,0.05) # we set the initial values.
res$prob<-prob
res$delta<-delta
res$mask<-mask
res$prob<-prob
# we transform the data in each site !
tdata<-data2gpdcensored(data=alldata,prob=prob,kappa=param[11])
# we transpose because we use the parallel version of the program
tdata<-t(tdata)
y<-as.numeric(tdata)
mask<-res$mask
delta<-c(delta.coords,delta.times)
maxit<-2000 # number of iterations for the optimization algorithm
method<-"Nelder-Mead" # method for optim.   See ?optim for details

#################################
# estimation of the parameters  #
#################################

res<-sptgamma.fit(y=y, xcoords = xy[,1], ycoords =xy[,2], tcoords=tcoords,  
                  init = res$param, delta=delta,
                  threshold = threshold, maxit=maxit,method = method,
                  mask=mask,ncores=ncores,
                  nsites =nsites,ntimes=ntimes)

##############
# comparison #
##############

# true parameters
print(param[mask])
# starting values of the parameters
print(res$start)
# estimated parameters
print(res$param[mask])
