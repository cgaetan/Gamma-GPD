# simulation data
rm(list=ls())
filename<-"meteo"  
std.boot<-TRUE
sendmail<-FALSE
ncores<-32
coords<-read.table("coords.txt",header =FALSE)[,-1]
alldata<-read.table("alldata.txt",header =FALSE)
prob.delta<-0.6
prob<-0.99
kappa<-1/(1-prob)-1

# starting value
# parameters

semiax1<-151.7477000
semiax2<- 64.2817000
semiaxes<-c(semiax1,semiax2)
phi<- pi/6
duration<-7
velocity<-c(0.2,0.4)
beta<-1
alpha<-1
xi<-1
sigma<-kappa+1
param<-c(xi,sigma,semiaxes,  phi, duration, velocity, 
         alpha,beta, kappa)
mask<-rep(FALSE,length(param))
mask[c(3:8)]<-TRUE
delta.times<-10 


# starting value

block.size<-1000
block.number<-500

filename<-paste(filename,"-gamma-rev",sep = "")
source("skeleton-gamma-rev-2.R")
