source('functions-gamma.R')
library("parallel")
library("doParallel")
registerDoParallel(cores=ncores)
res<-list()
res$ncores<-ncores
res$param<-param
res$prob<-prob
d<-dist(coords)
delta.coords<-as.numeric(quantile(d,prob=prob.delta))
npairs<-sum(d <=delta.coords)
res$delta<-delta<-c(delta.coords,delta.times)
res$threshold<-threshold<-0
res$mask<-mask
res$prob<-prob
# we transform the data in each site !

tdata<-data2gpdcensored(data=alldata,prob=prob,kappa=param[11])
# we transpose because we use the parallel version of the program
tdata<-t(tdata)
ntimes<-ncol(tdata)
nsites<-nrow(tdata)
y<-as.numeric(tdata)
# warning !!! we consider only the exceedances so we set threshold = 0

save(res,file=paste(filename,"-results-starting-values.Rdata",sep="")) 
xy<-rep(1,ntimes)%x%as.matrix(coords)
tcoords<-(1:ntimes)%x%rep(1,nsites)




  # selection of the parameters
  print("all estimation step")
  mask<-res$mask

  threshold<-res$threshold
  delta<-c(delta.coords,delta.times)
  
  ptm <- proc.time()
  res<-sptgamma.fit(y=y, xcoords = xy[,1], ycoords =xy[,2], tcoords=tcoords,  
                        init = res$param, delta=delta,
                        threshold = threshold, mask=mask,ncores=ncores,nsites =nsites,ntimes=ntimes)
  res$prob<-prob
  res$elapsed<-proc.time() - ptm

  sink(paste(filename,"-results-estimation.txt",sep=""))
  print(res)
  sink()
  save(res,file=paste(filename,"-results-3rd-step.Rdata",sep=""))  
if (std.boot) {
  
 
  print("gamma std boot step")  
  
  load(file = paste(filename,"-results-3rd-step.Rdata",sep=""))  
  delta<-res$delta
  mask<-res$mask
  h.hat<-res$hessian 
  theta<-res$param[mask]
  res$prob<-prob

  parnames<-c("xi","sigma","semiax1","semiax2","phi","duration","velocity1","velocity2","alpha","beta","kappa")	
  z <- as.numeric((y > threshold))
  y[z == 0]<-res$threshold
  print("likelihood evaluation") 
  f<-PLneg.spt.parallel.wrap(theta=theta,y=y,z=z,threshold= threshold,xcoords = xy[,1], ycoords =xy[,2], tcoords=tcoords,  
                             delta=delta, param=param, mask=mask, ncores=ncores,nsites=nsites, ntimes=ntimes)
  res$ntot<-f$ntot
  res$negplik<-f$val # negative likelihood value    
  ptm <- proc.time()
  res$block.size<-block.size
  res$block.number<-block.number
  res$ncores<-ncores
  print("Start of the bootstrap")
  j.hat<-j.hat.spt.gamma(param = res$param, mask = res$mask,
                         y=y,z=z, threshold=res$threshold, 
                         xcoords =xy[,1], ycoords =xy[,2], 
                         tcoords=tcoords, 
                         delta =res$delta, 
                         ncores=ncores,
                        nsites = nsites,ntimes=ntimes, 
                        block.number = block.number,
                        block.size = block.size)  
  print("End of the boostrap")
  res$h.hat<-h.hat
  inv.h.hat<-ginv(h.hat)
  res$j.hat<-j.hat
  colnames(res$j.hat)<-parnames[res$mask]
  rownames(res$j.hat)<-parnames[res$mask]
  pen<-sum(diag((res$ntot*j.hat)%*%inv.h.hat))
  res$pen<-pen
  var.hat<-(inv.h.hat%*%(res$ntot*j.hat)%*%inv.h.hat)
  res$var.hat<-var.hat	
  colnames(res$var.hat)<-parnames[res$mask]
  rownames(res$var.hat)<-parnames[res$mask]
  res$std.err<-sqrt(diag(res$var.hat))
  res$clic<-res$negplik+pen # + because we have considered the negative value of cl		
  res$ncores<-ncores
  res$elapsed<-proc.time() - ptm
  save(res,file=paste(filename,"-results.Rdata",sep=""))
  print(res)
  sink(paste(filename,"-results-std-boot.txt",sep=""))
  print(res)
  sink()
 
}
