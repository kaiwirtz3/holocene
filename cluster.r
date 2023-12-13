# create clusters
#
# kai wirtz (Hereon) 2023
#
# forms regions using an extended kmeans algorithm with optimized total number of clusters, based on SPD growth rates
#   saves clusti.mat
scdir='p3k14c/'
#acs=8 bcs=3 acs=2
acs=1
bcs=10
myargs = commandArgs(trailingOnly=TRUE)
print(paste('args:',length(myargs)))

if (length(myargs)>0) {
 XOUT=0
 ti0 = as.numeric(myargs[1])
 if (length(myargs)>1) {
   ti1 = as.numeric(myargs[2])
   } else { ti1=ti0 }
 if (ti1 < ti0) {ti1 = ti0}
 } else { #all sub-domains
 ti0=1
 ti1=18
 XOUT=1
 }
if (length(myargs)>2) {
  cscrit  = as.numeric(myargs[3])
  } else {
  cscrit= 120 # critical cluster size
  }
#for (cscrit in c(800 )) { #c(300,500,800,1000)
cc='Europe'

## XOUT=1
print(paste('time slice:',ti0,ti1,' XOUT=',XOUT,'cscrit=',cscrit))
# 1: run standalone (requires few library calls)
if (1) {
  library(sp)
  library(rworldmap)
  #library(rcarbon)
  library(cluster)
  library('R.matlab')
  base <- getMap(resolution="low") #extract basemap
  library("RColorBrewer")
  }
# read data
  file <- paste0(scdir,'bin/eurodat_all3','.bin')
  print(paste('read data from ',file))
  load(file)
# names of sub-continental patches
nclustsum=0

# statistical rigor
itern = 3000

#dt    = 9     # number of time segments
##print(max(pclust))
# loop over patches / data segments (entire set exceeds memory)
acs0=acs
lons=datc$lons
lats=datc$lats

 # maximal number of clusters
nn=length(lons)
print(paste('nn=',nn))
kmax=1+floor(nn/cscrit)
kmax=min(max(kmax,10),50)
print(paste('kmax=',kmax))
## TODO
clusti=rep(0,55000) # max number of dates

  # loop over time segments
for (ti in seq(ti0,ti1)) {
  tis=breaks[ti+1]
  print(tis)
  # withiness vector; cluster skill
  wi=rep(9,kmax+1)
  mcs=rep(0,kmax+1)
  mincsv=rep(0,kmax+1)

  # single cluster as reference
  cluster <- kmeans(datc[c(1:2,2+ti)], 1)
  wss0=cluster$tot.withinss

  # screen sequence of clusters
  #wi[1]=9E9

  for (k in 23+seq(kmax-23)) {
    wss=9E9
    # random iterations
  ##  for (it in seq(itern)){
    k0=which.min(wi)
    print(paste(wi[k-1],wi[k0]*3))
    if (wi[k-1]<wi[k0]*3){
    for (it0 in seq(itern)){
 # data are clustered by the k-means method:@
  # partition the points into k groups such that the sum of squares from points to the assigned cluster centres is minimized.
      cluster <- kmeans(datc[c(1:2,3+ti-ti0)], k)
      mclu=sort(cluster$size)
      mincs = mean(mclu[1:2]) # size of two smallest clusters

      # skill function depending on min cluster size
      #mf = 1 + exp(-(mincs*sqrt(sd(lon)*sd(lat))/(cscrit*acs)-1)*bcs)
      mf = 1 + exp(-(mincs/(cscrit*acs)-1)*bcs)
      ##mf = 1 + exp(-(mincs/(cscrit*8/sd(lon))-1)*3)

      # print(sd(lon)) print(mf)
      # lonstd(i)=std(lons(it));latstd(i)=std(lats(it));

      # compare with best partitioning among random samples
      ww=cluster$tot.withinss/wss0
      if (it0%%200==0) {
        X <- sprintf('%d %4d mc=%1.1f %1.5f %1.5f',k,it0,mincs,ww*mf,wss)
        print(X)}
      if (ww*mf < wss) {
        wi[k]=ww*mf
        mcs[k]=mf #
        mincsv[k]=mincs
        wss=wi[k]
        }
      } # end for it
     } # end if
    } # end for k
  # compute optimal cluster number
  k=which.min(wi)
  k=max(k,2)
  Xa <- paste(ti,'k=',k,'mc=',mincsv[k])
  print(Xa)
  wi[1]=wi[2]

  # plot clustering statistics
  if(XOUT==1) {
    x11(width = 9, height =8)
    } else {
    pfile <- paste0(scdir,'plots/cluster_withinss', ti,'_',cscrit,'.png')
    print(paste("plot to ",pfile))
    png(filename=pfile,width = 9, height =8, units = "cm", res=300 )
    }
  #par(oma = c(1, 1, 2, 1), mar = c(0.08, 0.08, 0.9, 0.3),mfrow=c(3,3), cex.lab=0.5,cex.sub=0.5, cex.main=1.,cex.axis=0.7)
  #print(wi)
  #tis=(breaks[ti]+breaks[ti+1])/2

  plot(wi,main=paste(cc[1],tis),log="y",xlim=c(10,kmax))
  abline(v = k, col="red", lwd=3)
  lines(mcs*0.5, col="blue", lwd=2)

  wss=9E9
  # repeat kmeans for optimal cluster number
  for (it in seq(itern)){
    cluster <- kmeans(datc[c(1:2,3+ti-ti0)], k)
    mclu = sort(cluster$size)
    ##print(mclu[1:4])
    mincs = mean(mclu[1:2])
    ##mf = 1 + exp(-(mincs/cscrit-1)*4)
    mf = 1 + exp(-(mincs/(cscrit*acs)-1)*bcs)

    if (cluster$tot.withinss/wss0*mf < wss) {
      wss=cluster$tot.withinss/wss0*mf
      cluc=cluster$cluster
      mincs0=mincs
      cluster0=cluster
      }
    } # end for it
    # spatial plot of clusters
   # xrange <- bbox(locs)[1,]
    #yrange <- bbox(locs)[2,]
    xrange <- c(min(lons)+3,max(lons)-3)
    yrange <- c(min(lats)+1,max(lats)-2)
    #print(k)
    cmap=palette(rainbow(k+4))
    cmap=palette(rainbow(k+4)) # needs to be invoked 2 times : WEIRD!
    #print(length(cmap))

    if(XOUT==1) {
      x11(width=11,height=13)
      } else {
      pfile <- paste0(scdir,'plots/clusters_', cc[1],tis,'.png')
      print(paste("plot to ",pfile))
      png(filename=pfile,width = 11, height = 13, units = "cm", res=600 )
      }
    par(oma = c(1, 1, 2, 1), mar = c(0.08, 0.08, 0.9, 0.3),mfrow=c(1,1), cex.lab=0.5,cex.sub=0.5, cex.main=1.,cex.axis=0.7)
   plot(base,col=rgb(0.95, 0.95, 0.95),border="antiquewhite3",xlim=xrange,ylim=yrange,main=paste(tis,"kBP"))

    if (0) {
      for (i in -2:2) {
          ii <- which(p == i)
         # length(ii)
          points(lons[ii],lats[ii], pch = 19, col = cmap[3+i], cex = .8)    }
      } # if 0
  sprintf('t=%d cluc= %d\n',tis,length(cluc))

  # plot all sites of a clustered region in a specific color
  for (i in seq(k)) {
    ii <- which(cluc == i)
    rgrm = mean(rmm[ii,1+ti-ti0], na.rm=TRUE)
    X <- sprintf('%d %d %1.3f',i,length(ii),rgrm)
    ##print(X)
    X <- sprintf('%d %1.1f',i,rgrm*10)

    points(lons[ii],lats[ii], pch = 19, col = cmap[i], cex = .3)
    text(mean(lons[ii])-1,mean(lats[ii]),X,cex=0.5+XOUT*1.)
    #print(length(pclust[ii]))
    clusti[pcl[ii]]=i
    nclustsum=nclustsum+length(ii)
    } # end for i
  text(-6,68,Xa,cex=1+XOUT*1.)

  # store results for each data segment in R binary file
  file <- paste0(scdir,'bin/clusti3_', cc[1],tis,'_',cscrit,'.bin')
  save(file=file,clusti,k,wi,cluc)
  writeMat(paste0(scdir,'mat/clusti3_',tis,'_',cscrit,'.mat'),clusti=clusti) #'_',cscrit,
  # export plot as PDF
  if(XOUT==1) {
    file <- paste0(scdir,'plots/clusters3_', cc[1],tis,'.pdf')
    print(paste("plot to",file))
    dev.copy2pdf(file=file, out.type = "pdf")
    } else {
    dev.off()
    }
  } # end for ti
