# create region clusters of C14 sites
#
# kai wirtz (Hereon) 2023
#
# forms regions using an extended kmeans algorithm with optimized total number of clusters, based on SPD growth rates
#   saves clusti.mat
scdir='out/'

# checks for input arguments that control partial computation
#    in a parallel mode
# Examples: "Rscript cluster.r $SLURM_ARRAY_TASK_ID"
#                within a slurm.sh
#           "Rscript cluster.r 1 3"
#                to process only time segments 1 to 3

# retrieve time segments from arguments or sets default range
myargs  <-  commandArgs(trailingOnly=TRUE)
#print(paste('args:',length(myargs)))
if (length(myargs)>0) {
 ti0  <-  as.numeric(myargs[1])
 if (length(myargs)>1) {
   ti1  <-  as.numeric(myargs[2])
   } else { ti1 <- ti0 }
 if (ti1 < ti0) {ti1  <-  ti0}
 } else { #all sub-domains
 ti0 <- 1
 ti1 <- 18
 }

# set clustering parameter cscrit
if (length(myargs)>2) {
  cscrit   <-  as.numeric(myargs[3])
  } else {
  cscrit <-  120 # critical cluster size
  }

print(paste('time slice:',ti0,ti1,'cscrit=',cscrit))

cc   <- 'Europe'
XOUT <- 0  # no output to window

# 1: run standalone (requires few library calls)
if (1) {
  library(sp)
  library(rworldmap)
  library(stats)
#  library(cluster)
  library('R.matlab')
  base <- getMap(resolution="low") #extract basemap
  library("RColorBrewer")
  }

# read data
file <- paste0(scdir,'bin/eurodat_all','.bin')
print(paste('read data from ',file))
load(file)

#
nclustsum <- 0
bcs    <- 10 # cluster parameter
itern  <- 1000 # statistical rigor, should be >1000 and <10000, depends on computing resource

lons <- datc$lons
lats <- datc$lats

 # maximal number of clusters
nn   <- length(lons)
kmax <- 1+floor(nn/cscrit)
kmax <- min(max(kmax,10),50)
## TODO
clusti <- rep(0,92000) # max number of dates

# loop over time segments
for (ti in seq(ti0,ti1)) {
  tis <-breaks[ti+1]
  print(tis)
  # withiness vector; cluster skill
  wi    <- rep(9,kmax+1)
  mcs   <- rep(0,kmax+1)
  mincsv<- rep(0,kmax+1)

  # single cluster as reference
  cluster <- kmeans(datc[c(1:2,2+ti)], 1)
  wss0 <- cluster$tot.withinss

  # screen sequence of clusters
  for (k in 26+seq(kmax-23)) {
    wss <- 9E9
    # random iterations
    k0 <- which.min(wi)
    #print(paste(wi[k-1],wi[k0]*3))
    if (wi[k-1]<wi[k0]*3){

    # loop over random iterations
      for (it0 in seq(itern)){

    # data are clustered by the k-means method
    # partition the points into k groups such that
    #  the sum of squares from points to the assigned cluster centres is minimized.
        cluster <- kmeans(datc[c(1:2,3+ti-ti0)], k)
        mclu    <- sort(cluster$size)
        mincs   <- mean(mclu[1:2]) # size of two smallest clusters

        # skill function depending on minimal cluster size
        mf  <- 1 + exp(-(mincs/(cscrit)-1)*bcs)

        # compare with best partitioning among random samples
        ww  <- cluster$tot.withinss/wss0

        # interim output
        if (it0%%150==0) {
          X <- sprintf('%d %4d mc=%1.1f %1.3f %1.3f',k,it0,mincs,ww*mf,wss)
          print(X)
          }
        # store best cluster layout
        if (ww*mf < wss) {
          wi[k] <- ww*mf
          mcs[k]<- mf #
          mincsv[k] <- mincs
          wss   <- wi[k]
          }
        } # end for iteration
      } # end if
    } # end for k
  # compute optimal cluster number
  k <- which.min(wi)
  k <- max(k,2)
  Xa<- paste(ti,'k=',k,'mc=',mincsv[k])
  print(Xa)
  wi[1] <- wi[2]

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
  plot(wi,main=paste(cc[1],tis),log="y",xlim=c(10,kmax))
  abline(v = k, col="red", lwd=3)
  lines(mcs*0.5, col="blue", lwd=2)

  wss <-9E9

  # repeat kmeans for optimal cluster number

  for (it in seq(itern)){
    cluster <- kmeans(datc[c(1:2,3+ti-ti0)], k)
    mclu  <- sort(cluster$size)
    mincs  <- mean(mclu[1:2])
    mf  <- 1 + exp(-(mincs/(cscrit)-1)*bcs)

    if (cluster$tot.withinss/wss0*mf < wss) {
      wss <-cluster$tot.withinss/wss0*mf
      cluc <-cluster$cluster
      mincs0 <-mincs
      cluster0 <-cluster
      }
    } # end for it
    # spatial plot of clusters
    xrange <- c(min(lons)+3,max(lons)-3)
    yrange <- c(min(lats)+1,max(lats)-2)

    #print(k)
    cmap <-palette(rainbow(k+4))
    cmap <-palette(rainbow(k+4)) # needs to be invoked 2 times : WEIRD!

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
          points(lons[ii],lats[ii], pch = 19, col = cmap[3+i], cex = .8)    }
      } # if 0
   sprintf('t=%d cluc= %d\n',tis,length(cluc))

  # plot all sites of a clustered region in a specific color
  for (i in seq(k)) {
    ii <- which(cluc == i)
    rgrm  <- mean(rmm[ii,1+ti-ti0], na.rm=TRUE)
    X <- sprintf('%d %d %1.3f',i,length(ii),rgrm)
    ##print(X)
    X <- sprintf('%d %1.1f',i,rgrm*10)

    points(lons[ii],lats[ii], pch = 19, col = cmap[i], cex = .3)
    text(mean(lons[ii])-1,mean(lats[ii]),X,cex=0.5+XOUT*1.)
    #print(length(pclust[ii]))
    clusti[pcl[ii]] <-i
    nclustsum <-nclustsum+length(ii)
    } # end for i
  text(-6,68,Xa,cex=1+XOUT*1.)

  # store results for each data segment in R binary file
  file <- paste0(scdir,'bin/clusti_', cc[1],tis,'_',cscrit,'.bin')
  save(file=file,clusti,k,wi,cluc)

  # store cluster index of dates for each data segment in Matlab binary file
  writeMat(paste0(scdir,'mat/clusti_',tis,'_',cscrit,'.mat'),clusti=clusti) #'_',cscrit,

  # export plot as PDF
  if(XOUT==1) {
    file <- paste0(scdir,'plots/clusters_', cc[1],tis,'.pdf')
    print(paste("plot to",file))
    dev.copy2pdf(file=file, out.type = "pdf")
    } else {
    dev.off()
    }
  } # end for ti
