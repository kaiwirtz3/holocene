# collect grid cells and prepares clustering
#
# kai wirtz (Hereon) 2023
#
# forms regions using an extended kmeans algorithm with optimized total number of clusters, based on SPD growth rates
#   saves datc.bin
library(sp)

scdir='out/'  # main data IO directory

# range of grid patches
ri0=1
ri1=64

# statistical rigor
plim  = 0.05

# clean variables
lons=c()
lats=lons;pps=lons;rmm=lons;pcl=lons;sit=lons
s0=0

# loop over grid patches
for (ri in ri0:ri1)
  {
  file <- paste0(scdir,'bin/eurogrid_',ri,'.bin')
  if file.exists(file){
    load(file)
    print(paste(ri,length(pcl),length(pclust),'pclust=',pclust[1],pclust[length(pclust)]))

    # checks if patch contains sites/dates
    if(length(pclust)>0){

      pcl=c(pcl,pclust) # merges patch data

      locs=coordinates(eurospatial$locations)

      # position of sites
      lon = as.vector(locs[,1])
      lat = as.vector(locs[,2])

      lons= c(lons,lon) # merges sites' georefs
      lats= c(lats,lat)
      sit = c(sit,s0+sites)
      s0  = s0+sites[length(sites)]

      nn=length(eurospatial$pvalHi[,1])
      dt=length(eurospatial$pvalHi[1,])

      ##print(c(nn,length(lon),length(lons)))
      # create growth matrix
      pp=matrix(nrow = nn, ncol = dt)
      pp0=pp
      j=1

      # loop over time segments
      for (ti in seq(1,dt)) {

        pvalHi=eurospatial$pvalHi[,ti]
        pvalLo=eurospatial$pvalLo[,ti]
        qvalHi=eurospatial$qvalHi[,ti]
        qvalLo=eurospatial$qvalLo[,ti]

        # classify according p-value
        iHi=which(pvalHi<=plim)
        iLo=which(pvalLo<=plim)

        # classify also according q-value
        iHi2=which(qvalHi[iHi]<=plim)
        iLo2=which(qvalLo[iLo]<=plim)

        # mean rate of change
        #rt=round(rgr[ti]/2.5E-3)
        rt=0.5*asinh(rgr[ti]/2.E-3)
        if (is.na(rt) | is.nan(rt)) { rt=0 }
        #if (rt>2) { rt=2 }
        #if (rt<-2) { rt=-2 }

        fac=1
        p=rep(rt,nn)
        p0=p
        p[iHi]=p[iHi]+1*fac
        p[iHi[iHi2]]=p[iHi[iHi2]]+2*fac
        p[iLo]=p[iLo]-1*fac
        p[iLo[iLo2]]=p[iLo[iLo2]]-2*fac
        pp[,j]=p # discrete growth levels into matrix
        pp0[,j]=p0 # discrete growth levels into matrix
       ##   print(paste(ti,length(which(is.na(p0))),length(p0),length(which(is.na(p))),length(p)))

        j=j+1
        } # end for ti
      # generate data frame
      #print(paste(ri,length(which(is.na(pp0))),length(pp0),length(which(is.na(pp0))),length(pp0)))
      pps=rbind(pps,pp)
      rmm=rbind(rmm,pp0)
      } # if(length(pclust)>0)
    } # if file
  else {
    print(paste(file,'not found!'))
    }
  } # for ri

  # write binary data
  datc=data.frame(lons,lats,pps)
  file <- paste0(scdir,'bin/eurodat_all','.bin')
  print(paste('data ready to ',file))
  save(file=file,datc,rmm,pcl,breaks,sit)
