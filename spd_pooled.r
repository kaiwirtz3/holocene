# calculates Summed Probability Density (SPD) and related growth for pooled dates
#   using an array of methods
#
# kai wirtz (Hereon) 2023
#
rm(list=ls())
library(rcarbon)  # RCARBON by Crema2017
library('R.matlab')
#library('doSNOW')
source("movavg.r")

# input/output directory
scdir='out/'
# list of SPD methods
#binv=c(50,200,100,100)

# pooled continent/region
cc='europe0'; news='_EU'
#cc='EA';news='_EA' cc='NIreland';news='_Ire'

# exponetial filter for smoothing
dt = 20
ff = exp(-2*(seq(-dt,dt)/dt)**2)
fil= ff/sum(ff)

breaks=seq(3000,9800,400)

# read C14 dates
tag <- paste0(news,'_NoNorm_Bin100') #(,'_Norm_Bin100','_NoNorm50','_NoNorm150'
dat <- readMat(paste0(scdir,'mat/C14_',cc,'_neo.mat'))
        # variables: 'lonsn','latsn','C14agesn','C14SDsn','SiteIDsn','datIDsn'

dt  = 10; dt2=2*dt
lat=round(dat$lat,digits=3)
lon=round(dat$lat,digits=3)

# prepare graphical output
file = paste0(scdir,"plots/spd_all3",tag,".png")
print(paste('new figure',file))
png(file,width=1200,height=840,units='px')
par(oma = c(1, 0, 1, 1), mar = c(1.2, 1, 0.14, 0.5),  cex.lab=1.5,cex.sub=1, cex.main=1.5,cex.axis=2)

ci=c()
sites0 <- data.frame(lat=lat,lon=lon)
usite  <- unique(sites0[,c('lat','lon')])
usitei <- as.numeric(rownames(usite))
usitei0<- usitei
usitei <- c(usitei,length(lat)+1)
sites  <- NULL
n0     <- length(usitei)

# loop over dates to fill site vector
for (ui in seq(n0-1)) {
    jj = usitei[ui]
    ni = which(lat==lat[jj] & lon==lat[jj])
    sites=c(sites,rep(ui,usitei[ui+1]-usitei[ui]))
    ci = c(ci,ni)
    }
print(paste('sites',length(sites),sites[1],sites[length(sites)],ui))

# filter NaNs
i2 <- which(!is.na(dat$C14agesn[ci]))
ii = ci[i2]
sites=sites[i2]
ndates=length(i2)

# print number of valid C14 dates per region
X <- sprintf('total: %d valid dates: %d ',length(lat),ndates)
print(X)

#  calibration using intcal20
#  run calibration over 4 cores
eurodates <- calibrate(dat$C14agesn[ii],dat$C14SDsn[ii],calCurves='intcal20',ncores=4, normalised=FALSE)
print(paste('eurodates ready ...'))

# create time bins
eurobins <- binPrep(sites=sites,ages=dat$C14agesn[ii],h=100)

# set time vector and calculate SPD
timeRange<- c(10000,3000) # time range calBP, older date first  c(10200,3000)
steps    <- seq(timeRange[1],timeRange[2],-dt) #
tirgr    <- steps[2:(length(steps)-1)]
clu_spd  <- spd(x = eurodates,bins=eurobins,timeRange=timeRange)

# calculate relative change rate (RGR) and smooth
clu_rgra <- spd2rc(clu_spd,breaks = steps)
rgr      <- clu_rgra$roca
rgrm     <- filter(rgr, filter = fil, method = 'convolution', circular = TRUE,sides=2)

# retrieve and smooth density
y  = clu_spd$grid$PrDens
ym = movavg(y, 20, type="n")

# store into matrix
mi = seq(round(0.5*dt2),floor(length(clu_spd$grid$calBP)/ dt2)* dt2, dt2)
ym = ym[mi]
tm = clu_spd$grid$calBP #+ timscale/2
tm2= tm[mi] #+ timscale/2

# save population data as Matlab binary
file = paste0(scdir,"mat/AllPop",tag,'_all.mat')
print(paste0("write data to ", file))
writeMat(file, poptime=tm2, ym=ym,trgr=tirgr,rgr=rgr,nsites=n0,ndates=ndates)

# plot SPD (with own method) and smoothed variants
plot(clu_spd,type="simple",col="indianred",lwd=1,lty=2,xlim=timeRange)
lines(tm,y,col="indianred",lwd=1,lty=1)
lines(tm2,ym,col="red",lwd=2,lty=1)

# plot RGR
fac=sd(ym)/sd(rgr)
lines(tirgr,mean(ym)+fac*rgr,col="blue",lwd=1,lty=2)
lines(tirgr,mean(ym)+fac*rgrm,col="blue",lwd=2,lty=1)
lines(tirgr,mean(ym)+0*rgr,col="blue",lwd=1,lty=2)

text(9000, mean(y)*0.05, labels = length(ii), cex = 2, col = NULL)
#text(ti+dtl*0.5, mean(y)*0.05, labels = c(as.character(cd[i,4]),':',as.character(cd[i,5])) , cex = 2, col = NULL)
dev.off()
