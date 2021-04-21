##Test using Bivariate Conditional Bias Correction with one cell of data

library(scales)
library(ncdf4)
library(PCICt)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----------------------------------------------------
##Extract a time series of BCCI 


extract_series <- function(nc,var.name,lonc,latc) {

   lon <- ncvar_get(nc,'lon')   
   lat <- ncvar_get(nc,'lat')   

   lon.ix <- which.min(abs(lonc-lon))
   lat.ix <- which.min(abs(latc-lat))

   cell <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))

   time.series <- netcdf.calendar(nc)

   rv <- list(time=time.series,data=cell)
   return(rv)
}

##----------------------------------------------------
##****************************************************

##Kelowna
##lonc <- -119.4
##latc <- 49.9

##Prince George
site <- 'Prince-George'
lonc <- -122.75
latc <- 54.92

var.name <- 'tasmin'
plot.all <- FALSE

save.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/save_dir/" 
bcci.dir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/mvbc/BCCI/"
bccaq.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"
gcm.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/gcm_grid/"
qdm.dir <-  "/storage/data/climate/downscale/BCCAQ2+PRISM/mvbc/QDM/"

##Plot time series of max temp with annual means
plot.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/plots/"


##------------------------------------------------------
##CanRCM4 Acting as observations

rcm.file <- paste0(bccaq.dir,var.name,"_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(rcm.file)

tas.obs.series <- extract_series(nc,var.name,lonc,latc)
##save.file <- paste0(save.dir,var.name,"_CanRCM4_BC_RCM-Grid_Kelowna_cell_1950-2100.RData")
##save(series,file=save.file)

ost <- grep('1951-01-01',tas.obs.series$time)
oen <- grep('2005-12-31',tas.obs.series$time)

tas.obs <- tas.obs.series$data[ost:oen]
tas.bins <- quantile(tas.obs,c(0.2,0.4,0.6,0.8))

nc_close(nc)
tas.cal.index <- ost:oen
tas.cal.index[ tas.obs <= tas.bins[1]] <- 1
tas.cal.index[ tas.obs > tas.bins[1] & tas.obs <= tas.bins[2]] <- 2
tas.cal.index[ tas.obs > tas.bins[2] & tas.obs <= tas.bins[3]] <- 3
tas.cal.index[ tas.obs > tas.bins[3] & tas.obs <= tas.bins[4]] <- 4
tas.cal.index[ tas.obs > tas.bins[4]] <- 5

tas.cal.factor <- as.factor(tas.cal.index)

var.name <- 'pr'
rcm.file <- paste0(bccaq.dir,var.name,"_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(rcm.file)

pr.canrcm4.series <- extract_series(nc,var.name,lonc,latc)
pr.canrcm4.series$data <- pr.canrcm4.series$data * 86400
pr.obs <- pr.canrcm4.series$data[ost:oen]


##save.file <- paste0(save.dir,var.name,"_CanRCM4_BC_RCM-Grid_Kelowna_cell_1950-2100.RData")
##save(series,file=save.file)

nc_close(nc)

##------------------------------------------------------
var.name <- 'pr'
bccaq.file <- paste0(bccaq.dir,var.name,"_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")

nc <- nc_open(bccaq.file)

pr.bccaq.series <- extract_series(nc,var.name,lonc,latc)
###save.file <- paste0(save.dir,var.name,"_BCCAQv2_BC_RCM-Grid_Kelowna_cell_1950-2100.RData")
##save(tasmax.series,file=save.file)
nc_close(nc)

var.name <- 'tasmin'
##bccaq.file <- paste0(bcci.dir,var.name,"_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
tas.ci.file <- paste0(bcci.dir,var.name,"_BCCI_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(tas.ci.file)

tas.ci.series <- extract_series(nc,var.name,lonc,latc)
save.file <- paste0(save.dir,var.name,"_CI_BC_RCM-Grid_",site,"_cell_1950-2100.RData")
save(tas.ci.series,file=save.file)
nc_close(nc)

var.name <- 'pr'
bcci.file <- paste0(bcci.dir,var.name,"_BCCI_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(bcci.file)
pr.ci.series <- extract_series(nc,var.name,lonc,latc)
save.file <- paste0(save.dir,var.name,"_BCCI_BC_RCM-Grid_",site,"_cell_1950-2100.RData")
save(pr.ci.series,file=save.file)

gcm.file <- paste0(gcm.dir,var.name,"_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(gcm.file)
pr.gcm.series <- extract_series(nc,var.name,lonc,latc)
pr.gcm.series$data <- pr.gcm.series$data*86400
save.file <- paste0(save.dir,var.name,"_BC_GCM-Grid_",site,"_cell_1950-2100.RData")
save(pr.gcm.series,file=save.file)

nc_close(nc)

##GCM Series
tas.gcm.file <- paste0(gcm.dir,"tasmax_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(tas.gcm.file)
tx.gcm.series <- extract_series(nc,'tasmax',lonc,latc)
nc_close(nc)

tn.gcm.file <- paste0(gcm.dir,"tasmin_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(tn.gcm.file)
tn.gcm.series <- extract_series(nc,'tasmin',lonc,latc)
nc_close(nc)

##QDM Series

qdm.file <- paste0(qdm.dir,"tasmax_QDM_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(qdm.file)
tx.qdm.series <- extract_series(nc,'tasmax',lonc,latc)
nc_close(nc)

qdm.file <- paste0(qdm.dir,"tasmin_QDM_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
nc <- nc_open(qdm.file)
tn.qdm.series <- extract_series(nc,'tasmin',lonc,latc)
nc_close(nc)


##-----------------------------------------------------------------------
##***********************************************************************
tas.series <- tx.qdm.series
pr.series <- pr.gcm.series


##***********************************************************************
##-----------------------------------------------------------------------
##Plot annual maxima, annual total, annual mean precip from BCCI, RCM, GCM

year.fac <- as.factor(format(pr.gcm.series$time,'%Y'))
 
##Annual mean comparison
plot.file <- paste0(plot.dir,'pr_annual_mean_RCM-GCM-CI_',site,'_cell_1950-2100.png')

png(file=plot.file,height=4,width=6,res=150,units='in')

pcol <- 'gray22'
years <- levels(year.fac)

pr.gcm.years <- tapply(pr.gcm.series$data,year.fac,mean)
pr.ci.years <-  tapply(pr.ci.series$data,year.fac,mean)
pr.rcm.years <-  tapply(pr.canrcm4.series$data,year.fac,mean)

plot(years,pr.gcm.years,type='l',col='blue',lwd=2,ylim=c(1.5,5),
     col.axis=pcol,col.lab=pcol,xlab='Date',ylab='Ann. Pr. (mm)',
     main='Annual Mean Precip Time Series',col.main=pcol)
abline(h=0:6,col='lightgray',lwd=0.5,lty=2)
lines(years,pr.gcm.years,col='blue',lwd=2)
lines(years,pr.ci.years,col='red',lwd=1.5)
lines(years,pr.rcm.years,col='darkgreen',lwd=1)
legend('topleft',leg=c('GCM','CI','RCM'),col=c('blue','red','darkgreen'),pch=15)

box(which='plot',col=pcol)

dev.off()

##Annual max comparison
plot.file <- paste0(plot.dir,'pr_annual_max_RCM-GCM-CI_',site,'_cell_1950-2100.png')

png(file=plot.file,height=4,width=6,res=150,units='in')

pcol <- 'gray22'
years <- levels(year.fac)

pr.gcm.years <- tapply(pr.gcm.series$data,year.fac,max)
pr.ci.years <-  tapply(pr.ci.series$data,year.fac,max)
pr.rcm.years <-  tapply(pr.canrcm4.series$data,year.fac,max)

plot(years,pr.gcm.years,type='l',col='blue',lwd=2,
     col.axis=pcol,col.lab=pcol,xlab='Date',ylab='Ann. Pr. (mm)',
     main='Annual Max Precip Time Series',col.main=pcol)
abline(h=seq(0,100,20),col='lightgray',lwd=0.5,lty=2)
lines(years,pr.gcm.years,col='blue',lwd=2)
lines(years,pr.ci.years,col='red',lwd=1.5)
lines(years,pr.rcm.years,col='darkgreen',lwd=1)
legend('topleft',leg=c('GCM','CI','RCM'),col=c('blue','red','darkgreen'),pch=15)

box(which='plot',col=pcol)

dev.off()

##Annual Total comparison
plot.file <- paste0(plot.dir,'pr_annual_total_RCM-GCM-CI_',site,'_cell_1950-2100.png')

png(file=plot.file,height=4,width=6,res=150,units='in')

pcol <- 'gray22'
years <- levels(year.fac)

pr.gcm.years <- tapply(pr.gcm.series$data,year.fac,sum)
pr.ci.years <-  tapply(pr.ci.series$data,year.fac,sum)
pr.rcm.years <-  tapply(pr.canrcm4.series$data,year.fac,sum)

plot(years,pr.gcm.years,type='l',col='blue',lwd=2,ylim=c(600,1700),
     col.axis=pcol,col.lab=pcol,xlab='Date',ylab='Ann. Pr. (mm)',
     main='Annual Total Precip Time Series',col.main=pcol)
abline(h=seq(400,1800,200),col='lightgray',lwd=0.5,lty=2)
lines(years,pr.gcm.years,col='blue',lwd=2)
lines(years,pr.ci.years,col='red',lwd=1.5)
lines(years,pr.rcm.years,col='darkgreen',lwd=1)
legend('topleft',leg=c('GCM','CI','RCM'),col=c('blue','red','darkgreen'),pch=15)

box(which='plot',col=pcol)

dev.off()




##-----------------------------------------------------------------------




year.fac <- as.factor(format(tas.series$time,'%Y'))
tas.ann <- tapply(tas.series$data,year.fac,mean)
y.ix <- grep('*-07-15',tas.series$time)

dates <- as.Date(as.character(tas.series$time))
plot.file <- paste0(plot.dir,var.name,'_',site,'_cell_1950-2100.png')

if (plot.all) {

png(file=plot.file,height=4,width=6,res=150,units='in')

pcol <- 'gray22'

plot(dates,tas.series$data,pch=16,col=alpha('darkorange',0.35),
     col.axis=pcol,col.lab=pcol,xlab='Date',ylab='Max Temp (degC)',
     main='TAS Time Series',col.main=pcol)
lines(dates[y.ix],tas.ann,col='red',lwd=2)
abline(h=seq(-50,50,10),col='lightgray',lty=2,lwd=0.5)
box(which='plot',col=pcol)

dev.off()
}

##Detrend Tasmax future 2006-2100

yst <- grep('2005-01-01',tas.series$time)
yen <- length(tas.series$time)

tas.df <- data.frame(x=dates[yst:yen],y=tas.series$data[yst:yen])
tas.fit <- lm(y~x,tas.df)

past.mean <- mean(tas.series$data[1:(yst-1)])
future.trend <- tas.fit$fitted.values
base <- c(rep(past.mean,length(1:(yst-1))),future.trend)

tas.detrended <- tas.series$data - base
tas.dt.ann <- tapply(tas.detrended,year.fac,mean)

plot.file <- paste0(plot.dir,'tasmin_detrended_',site,'_cell_1950-2100.png')
if (plot.all) {
png(file=plot.file,height=4,width=6,res=150,units='in')

plot(dates,tas.detrended,pch=16,col=alpha('goldenrod',0.35),
     col.axis=pcol,col.lab=pcol,xlab='Date',ylab='Max Temp (degC)',
     main='TAS Detrended Time Series',col.main=pcol)
lines(dates[y.ix],tas.dt.ann,col='black',lwd=2)
lines(dates[y.ix],tas.ann-rep(past.mean,151),col='red',lwd=2)
abline(h=seq(-50,50,10),col='lightgray',lty=2,lwd=0.5)
box(which='plot',col=pcol)

dev.off()
}
##Divide Tasmax into 5 bins 

tas.bins <- quantile(tas.detrended,c(0.2,0.4,0.6,0.8))

##Find index for each bin
tas.qt.index <- 1:yen
tas.qt.index[ tas.detrended <= tas.bins[1]] <- 1
tas.qt.index[ tas.detrended > tas.bins[1] & tas.detrended <= tas.bins[2]] <- 2
tas.qt.index[ tas.detrended > tas.bins[2] & tas.detrended <= tas.bins[3]] <- 3
tas.qt.index[ tas.detrended > tas.bins[3] & tas.detrended <= tas.bins[4]] <- 4
tas.qt.index[ tas.detrended > tas.bins[4]] <- 5

##Histograms for Tas bins

plot.file <- paste0(plot.dir,'tasmin_detrended_',site,'_cell_5bin_histograms.png')
if (plot.all) {
png(file=plot.file,height=6,width=6,res=150,units='in')

par(mfrow=c(3,2))

for (i in 1:5) {
   tas.bin <- tas.detrended[tas.qt.index==i]
   breaks <- seq(floor(min(tas.bin)/5)*5,ceiling(max(tas.bin)/5)*5,1) ##pretty(tas.bin) ##
   print(breaks)
   hb <- hist(tas.bin,breaks=breaks,plot=F)

   plot(hb,freq=F,border='red',col=alpha('red',0.35),
     col.axis=pcol,col.lab=pcol,xlab='Max Temp (degC)',ylab='Density',
     main=paste0('TASMIN Bin ',i),col.main=pcol)

   box(which='plot',col=pcol)
} 
dev.off()
}

##Histograms of corresponding Pr bins

plot.file <- paste0(plot.dir,'pr_NONZERO_',site,'_cell_5bin_histograms.png')
##plot.file <- paste0(plot.dir,'pr_',site,'_cell_5bin_histograms.png')
if (plot.all) {
png(file=plot.file,height=6,width=6,res=150,units='in')

par(mfrow=c(3,2))

for (i in 1:5) {
   pr.bin <- pr.series$data[tas.qt.index==i]
   pr.nz <- pr.bin[pr.bin > 1]
   pr.bin <- pr.nz
   breaks <- seq(floor(min(pr.bin)/5)*5,ceiling(max(pr.bin)/5)*5,1) ##pretty(pr.bin,n=10) 
   print(breaks)
   hb <- hist(pr.bin,breaks=breaks,plot=F)

   plot(hb,freq=F,border='blue',col=alpha('blue',0.35),
     col.axis=pcol,col.lab=pcol,xlab='Precipitation (mm)',ylab='Density',
     main=paste0('PRECIP Bin ',i),col.main=pcol)

   box(which='plot',col=pcol)
} 
dev.off()
}

##------------------------------------------------------------
##Set up QDM to operate on each bin - no time division
library(ClimDown)

cstart <- getOption('calibration.start')
cend <- getOption('calibration.end')

print("Opening the input files and reading metadata")
gcm <- nc_open(gcm.file)
obs <- nc_open(rcm.file)

lat <- gcm$dim$lat$vals
lon <- gcm$dim$lon$vals

gcm.time <- compute.time.stats(gcm)
obs.time <- compute.time.stats(obs, cstart, cend)

gcm.obs.subset.i <- gcm.time$vals > as.PCICt(cstart, attr(gcm.time$vals, 'cal')) & gcm.time$vals < as.PCICt(cend, attr(gcm.time$vals, 'cal'))

varname <- 'pr'
time.factors <- mk.factor.set(obs.time$vals,
                                  gcm.time$vals[gcm.obs.subset.i],
                                  gcm.time$vals,
                                  multiyear=getOption('multiyear'),
                                  seasonal=getOption('seasonal')[[varname]],
                                  n.multiyear=getOption('multiyear.window.length'),
                                  expand.multiyear=getOption('expand.multiyear')
                                  )

mp.years <- substr(as.character(time.factors$mp),0,4)

bcbc.mp <- as.factor(paste0(mp.years,'-',tas.qt.index))

bcbc.factors <- list(oc=tas.cal.factor,
                     mc=as.factor(tas.qt.index[gcm.obs.subset.i]),
                     mp=bcbc.mp)

##Test QDM

varname <- 'pr'
pr.qdm.test <- tQDM(o.c=pr.obs, m.c=pr.series$data[gcm.obs.subset.i], m.p=pr.series$data,
                  bcbc.factors$oc,
                  bcbc.factors$mc,
                  bcbc.factors$mp,
                  ratio=getOption('ratio')[[varname]],
                  trace=getOption('trace'),
                  jitter.factor=getOption('jitter.factor'),
                  n.tau=getOption('tau')[[varname]])

##Annual comparison
plot.file <- paste0(plot.dir,'pr_tx_annual_mean_BCBC-from-GCM_QDM_',site,'_cell_1950-2100.png')

png(file=plot.file,height=4,width=6,res=150,units='in')

pcol <- 'gray22'
years <- levels(year.fac)

pr.bccaq.years <- tapply(pr.bccaq.series$data,year.fac,mean)
pr.bcbc.years <-  tapply(pr.qdm.test,year.fac,mean)
pr.rcm.years <-  tapply(pr.canrcm4.series$data,year.fac,mean)

plot(years,pr.bccaq.years,type='l',col='blue',lwd=2,
     col.axis=pcol,col.lab=pcol,xlab='Date',ylab='Ann. Pr. (mm)',
     main='Annual Mean Precip (Pr|Tasmax) from GCM Time Series',col.main=pcol)
lines(years,pr.bcbc.years,col='red',lwd=1)
lines(years,pr.rcm.years,col='darkgreen',lwd=1)
legend('topleft',leg=c('BCCAQv2','BCBC','RCM'),col=c('blue','red','darkgreen'),pch=15)

box(which='plot',col=pcol)

dev.off()

##Plot one year of precip
year <- 2099
##plot.file <- paste0(plot.dir,'pr_one_year_',year,'_BCCAQv2_BCBC_Kelowna_cell_1950-2100.png')
##plot.file <- paste0(plot.dir,'pr_one_year_',year,'_BCCAQv2_CanRCM4_Kelowna_cell_1950-2100.png')

plot.file <- paste0(plot.dir,'pr_one_year_',year,'_BCBC-from-GCM_CanRCM4_',site,'_cell_1950-2100.png')

png(file=plot.file,height=4,width=6,res=150,units='in')

y.ix <- grep(year,format(dates,'%Y'))
ydates <- dates[y.ix]

pr.bccaq.year <- pr.bccaq.series$data[y.ix]
pr.bcbc.year <-  pr.qdm.test[y.ix]
pr.can.year <- (pr.canrcm4.series$data[y.ix])

plot(ydates,pr.bccaq.year,type='l',col='white',lwd=2,ylim=c(0,35),
     col.axis=pcol,col.lab=pcol,xlab='Date',ylab='Daily Pr. (mm)',
     main=paste0('Precip Time Series: ',year),col.main=pcol)
lines(ydates,pr.bcbc.year,col='red',lwd=2)
lines(ydates,pr.can.year,col='green',lwd=1)
legend('topleft',leg=c('BCCAQv2','BCBC','CanRCM4'),col=c('blue','red','green'),pch=18)

box(which='plot',col=pcol)

dev.off()





