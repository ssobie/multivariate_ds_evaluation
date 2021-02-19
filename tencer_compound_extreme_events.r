##Apply the evaluation measures in Francois et al 2020 to
##PNWNAmet and ANUSPLIN as examples 
##

##-----------------------------------------------------------

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##Functions

##------------------------------------------------------------
add_attributes_ncdf <- function(var.name, write.nc, gcm.nc) {

  lon.atts <- ncatt_get(gcm.nc,'lon')
  print('Lon names')
  lon.names <- names(lon.atts)
  for (j in 1:length(lon.atts))
    ncatt_put(write.nc,varid='lon',attname=lon.names[j],attval=lon.atts[[j]])
  print('Lat names')
  lat.atts <- ncatt_get(gcm.nc,'lat')
  lat.names <- names(lat.atts)
  for (j in 1:length(lat.atts))
    ncatt_put(write.nc,varid='lat',attname=lat.names[j],attval=lat.atts[[j]])
  print('Var names')
  ncatt_put(write.nc,varid=var.name,attname='standard_name',attval='spearman_correlation')
  ncatt_put(write.nc,varid=var.name,attname='long_name',attval='Spearman Correlation')
  ncatt_put(write.nc,varid=var.name,attname='missing_value',attval=1.e+20)
  ncatt_put(write.nc,varid=var.name,attname='_FillValue',attval=1.e+20)

  print('Time atts')
  ##Time attributes
  ncatt_put(write.nc,varid='time',attname='units',attval='days since 2000-01-01 00:00:00')
  ncatt_put(write.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(write.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(write.nc,varid='time',attname='calendar','standard')
  print('Global atts')
  ##Global Attributes
  global.atts <- ncatt_get(gcm.nc,varid=0)
  global.names <- names(global.atts)
  for (g in 1:length(global.atts))
    ncatt_put(write.nc,varid=0,attname=global.names[g],attval=global.atts[[g]])
  ##ncatt_put(write.nc,varid=0,attname='blank',attval='Further notes')

  ##Clear extraneous history
  ncatt_put(write.nc,varid=0,attname='history',attval='')
}

##------------------------------------------------------------

get_var_units <- function(var.name) {
   rv <- switch(var.name,
                pr='',
                tasmax='',
                tasmin='',
                sp_cor='')
   return(rv)
}


##-----------------------------------------------------------

make_empty_threshold_file <- function(var.name,input.file,pctl,tmp.dir) {

  thresh.file <- gsub(var.name,paste0(var.name,'_',pctl,'_thresholds_'),input.file)

  work <- paste0("cdo -s daymean ",tmp.dir,input.file," ",tmp.dir,"tmp.file.nc")
  system(work)
  Sys.sleep(1)
  work <- paste0("cdo -s mulc,0 ",tmp.dir,"tmp.file.nc ",tmp.dir,thresh.file)
  system(work)
  Sys.sleep(1)

  nc <- nc_open(paste0(tmp.dir,thresh.file),write=TRUE)

  ncatt_put(nc,var.name,'long_name',paste0(var.name,'_',pctl,'_threshold'))
  ncatt_put(nc,var.name,'standard_name',paste0(var.name,'_',pctl,'_threshold'))

  nc_close(nc)
  file.remove(paste0(tmp.dir,"tmp.file.nc"))

  return(thresh.file)

}

##-----------------------------------------------------------
##Precipitation thresholds

calculate_extreme_precip_thresholds <- function(pr.file,pctl,t.bnds,tmp.dir) {

   nc <- nc_open(paste0(tmp.dir,pr.file))
   scale <- 1
   if (nc$var$pr$units=='kg m-2 s-1') {scale <- 86400}
   
   time <- netcdf.calendar(nc)
   yst <- head(grep(t.bnds[1],time),1)
   yen <- tail(grep(t.bnds[2],time),1)
   yct <- yen-yst+1

   base.dates <- time[yst:yen]
   day.fac <- as.factor(format(time,'%m-%d'))
   ydays <- levels(day.fac)
   ndays <- length(ydays)   

   day.ix <- vector(mode='list',length=length(ydays))

   for (k in 1:length(ydays)) {
      ix <- grep(ydays[k],format(base.dates,'%m-%d'))
      all.ix <- which(time %in% base.dates[ix])
      st <- all.ix-14
      en <- all.ix+14
      day.ix[[k]] <- as.vector((mapply(":",st,en)))
  }
   
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  

  nlon <- length(lon)
  nlat <- length(lat)
   
  thresh.day <- array(NA,c(nlon,nlat,ndays))

   for (i in 1:nlon) {
      print(paste0('Lon: ',i,' of ',nlon))
      pr.1 <- pr <- ncvar_get(nc,'pr',start=c(i,1,1),count=c(1,-1,-1))*scale
      pr.1[pr <= 0.1] <- NA
      for (j in 1:nlat) {       
         for (k in 1:length(ydays)) {
            dx <- day.ix[[k]]
            pr.day <- pr.1[j,dx]
            thresh.day[i,j,k] <- quantile(pr.day,pctl,na.rm=TRUE,names=FALSE)    
         }
      }     
   }
   nc_close(nc)


}

##-----------------------------------------------------------
##Days over thresholds

days_over_under_thresholds <- function() {
}

##-----------------------------------------------------------



##Compare difference in Spearman Correlation between PNWNAmet and ANUSPLIN
##over the entire record, and for two subsets

##-----------------------------------------------------------
clim_file <- function(input.file,output.file,tmp.dir) {

   work <- paste0("cdo timmean ",tmp.dir,input.file," ",tmp.dir,output.file)
   system(work)
   Sys.sleep(1)   
   print('Replace with creation of a new file') 
   return(output.file)
}

##-----------------------------------------------------------

tmp.dir <- "/local_temp/ssobie/tencer_extremes/"

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir)
}

obs <- 'CanRCM4'
t.bnds <- c(1951,1980)
interval <- paste0(t.bnds,collapse="-")

if (obs=='PNWNAmet') {
   obs.dir <- "/storage/data/climate/observations/gridded/PNWNAmet/"
   pr.file <- "pr_allBC_PNWNAmet_1945-2012.nc"
   tasmax.file <- "tasmax_allBC_PNWNAmet_1945-2012.nc"
   tasmin.file <- "tasmin_allBC_PNWNAmet_1945-2012.nc"

}

if (obs=='ANUSPLIN') {
   obs.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/"
   pr.file <- "pr_BC_RCM-Grid_ANUSPLIN_day_19500101-20101231.nc"
   tasmax.file <- "tasmax_BC_RCM-Grid_ANUSPLIN_day_19500101-20101231.nc"
   tasmin.file <- "tasmin_BC_RCM-Grid_ANUSPLIN_day_19500101-20101231.nc"

}

if (obs=='CanRCM4') {
   obs.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"
   pr.file <- "pr_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmax.file <- "tasmax_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmin.file <- "tasmin_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   write.dir <- paste0(obs.dir,"CanRCM4_Derived/tencer_compound/")
}



##-----------------------------------------------
##Precipitation Thresholds

file.copy(from=paste0(obs.dir,pr.file),to=tmp.dir)

pr.thresh.file <- make_empty_threshold_file(var.name='pr',input.file=pr.file,pctl=0.75,tmp.dir)

pr.thresh.values <- calculate_extreme_precip_thresholds(pr.file,pctl=0.75,t.bnds,tmp.dir)

tnc <- nc_open(paste0(tmp.dir,pr.thresh.file),write=TRUE)
ncvar_put(tnc,'pr',pr.thresh.values)
nc_close(tnc)

file.copy(from=paste0(tmp.dir,pr.thresh.file),to=write.dir)

##Precipitation above thresholds

browser()

##-----------------------------------------------





file.copy(from=paste0(obs.dir,tasmax.file),to=tmp.dir)
file.copy(from=paste0(obs.dir,tasmin.file),to=tmp.dir)
print('Done copying')

##Make file to store coefficients

pr.tx.cor.file <- paste0(season,"_spearman_corr_pr_tasmax_BC_",obs,"_",interval,".nc") 
make_empty_spearman_file(var.name='sp_cor',pr.file,pr.tx.cor.file,tmp.dir)
pr.tx.nc <- nc_open(paste0(tmp.dir,pr.tx.cor.file),write=TRUE)

pr.tn.cor.file <- paste0(season,"_spearman_corr_pr_tasmin_BC_",obs,"_",interval,".nc")
make_empty_spearman_file(var.name='sp_cor',pr.file,pr.tn.cor.file,tmp.dir)
pr.tn.nc <- nc_open(paste0(tmp.dir,pr.tn.cor.file),write=TRUE)

pr.ts.cor.file <- paste0(season,"_spearman_corr_pr_tas_BC_",obs,"_",interval,".nc")
make_empty_spearman_file(var.name='sp_cor',pr.file,pr.ts.cor.file,tmp.dir)
pr.ts.nc <- nc_open(paste0(tmp.dir,pr.ts.cor.file),write=TRUE)

tx.tn.cor.file <- paste0(season,"_spearman_corr_tasmax_tasmin_BC_",obs,"_",interval,".nc")
make_empty_spearman_file(var.name='sp_cor',pr.file,tx.tn.cor.file,tmp.dir)
tx.tn.nc <- nc_open(paste0(tmp.dir,tx.tn.cor.file),write=TRUE)

pnc <- nc_open(paste0(tmp.dir,pr.file))
txc <- nc_open(paste0(tmp.dir,tasmax.file))
tnc <- nc_open(paste0(tmp.dir,tasmin.file))

lon <- ncvar_get(pnc,'lon')
nlon <- length(lon)

lat <- ncvar_get(pnc,'lat')
nlat <- length(lat)

time <- netcdf.calendar(pnc)
file.cal <- attr(time,'cal')

seas.ix <- grep(seas.grep,format(time,'%m'))

##Clip to common interval
time.sub <- time[seas.ix]
time.clip <- time.sub <= as.PCICt(paste0(t.bnds[2],'-12-31 23:59:59'),cal=file.cal) &
             time.sub >= as.PCICt(paste0(t.bnds[1],'-01-01 00:00:00'),cal=file.cal) 
time.seas <- seas.ix[time.clip]
ntime <- length(time.seas)

for (i in 1:length(lat)) {
   print(paste0(i," of ",length(lat)))

   pr <- ncvar_get(pnc,'pr',start=c(1,i,1),count=c(-1,1,-1))[,time.seas]
   pr.list <- lapply(seq_len(nrow(pr)), function(k) pr[k,])

   tasmax <- ncvar_get(txc,'tasmax',start=c(1,i,1),count=c(-1,1,-1))[,time.seas]
   tasmax.list <- lapply(seq_len(nrow(tasmax)), function(k) tasmax[k,])

   tasmin <- ncvar_get(tnc,'tasmin',start=c(1,i,1),count=c(-1,1,-1))[,time.seas]
   tasmin.list <- lapply(seq_len(nrow(tasmin)), function(k) tasmin[k,])
   tas <- (tasmax + tasmin) / 2
   tas.list <- lapply(seq_len(nrow(tas)), function(k) tas[k,])

   
   pr.tx.cor <- mapply(FUN=spearman_corr,pr.list,tasmax.list)   
   ncvar_put(pr.tx.nc,'sp_cor',pr.tx.cor,start=c(1,i,1),count=c(nlon,1,1))

   pr.tn.cor <- mapply(FUN=spearman_corr,pr.list,tasmin.list)  
   ncvar_put(pr.tn.nc,'sp_cor',pr.tn.cor,start=c(1,i,1),count=c(nlon,1,1)) 

   pr.ts.cor <- mapply(FUN=spearman_corr,pr.list,tas.list)   
   ncvar_put(pr.ts.nc,'sp_cor',pr.ts.cor,start=c(1,i,1),count=c(nlon,1,1))
   
   tx.tn.cor <- mapply(FUN=spearman_corr,tasmax.list,tasmin.list)   
   ncvar_put(tx.tn.nc,'sp_cor',tx.tn.cor,start=c(1,i,1),count=c(nlon,1,1))

   ##pr.lag1 <- unlist(lapply(pr.list,lag1_autocorrelation))
   ##ncvar_put(pr.lag1.nc,'sp_cor',pr.lag1,start=c(1,i,1),count=c(nlon,1,1))

}

nc_close(pr.tx.nc)
nc_close(pr.tn.nc)
nc_close(pr.ts.nc)
nc_close(tx.tn.nc)

nc_close(pnc)
nc_close(tnc)
nc_close(txc)

file.copy(from=paste0(tmp.dir,pr.tx.cor.file),to=sp.dir,overwrite=TRUE)
file.copy(from=paste0(tmp.dir,pr.tn.cor.file),to=sp.dir,overwrite=TRUE)
file.copy(from=paste0(tmp.dir,pr.ts.cor.file),to=sp.dir,overwrite=TRUE)
file.copy(from=paste0(tmp.dir,tx.tn.cor.file),to=sp.dir,overwrite=TRUE)


