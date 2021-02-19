##Apply the evaluation measures in Francois et al 2020 to
##PNWNAmet and ANUSPLIN as examples 
##

##-----------------------------------------------------------

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##Functions

##------------------------------------------------------------
add_attributes_ncdf <- function(var.name='sp_cor', write.nc, gcm.nc) {

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

make_empty_spearman_file <- function(var.name='sp_cor',base.file,sp.file,dir) {

  nc <- nc_open(paste0(dir,base.file),write=FALSE)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', 'days since 2000-01-01 00:00:00', 1,
                      unlim=FALSE, calendar='standard')

  var.geog <- ncvar_def(var.name, units=get_var_units(var.name),
                        dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  new.nc <- nc_create(paste(dir,sp.file,sep=''),var.geog)

  add_attributes_ncdf(var.name, new.nc,nc)
  ncvar_put(new.nc,'lon',lon)
  ncvar_put(new.nc,'lat',lat)

  nc_close(nc)
  nc_close(new.nc)
  return(sp.file)

}

##Spearman (rank) correction of two variables

spearman_corr <- function(x,y){

   spc <- cor(x,y,method='spearman')

}

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

tmp.dir <- "/local_temp/ssobie/spear_corr/"

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir)
}

obs <- 'CanRCM4'

seasons <- c('winter','spring','summer','fall','annual')

for (season in seasons) {

t.bnds <- c(1951,1980)
interval <- paste0(t.bnds,collapse="-")


seas.grep <- switch(season,
                    winter='(12|01|02)', ##Winter
                    spring='(03|04|05)', ##Spring
                    summer='(06|07|08)', ##Summer
                    fall='(09|10|11)', ##Fall
                    annual='*') ##Annual

if (obs=='PNWNAmet') {
   obs.dir <- "/storage/data/climate/observations/gridded/PNWNAmet/"
   pr.file <- "pr_allBC_PNWNAmet_1945-2012.nc"
   tasmax.file <- "tasmax_allBC_PNWNAmet_1945-2012.nc"
   tasmin.file <- "tasmin_allBC_PNWNAmet_1945-2012.nc"
   sp.dir <- paste0(obs.dir,'Derived/spearman/')   
   ##pr.file <- "pr_VW_PNWNAmet_1945-2012.nc"
   ##tasmax.file <- "tasmax_VW_PNWNAmet_1945-2012.nc"
   ##tasmin.file <- "tasmin_VW_PNWNAmet_1945-2012.nc"

}

if (obs=='ANUSPLIN') {
   obs.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/"
   pr.file <- "pr_BC_RCM-Grid_ANUSPLIN_day_19500101-20101231.nc"
   tasmax.file <- "tasmax_BC_RCM-Grid_ANUSPLIN_day_19500101-20101231.nc"
   tasmin.file <- "tasmin_BC_RCM-Grid_ANUSPLIN_day_19500101-20101231.nc"
   sp.dir <- paste0(obs.dir,'Derived/spearman/')   
}

if (obs=='CanRCM4') {
   obs.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/"
   pr.file <- "pr_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmax.file <- "tasmax_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmin.file <- "tasmin_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   sp.dir <- paste0(obs.dir,'Derived/spearman/')   
}



file.copy(from=paste0(obs.dir,pr.file),to=tmp.dir)
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

}
