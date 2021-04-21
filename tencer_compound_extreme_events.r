##Apply the evaluation measures in Francois et al 2020 to
##PNWNAmet and ANUSPLIN as examples 
##

##-----------------------------------------------------------

library(ncdf4)
library(PCICt)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##Functions

get_seasonal_fac <- function(seas.dates) {
   seasons <-  c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA", "JJA", "JJA", "SON", "SON", "SON", "DJF")
   years <- as.numeric(format(seas.dates,'%Y'))
   uni.yrs <- unique(years)
   months <- as.numeric(format(seas.dates,'%m'))
   dec.ix <- grep(12,months)
   years[dec.ix] <- years[dec.ix] + 1
   dec.fix <- years %in% uni.yrs
   yearly.fac <- as.factor(years[dec.fix])
   monthly.fac <- as.factor(format(seas.dates[dec.fix],'%m'))
   seasonal.fac <- factor(seasons[monthly.fac], levels=c('DJF', 'MAM', 'JJA', 'SON'))
   avg.fac <- list(yearly.fac,seasonal.fac)

  return(list(fac=avg.fac,fix=dec.fix))
}


##------------------------------------------------------------
add_attributes_ncdf <- function(var.name, write.nc, gcm.nc,pctl,missval=1.e+20,long.name) {

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

  ncatt_put(write.nc,var.name,'long_name',long.name)
  ncatt_put(write.nc,var.name,'standard_name',long.name)
  ncatt_put(write.nc,varid=var.name,attname='missing_value',attval=missval)
  ncatt_put(write.nc,varid=var.name,attname='_FillValue',attval=missval)

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
                pr='kg m-2 s-1',
                tasmax='degC',
                tasmin='degC',
                compound='%')               
   return(rv)
}


##-----------------------------------------------------------

make_empty_threshold_file <- function(var.name,input.file,thresh.file,pctl,tmp.dir) {

  ##thresh.file <- gsub(var.name,paste0(var.name,'_q',round(pctl*100),'_thresholds'),input.file)
  ##thresh.file <- paste0(var.name,'_q',round(pctl*100),'_thresholds_BC_RCM-Grid_",obs,"_day_19500101-21001231.nc")

  nc <- nc_open(paste0(tmp.dir,input.file),write=FALSE)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  time <- netcdf.calendar(nc)

  day.fac <- as.factor(format(time,'%m-%d'))
  ndays <- length(levels(day.fac))  

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', 'days since 2000-01-01 00:00:00', 1:ndays ,
                      unlim=FALSE, calendar=attr(time,'cal'))

  var.geog <- ncvar_def(var.name, units=get_var_units(var.name),
                        dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  new.nc <- nc_create(paste(tmp.dir,thresh.file,sep=''),var.geog)

  print('Time atts')
  ##Time attributes
  ncatt_put(new.nc,varid='time',attname='units',attval='days since 2000-01-01 00:00:00')
  ncatt_put(new.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='calendar','standard')

  add_attributes_ncdf(var.name, new.nc,nc,pctl,long.name=paste0(var.name,'_',pctl,'_threshold'))
  ncvar_put(new.nc,'lon',lon)
  ncvar_put(new.nc,'lat',lat)

  nc_close(nc)
  nc_close(new.nc)

###  return(thresh.file)

}

 
##-----------------------------------------------------------

make_empty_events_file <- function(var.name,input.file,write.file,t.bnds,tmp.dir) {

  ##events.file <- gsub(var.name,paste0(var.name,'_days_over_q',round(pctl*100)),input.file)

  nc <- nc_open(paste0(tmp.dir,input.file),write=FALSE)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')

  time.values <- ncvar_get(nc,'time')
  time.units <- ncatt_get(nc,'time','units')$value  
  time.cal <- ncatt_get(nc,'time','calendar')$value  
  time.series <- netcdf.calendar(nc)

  yst <- head(grep(t.bnds[1],time.series),1)
  yen <- tail(grep(t.bnds[2],time.series),1)

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units,  time.values[yst:yen] ,
                      unlim=FALSE, calendar=time.cal)

  var.geog <- ncvar_def(var.name, units='',
                        dim=list(x.geog, y.geog,t.geog),
                        missval=32767,prec='integer')
  new.nc <- nc_create(paste(tmp.dir,write.file,sep=''),var.geog)

  print('Time atts')
  ##Time attributes
  ncatt_put(new.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(new.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='calendar',time.cal)

  add_attributes_ncdf(var.name, new.nc,nc,pctl,missval=32767,long.name=paste0(var.name,' Days over threshold'))
  ncvar_put(new.nc,'lon',lon)
  ncvar_put(new.nc,'lat',lat)

  nc_close(nc)
  nc_close(new.nc)

}

##-----------------------------------------------------------

make_empty_compound_file <- function(var.pair,input.file,compound.file,type='seasonal',pctl,tmp.dir) {

  nc <- nc_open(paste0(tmp.dir,input.file),write=FALSE)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  time <- netcdf.calendar(nc)
  time.values <- ncvar_get(nc,'time')
  time.units <- ncatt_get(nc,'time','units')$value
  time.cal <-  ncatt_get(nc,'time','calendar')$value

  if (type=='daily') {
    day.fac <- as.factor(format(time,'%m-%d'))
    ndays <- length(levels(day.fac))  
  }
  if (type=='seasonal') {
      seas.grep <- "(*-01-15|*-04-15|*-07-15|*-10-15)"
      seas.ix <- grep(seas.grep,time)
      seas.values <- time.values[seas.ix]
  }
  

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, seas.values ,
                      unlim=FALSE, calendar=time.cal)

  var.geog <- ncvar_def(var.pair, units="count",
                        dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  new.nc <- nc_create(paste(tmp.dir,compound.file,sep=''),var.geog)

  print('Time atts')
  ##Time attributes
  ncatt_put(new.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(new.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='calendar',time.cal)

  add_attributes_ncdf(var.pair, new.nc,nc,pctl,long.name=paste0(var.pair,' Compund Events'))
  ncvar_put(new.nc,'lon',lon)
  ncvar_put(new.nc,'lat',lat)

  nc_close(nc)
  nc_close(new.nc)

  return(compound.file)

}


##-----------------------------------------------------------

make_empty_pr_frequency_file <- function(input.file,freq.file,type='seasonal',tmp.dir) {

  nc <- nc_open(paste0(tmp.dir,input.file),write=FALSE)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  time <- netcdf.calendar(nc)
  time.values <- ncvar_get(nc,'time')
  time.units <- ncatt_get(nc,'time','units')$value
  time.cal <-  ncatt_get(nc,'time','calendar')$value

  if (type=='seasonal') {
      seas.grep <- "(*-01-15|*-04-15|*-07-15|*-10-15)"
      seas.ix <- grep(seas.grep,time)
      seas.values <- time.values[seas.ix]
  }
  

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', time.units, seas.values ,
                      unlim=FALSE, calendar=time.cal)

  var.geog <- ncvar_def('pr', units="%",
                        dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  new.nc <- nc_create(paste(tmp.dir,freq.file,sep=''),var.geog)

  print('Time atts')
  ##Time attributes
  ncatt_put(new.nc,varid='time',attname='units',attval=time.units)
  ncatt_put(new.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='calendar',time.cal)

  add_attributes_ncdf('pr', new.nc,nc,pctl,long.name='Pr Frequency (> 0.1mm)')
  ncvar_put(new.nc,'lon',lon)
  ncvar_put(new.nc,'lat',lat)

  nc_close(nc)
  nc_close(new.nc)

  return(freq.file)

}

##-----------------------------------------------------------

make_empty_ratio_file <- function(var.name,var.pair,input.file,ratio.file,pctl,tmp.dir) {

  ##ratio.file <- gsub(var.name,paste0(var.pair,'_',pctl,'_ratio'),input.file)

  nc <- nc_open(paste0(tmp.dir,input.file),write=FALSE)
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  time <- netcdf.calendar(nc)
  time.values <- grep("(2000-01-15|2000-04-15|2000-07-15|2000-10-15)",time)

  ##Create new netcdf file
  x.geog <- ncdim_def('lon', 'degrees_east', lon)
  y.geog <- ncdim_def('lat', 'degrees_north', lat)
  t.geog <- ncdim_def('time', 'days since 2000-01-01 00:00:00', time.values ,
                      unlim=FALSE, calendar=attr(time,'cal'))

  var.geog <- ncvar_def(var.pair, units="Ratio",
                        dim=list(x.geog, y.geog,t.geog),
                        missval=1.e+20)
  new.nc <- nc_create(paste(tmp.dir,ratio.file,sep=''),var.geog)

  print('Time atts')
  ##Time attributes
  ncatt_put(new.nc,varid='time',attname='units',attval='days since 2000-01-01 00:00:00')
  ncatt_put(new.nc,varid='time',attname='long_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='standard_name',attval='Time')
  ncatt_put(new.nc,varid='time',attname='calendar','standard')

  add_attributes_ncdf(var.pair, new.nc,nc,pctl,long.name=paste0(var.name,'_',pctl,'_ratio'))
  ncvar_put(new.nc,'lon',lon)
  ncvar_put(new.nc,'lat',lat)

  nc_close(nc)
  nc_close(new.nc)

  return(ratio.file)

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
   return(thresh.day)
}

##-----------------------------------------------------------
##Quantile thresholds

calculate_extreme_thresholds <- function(var.name,file,pctl,t.bnds,window=14,tmp.dir) {

   nc <- nc_open(paste0(tmp.dir,file))
   scale <- 1
   if (var.name=='pr') {if (nc$var$pr$units=='kg m-2 s-1') {scale <- 86400}}
   
   time <- netcdf.calendar(nc)
   yst <- head(grep(t.bnds[1],time),1)
   yen <- tail(grep(t.bnds[2],time),1)
   yct <- yen-yst+1

   base.dates <- time[yst:yen]
   day.fac <- as.factor(format(time,'%m-%d'))
   ydays <- levels(day.fac)
   ndays <- length(ydays)   

   day.ix <- vector(mode='list',length=length(ydays))

   ##Window = 14 for precipitation
   ##Window = 2 for temperature

   for (k in 1:length(ydays)) {
      ix <- grep(ydays[k],format(base.dates,'%m-%d'))
      all.ix <- which(time %in% base.dates[ix])
      st <- all.ix-window
      en <- all.ix+window
      day.ix[[k]] <- as.vector((mapply(":",st,en)))
  }
   
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
 

  nlon <- length(lon)
  nlat <- length(lat)
   
  thresh.day <- array(NA,c(nlon,nlat,ndays))

   for (i in 1:nlon) {
      print(paste0('Lon: ',i,' of ',nlon))
      values.1 <- values <- ncvar_get(nc,var.name,start=c(i,1,1),count=c(1,-1,-1))*scale
      if (var.name=='pr') {values.1[values <= 0.1] <- NA}
      for (j in 1:nlat) {       
         for (k in 1:length(ydays)) {
            dx <- day.ix[[k]]
            values.day <- values.1[j,dx]
            thresh.day[i,j,k] <- quantile(values.day,pctl,na.rm=TRUE,names=FALSE)    
         }
      }     
   }
   nc_close(nc)

   return(thresh.day)
}


##-----------------------------------------------------------
##Days over thresholds

days_over_under_thresholds <- function(var.name,var.file,var.thresh.file,var.ext.file,compare,t.bnds,tmp.dir) {

print(compare)
   nc <- nc_open(paste0(tmp.dir,var.file))
   tnc <- nc_open(paste0(tmp.dir,var.thresh.file))
   enc <- nc_open(paste0(tmp.dir,var.ext.file),write=TRUE)

   over.fxn <- function(x,y){return(x>y)}
   under.fxn <- function(x,y){return(x<y)}

   comp.fxn <- switch(compare,
                       over=over.fxn,
                       under=under.fxn)
      

   scale <- 1
   if (var.name=='pr') {if (nc$var$pr$units=='kg m-2 s-1') {scale <- 86400}}
   
   time <- netcdf.calendar(nc)
   yst <- head(grep(t.bnds[1],time),1)
   yen <- tail(grep(t.bnds[2],time),1)
   yct <- yen-yst+1

   time.sub <- time[yst:yen]

   day.fac <- as.factor(format(time.sub,'%m-%d'))
   ydays <- levels(day.fac)
   ndays <- length(ydays)   
   ntime <- length(time.sub)

   day.ix <- vector(mode='list',length=length(ydays))

   for (k in 1:length(ydays)) {
      ix <- grep(ydays[k],format(time.sub,'%m-%d'))
      day.ix[[k]] <- ix
  }
   
  lon <- ncvar_get(nc,'lon')
  lat <- ncvar_get(nc,'lat')
  
  nlon <- length(lon)
  nlat <- length(lat)
   
  thresh.values <- ncvar_get(tnc,var.name)

   for (i in 1:nlon) {
      print(paste0('Lon: ',i,' of ',nlon))
      values <- ncvar_get(nc,var.name,start=c(i,1,yst),count=c(1,-1,yct))*scale
      ext.values <- values*0
      for (j in 1:nlat) {       
         for (k in 1:length(ydays)) {
            day.thresh <- thresh.values[i,j,k]
            dx <- day.ix[[k]]
            ext.values[j,dx] <- comp.fxn(values[j,dx],day.thresh)
         }
      }
      ncvar_put(enc,var.name,ext.values,start=c(i,1,1),count=c(1,nlat,ntime))     
   }

   nc_close(nc)
   nc_close(tnc)
   nc_close(enc)

}

##-----------------------------------------------------------
##Count compound events

count_compound_extremes <- function(var.pair,pr.file,tas.file,compound.file,pctl,tmp.dir) {

   var.split <- strsplit(var.pair,'_')[[1]]
   pr.var <- var.split[1]
   tas.var <- var.split[2]
   
   pnc <- nc_open(paste0(tmp.dir,pr.file))
   tnc <- nc_open(paste0(tmp.dir,tas.file))
   cnc <- nc_open(paste0(tmp.dir,compound.file),write=TRUE)

   time <- netcdf.calendar(tnc)
   seas.fac <- get_seasonal_fac(time)
   seas.time <- netcdf.calendar(cnc)

  lon <- ncvar_get(tnc,'lon')
  lat <- ncvar_get(tnc,'lat')
  
  nlon <- length(lon)
  nlat <- length(lat)
  ntime <- length(seas.time)
     
  thresh.values <- ncvar_get(tnc,tas.var)
  compound.frac <- matrix(NA,nrow=nlat,ncol=length(seas.time))

   for (i in 1:nlon) {
      print(paste0('Lon: ',i,' of ',nlon))
      pr.values <- ncvar_get(pnc,pr.var,start=c(i,1,1),count=c(1,-1,-1))[,seas.fac$fix]
      tas.values <- ncvar_get(tnc,tas.var,start=c(i,1,1),count=c(1,-1,-1))[,seas.fac$fix]
      both.values <- pr.values + tas.values

      for (j in 1:nlat) {       
         ###pr.count <- tapply(both.values[j,],seas.fac$fac,function(x){sum(x==2)/length(x)*100})
         pr.count <- tapply(both.values[j,],seas.fac$fac,function(x){sum(x==2)})
         ###tas.count <- tapply(tas.values[j,],seas.fac$fac,function(x){sum(x)/length(x)})
         compound.frac[j,] <- as.vector(t(round(pr.count,3)))
      }
      ncvar_put(cnc,var.pair,compound.frac,start=c(i,1,1),count=c(1,nlat,ntime))     
   }

   nc_close(pnc)
   nc_close(tnc)
   nc_close(cnc)

}


##-----------------------------------------------------------
##Seasonal Precipitation Frequency

seasonal_precip_frequency <- function(pr.file,freq.file,tmp.dir) {

   pnc <- nc_open(paste0(tmp.dir,pr.file))
   fnc <- nc_open(paste0(tmp.dir,freq.file),write=TRUE)

   time <- netcdf.calendar(pnc)
   seas.fac <- get_seasonal_fac(time)
   seas.time <- netcdf.calendar(fnc)

  lon <- ncvar_get(pnc,'lon')
  lat <- ncvar_get(pnc,'lat')
  
  nlon <- length(lon)
  nlat <- length(lat)
  ntime <- length(seas.time)

   scale <- 1
   if (pnc$var$pr$units=='kg m-2 s-1') {scale <- 86400}
     
  
  pr.freq <- matrix(NA,nrow=nlat,ncol=length(seas.time))

   for (i in 1:nlon) {
      print(paste0('Lon: ',i,' of ',nlon))
      pr.values <- ncvar_get(pnc,'pr',start=c(i,1,1),count=c(1,-1,-1))[,seas.fac$fix]*scale
      for (j in 1:nlat) {       
         frac <- tapply(pr.values[j,],seas.fac$fac,function(x){sum(x>1)/length(x)*100})

         pr.freq[j,] <- as.vector(t(round(frac,3)))
      }
      ncvar_put(fnc,'pr',pr.freq,start=c(i,1,1),count=c(1,nlat,ntime))     
   }

   nc_close(pnc)
   nc_close(fnc)
}

##-----------------------------------------------------------

calculate_compound_event_ratios <- function(compound.file,freq.file,ratio.file,var.pair,t.bnds,tmp.dir,
                                            pr.pctl=0.25,temp.pctl=0.1) {

   cnc <- nc_open(paste0(tmp.dir,compound.file))
   fnc <- nc_open(paste0(tmp.dir,freq.file))
   rnc <- nc_open(paste0(tmp.dir,ratio.file),write=TRUE)
   
   time <- netcdf.calendar(cnc)
   seas.fac <- get_seasonal_fac(time)
   yst <- head(grep(t.bnds[1],time),1)
   yen <- tail(grep(t.bnds[2],time),1)
   yct <- yen-yst+1

   lon <- ncvar_get(cnc,'lon')
   lat <- ncvar_get(cnc,'lat')
  
   nlon <- length(lon)
   nlat <- length(lat)
   ntime <- length(time)
   
  event.ratio <- matrix(NA,nrow=nlat,ncol=4)

  for (i in 1:nlon) {
      print(paste0('Lon: ',i,' of ',nlon))
      cp.values <- ncvar_get(cnc,var.pair,start=c(i,1,yst),count=c(1,-1,yct)) / 100
      freq.values <- ncvar_get(fnc,'pr',start=c(i,1,yst),count=c(1,-1,yct)) / 100
      
      for (j in 1:4) {
         seas.ix <- format(time[yst:yen],'%m') == (c("01","04","07","10")[j])
         cp.clim <- apply(cp.values[,seas.ix],1,mean)
         indep.clim <- apply(freq.values[,seas.ix],1,mean) * pr.pctl * temp.pctl
         event.ratio[,j] <- round(cp.clim / indep.clim,3)
      }

      ncvar_put(rnc,var.pair,event.ratio,start=c(i,1,1),count=c(1,nlat,4))     
  }

   nc_close(cnc)
   nc_close(fnc)
   

}

##-----------------------------------------------------------

make_precip_threshold_file <- function(pr.file,tmp.dir,t.bnds,obs,pctl=0.75) {
   interval <- paste0(t.bnds,collapse="-")
   pr.thresh.file <- paste0("pr_q",round(pctl*100),"_thresholds_BC_RCM-Grid_",obs,"_day_",interval,".nc")
   make_empty_threshold_file(var.name='pr',input.file=pr.file,thresh.file=pr.thresh.file,pctl=0.75,tmp.dir)
   pr.thresh.values <- calculate_extreme_precip_thresholds(pr.file,pctl=0.75,t.bnds,tmp.dir)
   tnc <- nc_open(paste0(tmp.dir,pr.thresh.file),write=TRUE)
   ncvar_put(tnc,'pr',pr.thresh.values)
   nc_close(tnc)
   file.copy(from=paste0(tmp.dir,pr.thresh.file),to=write.dir,overwrite=TRUE)

   return(pr.thresh.file)
}

make_tas_threshold_file <- function(var.name,tas.file,tmp.dir,t.bnds,pctl,obs) {
   interval <- paste0(t.bnds,collapse="-")
   tas.thresh.file <- paste0(var.name,"_q",round(pctl*100),"_thresholds_BC_RCM-Grid_",obs,"_day_",interval,".nc")
   make_empty_threshold_file(var.name=var.name,input.file=tas.file,thresh.file=tas.thresh.file,pctl=pctl,tmp.dir)
   thresh.values <- calculate_extreme_thresholds(var.name=var.name,file=tas.file,
                                                 pctl=pctl,t.bnds=t.bnds,window=2,tmp.dir=tmp.dir)

   tnc <- nc_open(paste0(tmp.dir,tas.thresh.file),write=TRUE)
   ncvar_put(tnc,var.name,thresh.values)
   nc_close(tnc)
   file.copy(from=paste0(tmp.dir,tas.thresh.file),to=write.dir,overwrite=TRUE)
   return(tas.thresh.file)
}


make_compound_count_file <- function(var.name,tas.file,type,tmp.dir,t.bnds,pctl) {

   thresh.file <- make_compund_threshold_file(var.name=var.name,input.file=tas.file,type=type,tmp.dir)

   thresh.values <- calculate_extreme_thresholds(var.name=var.name,file=tas.file,
                                                 pctl=pctl,t.bnds=t.bnds,window=2,tmp.dir=tmp.dir)

   tnc <- nc_open(paste0(tmp.dir,thresh.file),write=TRUE)
   ncvar_put(tnc,var.name,thresh.values)
   nc_close(tnc)
   file.copy(from=paste0(tmp.dir,thresh.file),to=write.dir,overwrite=TRUE)
   return(thresh.file)
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

tmp.dir <- "/local_temp/ssobie/tencer_extremes/"

if (!file.exists(tmp.dir)) {
   dir.create(tmp.dir)
}

obs.list <- c('CanRCM4','BCCAQv2','MBCn')
t.list <- list(c(1951,1980),
               c(1981,2010),
               c(2011,2040),
               c(2041,2070),
               c(2070,2099))

for (obs in obs.list) {
  print(obs)
  for (t.bnds in t.list) {
     interval <- paste0(t.bnds,collapse="-")
     print(interval)

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

if (obs=='BCCAQv2') {
   obs.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"
   pr.file <- "pr_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmax.file <- "tasmax_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmin.file <- "tasmin_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   write.dir <- paste0(obs.dir,"BCCAQv2_Derived/tencer_compound/")
}

if (obs=='MBCn') {
   obs.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"
   pr.file <- "pr_MBCn_iterated_trace_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmax.file <- "tasmax_MBCn_iterated_trace_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   tasmin.file <- "tasmin_MBCn_iterated_trace_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
   write.dir <- paste0(obs.dir,"MBCn_Derived/tencer_compound/")
}

if (!file.exists(write.dir)) {
   dir.create(write.dir,recursive=TRUE)
}



##-----------------------------------------------
##Precipitation Thresholds

if (1==1) {
 
file.copy(from=paste0(obs.dir,pr.file),to=tmp.dir,overwrite=TRUE)
pr.thresh.file <- make_precip_threshold_file(pr.file,tmp.dir,t.bnds,obs)
file.copy(from=paste0(write.dir,pr.thresh.file),to=tmp.dir,overwrite=TRUE)

print('Made threhold file')

##Precipitation above thresholds

ext.pr.file <- paste0("pr_days_over_q75_BC_RCM-Grid_",obs,"_day_",interval,".nc")

make_empty_events_file(var.name='pr',input.file=pr.file,ext.pr.file,t.bnds,tmp.dir)
days_over_under_thresholds('pr',pr.file,pr.thresh.file,ext.pr.file,'over',t.bnds,tmp.dir)
file.copy(from=paste0(tmp.dir,ext.pr.file),to=write.dir,overwrite=TRUE)


}


##-----------------------------------------------
##Temperature Thresholds

if (1==1) {
file.copy(from=paste0(obs.dir,tasmax.file),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(obs.dir,tasmin.file),to=tmp.dir,overwrite=TRUE)

##tasmax.90.file <- make_tas_threshold_file('tasmax',tasmax.file,tmp.dir,t.bnds,pctl=0.9,obs)
##tasmax.10.file <- make_tas_threshold_file('tasmax',tasmax.file,tmp.dir,t.bnds,pctl=0.1,obs)
##tasmin.90.file <- make_tas_threshold_file('tasmin',tasmin.file,tmp.dir,t.bnds,pctl=0.9,obs)
##tasmin.10.file <- make_tas_threshold_file('tasmin',tasmin.file,tmp.dir,t.bnds,pctl=0.1,obs)

var.names <- c('tasmax','tasmin')
pctls <- c(0.9,0.1)

var.names <- 'tasmin'
pctls <- 0.9

for (var.name in var.names) {

   tas.file <- switch(var.name,
                tasmax=tasmax.file,
                tasmin=tasmin.file)
   print(var.name)
   for (pctl in pctls) { 
     print(pctl)
     qval <- paste0("q",round(pctl*100))
     tas.thresh.file <- make_tas_threshold_file(var.name,tas.file,tmp.dir,t.bnds,pctl,obs)

     file.copy(from=paste0(write.dir,tas.thresh.file),to=tmp.dir,overwrite=TRUE)

     ##Temperature above and below thresholds
     compare <- switch(qval,
                       q90='over',
                       q10='under')
      ext.tas.file <- paste0(var.name,"_days_",compare,"_",qval,"_BC_RCM-Grid_",obs,"_day_",interval,".nc")

      make_empty_events_file(var.name=var.name,input.file=tas.file,ext.tas.file,t.bnds,tmp.dir)
      days_over_under_thresholds(var.name,tas.file,tas.thresh.file,ext.tas.file,compare,t.bnds,tmp.dir)
      file.copy(from=paste0(tmp.dir,ext.tas.file),to=write.dir,overwrite=TRUE)

      var.pair <- paste0('pr_',var.name)
      ext.pr.file <- paste0("pr_days_over_q75_BC_RCM-Grid_",obs,"_day_",interval,".nc")
      ext.tas.file <- paste0(var.name,"_days_",compare,"_",qval,"_BC_RCM-Grid_",obs,"_day_",interval,".nc")
      compound.file <- gsub(var.name,paste0(var.pair,'_seasonal_compound'),ext.tas.file)

      file.copy(from=paste0(write.dir,ext.pr.file),to=tmp.dir,overwrite=TRUE)
      file.copy(from=paste0(write.dir,ext.tas.file),to=tmp.dir,overwrite=TRUE)

      make_empty_compound_file(var.pair,ext.tas.file,compound.file,type='seasonal',pctl=qval,tmp.dir)

      count_compound_extremes(var.pair,ext.pr.file,ext.tas.file,compound.file,pctl,tmp.dir) 
      file.copy(from=paste0(tmp.dir,compound.file),to=write.dir,overwrite=TRUE)

      compound.avg.file <- gsub("seasonal_compound","seasonal_average_compound",compound.file)

      work <- paste0("cdo yseasmean ",tmp.dir,compound.file," ",tmp.dir,compound.avg.file) 
      system(work)
      Sys.sleep(1)

      file.copy(from=paste0(tmp.dir,compound.avg.file),to=write.dir,overwrite=TRUE)

      compound.total.file <- gsub("seasonal_compound","seasonal_total_compound",compound.file)

      work <- paste0("cdo yseassum ",tmp.dir,compound.file," ",tmp.dir,compound.total.file) 
      system(work)
      Sys.sleep(1)

      file.copy(from=paste0(tmp.dir,compound.total.file),to=write.dir,overwrite=TRUE)
   }
}

} 
  

##-----------------------------------------------
if (1==0) {
##Count compound extreme events - Moved this section into the loops just above
## cmp.ext <- 100 * n_both / n_temp_ext

var.pair <- 'pr_tasmax'
ext.pr.file <- paste0("pr_days_over_q75_BC_RCM-Grid_",obs,"_day_",interval,".nc")
ext.tas.file <- paste0("tasmax_days_over_q90_BC_RCM-Grid_",obs,"_day_",interval,".nc")
compound.file <- gsub("tasmax",paste0(var.pair,'_seasonal_compound'),ext.tas.file)

file.copy(from=paste0(write.dir,ext.pr.file),to=tmp.dir,overwrite=TRUE)
file.copy(from=paste0(write.dir,ext.tas.file),to=tmp.dir,overwrite=TRUE)

make_empty_compound_file(var.pair,ext.tas.file,compound.file,type='seasonal',pctl='q90',tmp.dir)

count_compound_extremes(var.pair,ext.pr.file,ext.tas.file,compound.file,pctl,tmp.dir) 
file.copy(from=paste0(tmp.dir,compound.file),to=write.dir,overwrite=TRUE)

compound.avg.file <- gsub("seasonal_compound","seasonal_average_compound",compound.file)

work <- paste0("cdo yseasmean ",tmp.dir,compound.file," ",tmp.dir,compound.avg.file) 
system(work)
Sys.sleep(1)

file.copy(from=paste0(tmp.dir,compound.avg.file),to=write.dir,overwrite=TRUE)

compound.total.file <- gsub("seasonal_compound","seasonal_total_compound",compound.file)

work <- paste0("cdo yseassum ",tmp.dir,compound.file," ",tmp.dir,compound.total.file) 
system(work)
Sys.sleep(1)

file.copy(from=paste0(tmp.dir,compound.total.file),to=write.dir,overwrite=TRUE)


}

##-----------------------------------------
##Precipitation Frequency (Pr > 0.1mm) by season
if (1==0) {

freq.file <- paste0("pr_seasonal_frequency_BC_RCM-Grid_",obs,"_day_19500101-21001231.nc")
file.copy(from=paste0(obs.dir,pr.file),to=tmp.dir,overwrite=TRUE)

make_empty_pr_frequency_file(pr.file,freq.file,type='seasonal',tmp.dir)

seasonal_precip_frequency(pr.file,freq.file,tmp.dir) 

file.copy(from=paste0(tmp.dir,freq.file),to=write.dir,overwrite=TRUE)

}


##-----------------------------------------
##Ratio of actual compound event frequency compared to
##expected event frequency

##Subset compound count to interval and calculate climatology
##Subset pr frequency to interval and calculate climatology
##Multiply pr freq by 0.25 (pr) and 0.1 (temp)
##Calculate ratio of compound frequency to expected

if (1==0) {

  var.pair <- 'pr_tasmin'
  pctl <- 'q10'
  compare <- switch(pctl,
             q90='over',
             q10='under')

  compound.file <- paste0(var.pair,"_seasonal_compound_days_",compare,"_",pctl,"_BC_RCM-Grid_",obs,"_day_19500101-21001231.nc")  
  freq.file <- paste0("pr_seasonal_frequency_BC_RCM-Grid_",obs,"_day_19500101-21001231.nc")
  ratio.file <- paste0(var.pair,"_",pctl,"_ratio_BC_RCM-Grid_",obs,"_day_",interval,".nc")

  file.copy(from=paste0(obs.dir,pr.file),to=tmp.dir,overwrite=TRUE) 
  file.copy(from=paste0(write.dir,compound.file),to=tmp.dir,overwrite=TRUE) 
  file.copy(from=paste0(write.dir,freq.file),to=tmp.dir,overwrite=TRUE)
  ratio.file <- make_empty_ratio_file('pr',var.pair,pr.file,ratio.file,pctl,tmp.dir)
  
  calculate_compound_event_ratios(compound.file,freq.file,ratio.file,var.pair,t.bnds,tmp.dir)

  file.copy(from=paste0(tmp.dir,ratio.file),to=write.dir,overwrite=TRUE)  

}

}

}