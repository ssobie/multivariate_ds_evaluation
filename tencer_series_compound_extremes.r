##Apply the evaluation measures in Francois et al 2020 to
##PNWNAmet and ANUSPLIN as examples 
##

##-----------------------------------------------------------

library(PCICt)


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


##-----------------------------------------------------------
##Quantile thresholds

calculate_extreme_thresholds <- function(var.name,var.series,pctl,t.bnds,window=14,tmp.dir) {

   ##Window = 14 for precipitation
   ##Window = 2 for temperature

   if (var.name=='pr') {
      window <- 14
   } else {
      window <- 2
   }
   
   time <- var.series$time
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
      st <- all.ix-window
      en <- all.ix+window
      day.ix[[k]] <- as.vector((mapply(":",st,en)))
  }
   
  values.1 <- values <- var.series$data
  thresh.day <- rep(0,length(ydays))

  if (var.name=='pr') {values.1[values <= 0.1] <- NA}

  for (k in 1:length(ydays)) {
     dx <- day.ix[[k]]
     values.day <- values.1[dx]
     thresh.day[k] <- quantile(values.day,pctl,na.rm=TRUE,names=FALSE)    
  }        
  return(thresh.day)
}


##-----------------------------------------------------------
##Days over thresholds

days_over_under_thresholds <- function(var.name,var.series,var.thresh,t.bnds,compare) {

   print(compare)

   over.fxn <- function(x,y){return(x>y)}
   under.fxn <- function(x,y){return(x<y)}

   comp.fxn <- switch(compare,
                       over=over.fxn,
                       under=under.fxn)
   
   time <- var.series$time
   yst <- head(grep(t.bnds[1],time),1)
   yen <- tail(grep(t.bnds[2],time),1)
   
   time.sub <- time[yst:yen]
   series.sub <- var.series$data[yst:yen]

   day.fac <- as.factor(format(time.sub,'%m-%d'))
   ydays <- levels(day.fac)
   ndays <- length(ydays)   
   ntime <- length(time.sub)

   day.ix <- vector(mode='list',length=length(ydays))

   for (k in 1:length(ydays)) {
      ix <- grep(ydays[k],format(time.sub,'%m-%d'))
      day.ix[[k]] <- ix
  }
   
  thresh.values <- var.thresh

  ext.values <- series.sub*0
  for (k in 1:length(ydays)) {
     day.thresh <- thresh.values[k]
     dx <- day.ix[[k]]
     ext.values[dx] <- comp.fxn(series.sub[dx],day.thresh)
   }
   return(ext.values)

}

##-----------------------------------------------------------
##Count compound events

count_compound_extremes <- function(var.pair,var.pr,var.tas,t.bnds,time) {

   var.split <- strsplit(var.pair,'_')[[1]]
   pr.var <- var.split[1]
   tas.var <- var.split[2]  

   yst <- head(grep(t.bnds[1],time),1)
   yen <- tail(grep(t.bnds[2],time),1)
   time.sub <- time[yst:yen]

   seas.fac <- get_seasonal_fac(time.sub)

   pr.values <- var.pr[seas.fac$fix]
   tas.values <- var.tas[seas.fac$fix]
   both.values <- pr.values + tas.values

   pr.count <- tapply(both.values,seas.fac$fac,function(x){sum(x==2)/length(x)*100})
   pr.sum <- apply(tapply(both.values,seas.fac$fac,function(x){sum(x==2)}),2,sum)
   ###compound.frac <- as.vector(round(t(pr.count),3))
   rv <- list(total=pr.sum,frac=round(pr.count,3))
   return(rv)
}


##-----------------------------------------------------------
##Seasonal Precipitation Frequency

seasonal_precip_frequency <- function(pr.series,t.bnds) {

   time <- pr.series$time
   yst <- head(grep(t.bnds[1],time),1)
   yen <- tail(grep(t.bnds[2],time),1)
   time.sub <- time[yst:yen]
   pr.sub <- pr.series$data[yst:yen]

   seas.fac <- get_seasonal_fac(time.sub)

   frac <- tapply(pr.sub[seas.fac$fix],seas.fac$fac,function(x){sum(x>0.1)/length(x)*100})
   ###pr.freq <- as.vector(t(round(frac,3)))
   return(round(frac,3))

}

##-----------------------------------------------------------

calculate_compound_event_ratios <- function(var.pair,compound.frac,pr.freq,
                                            pr.pctl=0.25,temp.pctl=0.1) {

   cp.values <- compound.frac / 100
   freq.values <- pr.freq / 100
      
   for (j in 1:4) {
      seas.ix <- format(time[yst:yen],'%m') == (c("01","04","07","10")[j])
      cp.clim <- apply(cp.values[,seas.ix],1,mean)
      indep.clim <- apply(freq.values[,seas.ix],1,mean) * pr.pctl * temp.pctl
      event.ratio[,j] <- round(cp.clim / indep.clim,3)
   }

}

##-----------------------------------------------------------

pr_days_count <- function(all.series,type,t.bnds) {

   ##-----------------------------------------------
   ##Precipitation Thresholds - need to loop over intervals
         
   pr.q75 <-  calculate_extreme_thresholds(var.name='pr',var.series=all.series$pr[[type]],
                                            pctl=0.75,t.bnds=t.bnds)                                           

   pr.days <- days_over_under_thresholds(var.name='pr',var.series=all.series$pr[[type]],pr.q75,
                                         t.bnds=t.bnds,'over')
   rv <- list(days=pr.days,thresh=pr.q75)                                         
   return(rv) ##pr.days)  
}

temp_days_count <- function(all.series,type,var.name,pctl,t.bnds) {

   qval <- paste0('q',pctl)                

   compare <- switch(qval,
                     q90='over',
                     q10='under')

   temp.qt <-  calculate_extreme_thresholds(var.name=var.name,var.series=all.series[[var.name]][[type]],
                                            pctl=pctl/100,t.bnds=t.bnds)            
   temp.days <- days_over_under_thresholds(var.name=var.name,var.series=all.series[[var.name]][[type]],temp.qt,
                                         t.bnds=t.bnds,compare)

   rv <- list(days=temp.days,thresh=temp.qt)
   return(rv)                                         
}




##-----------------------------------------------------------

##****************************************************************

##interval <- paste0(t.bnds,collapse="-")
var.list <- c('pr','tasmax','tasmin','tas')

###
sites <- list(list(site.name="Fort_Nelson",lonc=-122.71,latc=58.79),
              list(site.name="Skagway",lonc=-135.0,latc=59.25),
              list(site.name="Prince_Rupert",lonc=-130.3,latc=54.3),
              list(site.name="Golden",lonc=-116.99,latc=51.3),
              list(site.name="Ucluelet",lonc=-125.5,latc=48.9),
              list(site.name="Prince_George",lonc=-122.75,latc=53.92),
              list(site.name="Kelowna",lonc=-119.5,latc=49.9),
              list(site.name="Victoria",lonc=-123.4,latc=48.6))


for (site in sites) {

   site.name <- site$site.name
   print(site.name)

   load.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/save_dir/",site.name,"/")
   plot.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/plots/",site.name,"/")
   save.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/RData/",site.name,"/")
   if (!file.exists(plot.dir)) {
      dir.create(plot.dir,recursive=TRUE)
   }
   if (!file.exists(save.dir)) {
      dir.create(save.dir,recursive=TRUE)
   }

   all.series <- vector(mode='list',length=4)
   names(all.series) <- var.list
   for (var.name in var.list) {

      ##Precipitation Series
      type.series <- vector(mode='list',length=3)
      names(type.series) <- c('rcm','bccaqv2','mbcn') ##,'bcci',

      load.file <- paste0(var.name,"_day_CanRCM4_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$rcm <- var.series
      rm(var.series)
      ##load.file <- paste0(var.name,"_day_BCCI_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      ##load(paste0(load.dir,load.file))
      ##type.series$bcci <- var.series
      load.file <- paste0(var.name,"_day_BCCAQv2_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$bccaqv2 <- var.series
      rm(var.series)
      load.file <- paste0(var.name,"_day_MBCn_iterated_trace_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$mbcn <- var.series
      all.series[[var.name]] <- type.series
      rm(type.series)
   }

   types <- c('rcm','bccaqv2','mbcn') ##,'bcci'
   intervals <- list(c(1951,1980),
                     c(1981,2010),
                     c(2011,2040),
                     c(2041,2070),
                     c(2071,2100))

   temp.pairs <- list(tx90=list(temp='tasmax',pctl=90),
                      tn10=list(temp='tasmin',pctl=10))
                      ##tx10=list(temp='tasmax',pctl=10),
                      ##tn90=list(temp='tasmin',pctl=10),
                      

   ##Iterate over types, intervals, variables pairs  
   season.mat <- matrix(NA,nrow=length(types),ncol=length(intervals))
   compounds <- list(winter=season.mat,spring=season.mat,
                     summer=season.mat,fall=season.mat)

   compound.counts <- vector(mode='list',length=length(temp.pairs))                         
   pr.thresholds <- vector(mode='list',length=length(types))
   temp.thresholds <- vector(mode='list',length=length(types))
   tx.thresholds <- vector(mode='list',length=length(types))
   tn.thresholds <- vector(mode='list',length=length(types))

   for (k in 1:2) {
     
      temp.pair <- temp.pairs[[k]]
      var.pair <- paste0("pr_",temp.pair$temp)
      print(paste0(var.pair,' q',temp.pair$pctl))

      for (i in seq_along(types)) {
         type <- types[i]
         print(type)
         pr.thresh <- temp.thresh <- matrix(NA,nrow=length(intervals),ncol=365)
         for (j in seq_along(intervals)) {
            t.bnds <- intervals[[j]]
            print(paste0(t.bnds,collapse='-'))
            pr.q.days <- pr_days_count(all.series,type,t.bnds)
            pr.days <- pr.q.days$days
         
            temp.q.days <- temp_days_count(all.series,type,
                                      var.name=temp.pair$temp,pctl=temp.pair$pctl,t.bnds)
            temp.days <- temp.q.days$days                                      

            pr.temp.compound <- count_compound_extremes(var.pair=var.pair,pr.days,temp.days,
                                                     t.bnds=t.bnds,all.series$pr[[type]]$time)
            compounds$winter[i,j] <- pr.temp.compound$total[1]                                                     
            compounds$spring[i,j] <- pr.temp.compound$total[2]  
            compounds$summer[i,j] <- pr.temp.compound$total[3]  
            compounds$fall[i,j] <- pr.temp.compound$total[4]              

            pr.thresh[j,] <- pr.q.days$thresh
            temp.thresh[j,] <- temp.q.days$thresh
         }
         pr.thresholds[[i]] <- pr.thresh
         temp.thresholds[[i]] <- temp.thresh
      }                                                  
      if (k==1) { tx.thresholds <- temp.thresholds }
      if (k==2) { tn.thresholds <- temp.thresholds }
      compound.counts[[k]] <- compounds
      names(compound.counts)[k] <- names(temp.pairs)[k]
   }
   save(compound.counts,file=paste0(save.dir,"tencer_seasonal_counts_",site.name,".trace.RData"))
   names(pr.thresholds) <- types 
   names(tx.thresholds) <- types  
   names(tn.thresholds) <- types 

   thresholds <- list(pr=pr.thresholds,
                      tx=tx.thresholds,
                      tn=tn.thresholds) 
   save(thresholds,file=paste0(save.dir,"tencer_thresholds_",site.name,".trace.RData"))

   if (1==0) {
   ##-----------------------------------------------
   ##Precipitation Thresholds - need to loop over intervals
         
   pr.q75 <-  calculate_extreme_thresholds(var.name='pr',var.series=all.series$pr$rcm,
                                            pctl=0.75,t.bnds=c(1951,1980))                                           
   pr.days <- days_over_under_thresholds(var.name='pr',var.series=all.series$pr$rcm,pr.q75,
                                         t.bnds=c(1951,1980),'over')
   
   tx.q90 <-  calculate_extreme_thresholds(var.name='tasmax',var.series=all.series$tasmax$rcm,
                                            pctl=0.90,t.bnds=c(1951,1980))            
   tx.days <- days_over_under_thresholds(var.name='tasmax',var.series=all.series$tasmax$rcm,tx.q90,
                                         t.bnds=c(1951,1980),'over')

   pr.tx.compound <- count_compound_extremes(var.pair="pr_tasmax",pr.days,tx.days,t.bnds=c(1951,1980),all.series$pr$rcm$time) 
   ##Use the count of events over each interval: pr.tx.compound$total
   ##pr.freq <- seasonal_precip_frequency(all.series$pr$rcm,t.bnds=c(1951,1980))
   }


}




