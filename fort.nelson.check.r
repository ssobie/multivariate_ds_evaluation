##Compare series at one cell for all the different datasets

library(scales)
library(ncdf4)
library(PCICt)
library(xtable)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----------------------------------------------------
##Comparisons to consider

##Consider stats over 1951-1980, 1981-2010, 2011-2040, 
##                      2041-2070, 2071-2100

##----------------------------------------------------

get_seasonal_factor <- function(seas.dates) {

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



##----------------------------------------------------
##Degree Day Values
dd <- function(tmean,tbase) {
  g <- tmean - tbase
  days <- sum(g[g>0], na.rm=T)
  return(round(days))
}

gdd<-function(data,fac){tapply(data,fac, dd, tbase=5)}   ##Growing degree days
cdd<-function(data,fac){tapply(data,fac, dd, tbase=18)}  ##Cooling degree days
hdd<-function(data,fac){tapply(-data,fac,dd, tbase=-18)} ##Heating degree days
fdd<-function(data,fac){tapply(-data,fac,dd, tbase=0)} ##Freezing degree days


##----------------------------------------------------
##Precipitation as Snow
pas_fxn <- function(pr,tas) {
   temp <- tas > 0
   pas <- pr
   pas[temp] <- 0
   return(pas)
}

compute_pas <- function(input.series) {

   rcm.pas <- pas_fxn(input.series$pr$rcm$data,input.series$tas$rcm$data)
   bcci.pas <- pas_fxn(input.series$pr$bcci$data,input.series$tas$bcci$data)
   bccaqv2.pas <- pas_fxn(input.series$pr$bccaqv2$data,input.series$tas$bccaqv2$data)
   mbcn.pas <- pas_fxn(input.series$pr$mbcn$data,input.series$tas$mbcn$data)
   rv <- list(rcm=list(time=input.series$pr$rcm$time,
                       data=rcm.pas),
              bcci=list(time=input.series$pr$bcci$time,
                       data=bcci.pas),
              bccaqv2=list(time=input.series$pr$bccaqv2$time,     
                       data=bccaqv2.pas),
              mbcn=list(time=input.series$pr$mbcn$time,
                       data=mbcn.pas))     
   return(rv)
}


##----------------------------------------------------
##Precipitation as Snow
ftd_fxn <- function(tasmax,tasmin) {
   ftd <- tasmax > 0 & tasmin < 0
   return(ftd)
}

compute_ftd <- function(input.series) {

   rcm.ftd <- ftd_fxn(input.series$tasmax$rcm$data,input.series$tasmin$rcm$data)
   bcci.ftd <- ftd_fxn(input.series$tasmax$bcci$data,input.series$tasmin$bcci$data)
   bccaqv2.ftd <- ftd_fxn(input.series$tasmax$bccaqv2$data,input.series$tasmin$bccaqv2$data)
   mbcn.ftd <- ftd_fxn(input.series$tasmax$mbcn$data,input.series$tasmin$mbcn$data)

   rv <- list(rcm=list(time=input.series$tasmax$rcm$time,
                       data=rcm.ftd),
              bcci=list(time=input.series$tasmax$bcci$time,
                       data=bcci.ftd),
              bccaqv2=list(time=input.series$tasmax$bccaqv2$time,     
                       data=bccaqv2.ftd),
              mbcn=list(time=input.series$tasmax$mbcn$time,
                       data=mbcn.ftd))     
   return(rv)
}

##----------------------------------------------------

make_time_series_plot <- function(x.vals,rcm.series,bcci.series,bccaqv2.series,mbcn.series,
                                  plot.style) {
                                  
   y.vals <- c(rcm.series,bcci.series,bccaqv2.series,mbcn.series)

   plot.col <- plot.style$bg
   cx <- plot.style$cex
   cols <- plot.style$cols

   plot(c(),xlim=c(1949,2101),ylim=range(pretty(y.vals)),axes=F,
        main="",xlab='',ylab='',cex=cx,cex.lab=cx,xaxs='i',
        col.axis=plot.col,col.main=plot.col,col.lab=plot.col)

   axis(1,at=pretty(x.vals),label=pretty(x.vals),cex.axis=cx,col.axis=plot.col)
   axis(2,at=pretty(y.vals),label=pretty(y.vals),cex.axis=cx,col.axis=plot.col)
   abline(h=pretty(y.vals),col='lightgray',lty=2,lwd=0.7)
   abline(v=seq(1950,2100,30),col=plot.col,lwd=0.7)
   lines(x.vals,bcci.series,lwd=2.3,col=cols$bcci)     
   lines(x.vals,rcm.series,lwd=2.0,col=cols$rcm)     
   lines(x.vals,bccaqv2.series,lwd=1.7,col=cols$bccaqv2)     
   lines(x.vals,mbcn.series,lwd=1.4,col=cols$mbcn)

   legend('topleft',leg=c('CI','CanRCM4','BCCAQv2','MBCn'),col=unlist(cols),box.col=plot.col,
          pt.cex=2.5,pch=15,cex=cx,text.col=plot.col,horiz=TRUE)
   box(which='plot',col=plot.col)

}

##----------------------------------------------------
##Boxplots

make_interval_boxplots <- function(x.vals,rcm.series,bcci.series,bccaqv2.series,mbcn.series,
                                   plot.style) {

   plot.col <- plot.style$bg
   cx <- plot.style$cex
   cols <- plot.style$cols

   intervals <- list(c(1951,1980),
                     c(1981,2010),
                     c(2011,2040),
                     c(2041,2070),
                     c(2071,2100))
   rcm.mean <- bccaqv2.mean <- mbcn.mean <- rep(0,length(intervals))
   rcm.sd <- bccaqv2.sd <- mbcn.sd <- rep(0,length(intervals))
   bccaqv2.bias <- mbcn.bias <- rep(0,length(intervals))
   bccaqv2.rmse <- mbcn.rmse <- rep(0,length(intervals))
                           
   for (j in 1:5) {
      bnds <- intervals[[j]]
      yst <- grep(bnds[1],x.vals)
      yen <- grep(bnds[2],x.vals)
      bcci.sub <- bcci.series[yst:yen]
      rcm.sub <- rcm.series[yst:yen]
      bccaqv2.sub <- bccaqv2.series[yst:yen]
      mbcn.sub <- mbcn.series[yst:yen]
    
      y.vals <- c(rcm.series,bcci.series,bccaqv2.series,mbcn.series)
      y.lim <- range(y.vals)
      y.lim[2] <- y.lim[2] + 0.05*y.lim[2]
      plot(c(),xlim=c(0.5,4.5),ylim=y.lim,main="",
               xlab='',ylab='',cex=cx,axes=F,xaxs='i',
               col.axis=plot.col,col.main=plot.col,col.lab=plot.col)
      boxplot(bcci.sub,at=1,border=cols$bcci,col=alpha(cols$bcci,0.3),add=T,axes=F)             
      boxplot(rcm.sub,at=2,border=cols$rcm,col=alpha(cols$rcm,0.3),add=T,axes=F)             
      boxplot(bccaqv2.sub,at=3,border=cols$bccaqv2,col=alpha(cols$bccaqv2,0.3),add=T,axes=F)             
      boxplot(mbcn.sub,at=4,border=cols$mbcn,col=alpha(cols$mbcn,0.3),add=T,axes=F)             
      axis(1,at=1:4,label=c('CI','CanRCM4','BCCAQv2','MBCn'),cex.axis=cx,col.axis=plot.col,las=2)
      if (j==1) {
        axis(2,at=pretty(y.vals),label=pretty(y.vals),cex.axis=cx,col.axis=plot.col)
      }
      box(which='plot',col=plot.col)
      text(2.5, 0.975*(y.lim[2]), paste0(bnds,collapse='-'),col=plot.col,cex=cx)

      ##Stats for each interval
      ##Mean, SD, Bias, RMSE

      rcm.mean[j] <- mean(rcm.sub)
      bccaqv2.mean[j] <- mean(bccaqv2.sub)
      mbcn.mean[j] <- mean(mbcn.sub)

      rcm.sd[j] <- sd(rcm.sub)
      bccaqv2.sd[j] <- sd(bccaqv2.sub)
      mbcn.sd[j] <- sd(mbcn.sub)

      ##Bias
      bccaqv2.bias[j] <- mean(bccaqv2.sub - rcm.sub)
      mbcn.bias[j] <- mean(mbcn.sub - rcm.sub)

      ##RMSE 
      bccaqv2.rmse[j] <- sqrt( mean((bccaqv2.sub - rcm.sub)^2))
      mbcn.rmse[j] <- sqrt( mean((mbcn.sub - rcm.sub)^2))

   }

   rv <- list(rcm=list(avg=rcm.mean,std=rcm.sd),
              bccaqv2=list(avg=bccaqv2.mean,std=bccaqv2.sd,bias=bccaqv2.bias,rmse=bccaqv2.rmse),
              mbcn=list(avg=mbcn.mean,std=mbcn.sd,bias=mbcn.bias,rmse=mbcn.rmse))

   return(rv)
}

##----------------------------------------------------

annual_quantity <- function(input.series,ann.fx) {

   dates <- input.series$rcm$time
   yr.fac <- as.factor(format(dates,'%Y'))
   x.vals <- as.numeric(levels(yr.fac))

   rcm.ann <- tapply(input.series$rcm$data,yr.fac,ann.fx)
   bcci.ann <- tapply(input.series$bcci$data,yr.fac,ann.fx)
   bccaqv2.ann <- tapply(input.series$bccaqv2$data,yr.fac,ann.fx)
   mbcn.ann <- tapply(input.series$mbcn$data,yr.fac,ann.fx)
   rv <- list(rcm=rcm.ann,bcci=bcci.ann,bccaqv2=bccaqv2.ann,mbcn=mbcn.ann,x=x.vals)

   return(rv)

}

##----------------------------------------------------


seasonal_quantity <- function(input.series,seas.fx,seas) {

   dates <- input.series$rcm$time      
   yr.fac <- as.factor(format(dates,'%Y'))
   x.vals <- as.numeric(levels(yr.fac))

   seas.factors <- get_seasonal_factor(dates)  
   s.fac <- seas.factors$fac
   s.fix <- seas.factors$fix
   
   rcm.seas <- tapply(input.series$rcm$data[s.fix],s.fac,seas.fx)
   bcci.seas <- tapply(input.series$bcci$data[s.fix],s.fac,seas.fx)
   bccaqv2.seas <- tapply(input.series$bccaqv2$data[s.fix],s.fac,seas.fx)
   mbcn.seas <- tapply(input.series$mbcn$data[s.fix],s.fac,seas.fx)

   seas.ix <- switch(seas,
                     Winter=1,Spring=2,Summer=3,Fall=4)

   rcm.ann <- rcm.seas[,seas.ix]
   bcci.ann <- bcci.seas[,seas.ix]
   bccaqv2.ann <- bccaqv2.seas[,seas.ix]
   mbcn.ann <- mbcn.seas[,seas.ix]

   rv <- list(rcm=rcm.ann,bcci=bcci.ann,bccaqv2=bccaqv2.ann,mbcn=mbcn.ann,x=x.vals)

   return(rv)

}

##----------------------------------------------------

dd_quantity <- function(input.series,dd.fx) {

   dates <- input.series$rcm$time
   yr.fac <- as.factor(format(dates,'%Y'))
   x.vals <- as.numeric(levels(yr.fac))

   rcm.ann <- dd.fx(input.series$rcm$data,yr.fac)
   bcci.ann <- dd.fx(input.series$bcci$data,yr.fac)
   bccaqv2.ann <- dd.fx(input.series$bccaqv2$data,yr.fac)
   mbcn.ann <- dd.fx(input.series$mbcn$data,yr.fac)

   rv <- list(rcm=rcm.ann,bcci=bcci.ann,bccaqv2=bccaqv2.ann,mbcn=mbcn.ann,x=x.vals)

   return(rv)

}

##----------------------------------------------------

annual_comparison <- function(input.series,
                              site.name,plot.file,
                              plot.title,y.lab) {

   rcm.ann <- input.series$rcm
   bcci.ann <- input.series$bcci
   bccaqv2.ann <- input.series$bccaqv2
   mbcn.ann <- input.series$mbcn

   ##Make time series plot above and boxplots with listed stats below for 
   ##each interval

   x.vals <- input.series$x

   cols <- list(bcci='blue',rcm='green',bccaqv2='orange',mbcn='red')
   plot.col='gray22'
   cx <- 1.5
   plot.style <- list(bg=plot.col,cex=cx,cols=cols)
   
   png(file=plot.file,width=5,height=3.5,units='in',res=250,pointsize=6,bg='white')
   par(oma=c(4,6,4,2))
   par(mar=c(4,0,0,0))   
   layout(rbind(c(1,1,1,1,1),2:6))

   make_time_series_plot(x.vals,rcm.ann,bcci.ann,bccaqv2.ann,mbcn.ann,
                        plot.style)
   i.stats <- make_interval_boxplots(x.vals,rcm.ann,bcci.ann,bccaqv2.ann,mbcn.ann,
                          plot.style) 
   
   mtext(plot.title,side=3,line=1,col=plot.col,outer=TRUE)
   mtext(y.lab,side=2,line=3,col=plot.col,outer=TRUE)
   dev.off()     

   return(i.stats)
}

##----------------------------------------------------



##****************************************************
var.list <- c('pr','tasmax','tasmin','tas')

sites <- list(list(site.name="Fort_Nelson",lonc=-122.71,latc=58.79),
              list(site.name="Skagway",lonc=-135.0,latc=59.25),
              list(site.name="Prince_George",lonc=-122.75,latc=53.92),
              list(site.name="Prince_Rupert",lonc=-130.3,latc=54.3),
              list(site.name="Golden",lonc=-116.99,latc=51.3),
              list(site.name="Ucluelet",lonc=-125.5,latc=48.9))

##sites <- list(list(site.name="Prince_George",lonc=-122.75,latc=53.92),
##              list(site.name="Kelowna",lonc=-119.5,latc=49.9),
##              list(site.name="Victoria",lonc=-123.4,latc=48.4))

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
      type.series <- vector(mode='list',length=4)
      names(type.series) <- c('rcm','bcci','bccaqv2','mbcn')

      load.file <- paste0(var.name,"_day_CanRCM4_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$rcm <- var.series
      load.file <- paste0(var.name,"_day_BCCI_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$bcci <- var.series
      load.file <- paste0(var.name,"_day_BCCAQv2_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$bccaqv2 <- var.series
      load.file <- paste0(var.name,"_day_MBCn_iterated_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$mbcn <- var.series
      all.series[[var.name]] <- type.series
   }


##Plot average daily precipitation  

   r.time <- all.series$pr$rcm$time
   years <- format(r.time,'%Y')
   bnds <- c(1951,1980)
   yst <- head(grep(bnds[1],years),1)
   yen <- tail(grep(bnds[2],years),1)
   time.sub <- r.time[yst:yen]
   jdays <- as.factor(format(time.sub,'%j'))
      
   rcm.sub <- all.series$pr$rcm$data[yst:yen]
   mbc.sub <- all.series$pr$mbcn$data[yst:yen]

   rcm.sub[rcm.sub < 0.1] <- NA
   mbc.sub[mbc.sub < 0.1] <- NA

   rcm.days <- tapply(rcm.sub,jdays,function(x){sum(!is.na(x))})
   mbc.days <- tapply(mbc.sub,jdays,function(x){sum(!is.na(x))})

   plot(rcm.days,type='l',col='green',lwd=3) 
   lines(mbc.days,col='red',lwd=2)
   box(which='plot')   
browser()

   pr.ann.file <- paste0(plot.dir,'annual_mean_pr_',site.name,'.png')   
   pr.ann.title <- paste0(gsub('_',' ',site.name)," Annual Mean Precipitation")
   pr.ann.y <- "Precipitation (mm)"   
   pr.ann.series <- annual_quantity(all.series$pr,ann.fx=mean)
   pr.ann <- annual_comparison(pr.ann.series,site.name,pr.ann.file,pr.ann.title,pr.ann.y)
   save(pr.ann,file=paste0(save.dir,'annual_mean_pr_',site.name,'.RData'))

if (1==1) {
   ##Precipitation Indices
   pr.ann.file <- paste0(plot.dir,'annual_total_pr_',site.name,'.png')   
   pr.ann.title <- paste0(gsub('_',' ',site.name)," Annual Total Precipitation")
   pr.ann.y <- "Precipitation (mm)"   
   pr.ann.series <- annual_quantity(all.series$pr,ann.fx=sum)
   pr.ann <- annual_comparison(pr.ann.series,site.name,pr.ann.file,pr.ann.title,pr.ann.y)
   save(pr.ann,file=paste0(save.dir,'annual_total_pr_',site.name,'.RData'))

   pr.max.file <- paste0(plot.dir,'annual_max_pr_',site.name,'.png')   
   pr.max.title <- paste0(gsub('_',' ',site.name)," Annual Maximum Precipitation")
   pr.max.y <- "Precipitation (mm)"   
   pr.max.series <- annual_quantity(all.series$pr,ann.fx=max)
   pr.max <- annual_comparison(pr.max.series,site.name,pr.max.file,pr.max.title,pr.max.y)
   save(pr.max,file=paste0(save.dir,'annual_maximum_pr_',site.name,'.RData'))

   ##Seasonal Total Precipitation 
   seasons <- c('Winter','Spring','Summer','Fall')
   pr.seasonal <- vector(mode='list',length=4)
   for (s in seq_along(seasons)) {
      season <- seasons[s]
      pr.seas.file <- paste0(plot.dir,'annual_',tolower(season),'_pr_',site.name,'.png')   
      pr.seas.title <- paste0(gsub('_',' ',site.name)," ",season," Total Precipitation")
      pr.seas.y <- "Precipitation (mm)"   
      pr.seas.total <- seasonal_quantity(all.series$pr,seas.fx=sum,seas=season)
      pr.seasonal[[s]] <- annual_comparison(pr.seas.total,site.name,pr.seas.file,pr.seas.title,pr.seas.y)
   }
   save(pr.seasonal,file=paste0(save.dir,'seasonal_total_pr_',site.name,'.RData'))


   ##------------------------------------------------------------------------------------------------------
   ##Maximum Temperature Indices
   tx.ann.file <- paste0(plot.dir,'annual_average_tasmax_',site.name,'.png')   
   tx.ann.title <- paste0(gsub('_',' ',site.name)," Annual Average Maximum Temperature")
   tx.ann.y <- "Maximum Temperature (\u00B0C)"   
   tx.ann.series <- annual_quantity(all.series$tasmax,ann.fx=mean)
   tx.ann <- annual_comparison(tx.ann.series,site.name,tx.ann.file,tx.ann.title,tx.ann.y)
   save(tx.ann,file=paste0(save.dir,'annual_average_tasmax_',site.name,'.RData'))

   tx.max.file <- paste0(plot.dir,'annual_max_tasmax_',site.name,'.png')   
   tx.max.title <- paste0(gsub('_',' ',site.name)," Annual Maximum Temperature")
   tx.max.y <- "Maximum Temperature (\u00B0C)"   
   tx.max.series <- annual_quantity(all.series$tasmax,ann.fx=max)
   tx.max <- annual_comparison(tx.max.series,site.name,tx.max.file,tx.max.title,tx.max.y)
   save(tx.max,file=paste0(save.dir,'annual_maximum_tasmax_',site.name,'.RData'))

   ##Seasonal Maximum 
   seasons <- c('Winter','Spring','Summer','Fall')
   tasmax.seasonal <- vector(mode='list',length=4)
   names(tasmax.seasonal) <- seasons 
   for (s in seq_along(seasons)) {
      season <- seasons[s]
      tx.seas.file <- paste0(plot.dir,'annual_',tolower(season),'_tasmax_',site.name,'.png')   
      tx.seas.title <- paste0(gsub('_',' ',site.name)," ",season," Maximum Temperature")
      tx.seas.y <- "Maximum Temperature (\u00B0C)"
      tx.seas.mean <- seasonal_quantity(all.series$tasmax,seas.fx=mean,seas=season)
      tasmax.seasonal[[s]] <- annual_comparison(tx.seas.mean,site.name,tx.seas.file,tx.seas.title,tx.seas.y)
   }
   save(tasmax.seasonal,file=paste0(save.dir,'seasonal_average_tasmax_',site.name,'.RData'))

   ##Minimum Temperature Indices
   tn.ann.file <- paste0(plot.dir,'annual_average_tasmin_',site.name,'.png')   
   tn.ann.title <- paste0(gsub('_',' ',site.name)," Annual Average Minimum Temperature")
   tn.ann.y <- "Minimum Temperature (\u00B0C)"   
   tn.ann.series <- annual_quantity(all.series$tasmin,ann.fx=mean)
   tn.ann <- annual_comparison(tn.ann.series,site.name,tn.ann.file,tn.ann.title,tn.ann.y)
   save(tn.ann,file=paste0(save.dir,'annual_average_tasmin_',site.name,'.RData'))

   tn.min.file <- paste0(plot.dir,'annual_min_tasmin_',site.name,'.png')   
   tn.min.title <- paste0(gsub('_',' ',site.name)," Annual Minimum Temperature")
   tn.min.y <- "Minimum Temperature (\u00B0C)"   
   tn.min.series <- annual_quantity(all.series$tasmin,ann.fx=min)
   tn.min <- annual_comparison(tn.min.series,site.name,tn.min.file,tn.min.title,tn.min.y)
   save(tn.min,file=paste0(save.dir,'annual_minimum_tasmin_',site.name,'.RData'))

   ##Seasonal Minimum 
   seasons <- c('Winter','Spring','Summer','Fall')
   tasmin.seasonal <- vector(mode='list',length=4)
   for (s in seq_along(seasons)) {
      season <- seasons[s]
      tn.seas.file <- paste0(plot.dir,'annual_',tolower(season),'_tasmin_',site.name,'.png')   
      tn.seas.title <- paste0(gsub('_',' ',site.name)," ",season," Minimum Temperature")
      tn.seas.y <- "Minimum Temperature (\u00B0C)"
      tn.seas.mean <- seasonal_quantity(all.series$tasmin,seas.fx=mean,seas=season)
      tasmin.seasonal[[s]] <- annual_comparison(tn.seas.mean,site.name,tn.seas.file,tn.seas.title,tn.seas.y)
   }
   save(tasmin.seasonal,file=paste0(save.dir,'seasonal_average_tasmin_',site.name,'.RData'))

   ##------------------------------------------------------------------------------------------------------
   ##Degree Days
   dd.list <- c('gdd','hdd','cdd','fdd')
   dd.fxns <- list(gdd=gdd,cdd=cdd,hdd=hdd,fdd=fdd)
   dd.types <- list(gdd='Growing Degree Days',cdd='Cooling Degree Days',
                    hdd='Heating Degree Days',fdd='Freezing Degree Days')

   for (d in 1:4) {
      dd.name <- dd.list[d]
      dd.type <- dd.types[[dd.name]]
      dd.fx <- dd.fxns[[dd.name]]
      dd.file <- paste0(plot.dir,'annual_average_',dd.name,'_',site.name,'.png')      
      dd.title <- paste0(gsub('_',' ',site.name)," Annual Average ",dd.type)
      dd.ann.y <- "Degree Days (DD)"   
      dd.ann.series <- dd_quantity(all.series$tas,dd.fx=dd.fx)
      dd.ann <- annual_comparison(dd.ann.series,site.name,dd.file,dd.title,dd.ann.y)
      save(dd.ann,file=paste0(save.dir,'annual_average_',dd.name,'_',site.name,'.RData'))
   }   

   ##PAS
   pas.ann.file <- paste0(plot.dir,'annual_total_pas_',site.name,'.png')   
   pas.ann.title <- paste0(gsub('_',' ',site.name)," Annual Total Precipitation-as-Snow")
   pas.ann.y <- "Precipitation-as-Snow (mm)"   
   pas.series <- compute_pas(all.series)
   pas.ann.series <- annual_quantity(pas.series,ann.fx=sum)
   pas.ann <- annual_comparison(pas.ann.series,site.name,pas.ann.file,pas.ann.title,pas.ann.y)
   save(pas.ann,file=paste0(save.dir,'annual_average_pas_',site.name,'.RData'))

   ##Seasonal Total Precipitation-as-Snow
   seasons <- c('Winter','Spring','Summer','Fall')
   pas.seasonal <- vector(mode='list',length=4)
   for (s in seq_along(seasons)) {
      season <- seasons[s]
      pas.seas.file <- paste0(plot.dir,'annual_',tolower(season),'_pas_',site.name,'.png')   
      pas.seas.title <- paste0(gsub('_',' ',site.name)," ",season," Total Precipitation-as-Snow")
      pas.seas.y <- "Precipitation-as-Snow (mm)"   
      pas.seas.total <- seasonal_quantity(pas.series,seas.fx=sum,seas=season)
      pas.seasonal[[s]] <- annual_comparison(pas.seas.total,site.name,pas.seas.file,pas.seas.title,pas.seas.y)
   }
   save(pas.seasonal,file=paste0(save.dir,'seasonal_total_pas_',site.name,'.RData'))

   ##Freeze-Thaw
   ftd.ann.file <- paste0(plot.dir,'annual_total_ftd_',site.name,'.png')   
   ftd.ann.title <- paste0(gsub('_',' ',site.name)," Annual Total Freeze-Thaw Days")
   ftd.ann.y <- "Freeze-Thaw Days (Days)"   
   ftd.series <- compute_ftd(all.series)
   ftd.ann.series <- annual_quantity(ftd.series,ann.fx=sum)
   ftd.ann <- annual_comparison(ftd.ann.series,site.name,ftd.ann.file,ftd.ann.title,ftd.ann.y)
   save(ftd.ann,file=paste0(save.dir,'annual_average_freeze_thaw_',site.name,'.RData'))

   ##Seasonal Total Freeze-Thaw
   seasons <- c('Winter','Spring','Summer','Fall')
   ftd.seasonal <- vector(mode='list',length=4)
   for (s in seq_along(seasons)) {
      season <- seasons[s]
      ftd.seas.file <- paste0(plot.dir,'annual_',tolower(season),'_ftd_',site.name,'.png')   
      ftd.seas.title <- paste0(gsub('_',' ',site.name)," ",season," Total Freeze-Thaw Days")
      ftd.seas.y <- "Freeze-Thaw (Days)"   
      ftd.seas.total <- seasonal_quantity(ftd.series,seas.fx=sum,seas=season)
      ftd.seasonal[[s]] <- annual_comparison(ftd.seas.total,site.name,ftd.seas.file,ftd.seas.title,ftd.seas.y)
   }
   save(ftd.seasonal,file=paste0(save.dir,'seasonal_total_freeze_thaw_',site.name,'.RData'))
   ##Climdex

}



   ##Spearman Correlation, Laggged Correlation

   ##

}
