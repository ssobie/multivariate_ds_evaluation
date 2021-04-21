##Compare series at one cell for all the different datasets

library(scales)
library(ncdf4)
library(PCICt)
library(xtable)

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----------------------------------------------------
##Comparisons to consider
##Monthly total precipitation with monthly average temperatures

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

monthly_quantity <- function(input.series,mon.fx) {

   dates <- input.series$rcm$time
   mon.fac <- as.factor(format(dates,'%Y-%m'))
   x.vals <- levels(mon.fac)

   rcm.mon <- tapply(input.series$rcm$data,mon.fac,mon.fx)
   bcci.mon <- tapply(input.series$bcci$data,mon.fac,mon.fx)
   bccaqv2.mon <- tapply(input.series$bccaqv2$data,mon.fac,mon.fx)
   mbcn.mon <- tapply(input.series$mbcn$data,mon.fac,mon.fx)
   rv <- list(rcm=rcm.mon,bcci=bcci.mon,bccaqv2=bccaqv2.mon,mbcn=mbcn.mon,x=x.vals)

   return(rv)

}

compound_months <- function(pr.series,temp.series,pr.pctl,temp.pctl,comp_fxn) {

   intervals <- list(c(1951,1980),
                     c(1981,2010),
                     c(2011,2040),
                     c(2041,2070),
                     c(2071,2100))

   types <- c('rcm','bcci','bccaqv2','mbcn')
   compound.events <- matrix(NA,nrow=length(types),ncol=length(intervals))

   for (j in 1:5) {
      bnds <- intervals[[j]]
      yst <- head(grep(bnds[1],pr.series$x),1)
      yen <- tail(grep(bnds[2],pr.series$x),1)

      for (k in seq_along(types)) {
         type <- types[k]
         pr.sub <- pr.series[[type]][yst:yen]
         temp.sub <- temp.series[[type]][yst:yen]
         pr.qval <- quantile(pr.sub,pr.pctl/100,names=FALSE)
         temp.qval <- quantile(temp.sub,temp.pctl/100,names=FALSE)
         test <- comp_fxn(pr.sub,pr.qval,temp.sub,temp.qval)
         compound.events[k,j] <- comp_fxn(pr.sub,pr.qval,temp.sub,temp.qval)

      }
   }

   return(round(compound.events*100,2))
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
      load.file <- paste0(var.name,"_day_MBCn_iterated_trace_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")
      load(paste0(load.dir,load.file))
      type.series$mbcn <- var.series
      all.series[[var.name]] <- type.series
   }

   pr.mon.series <- monthly_quantity(all.series$pr,mon.fx=sum)
   tx.mon.series <- monthly_quantity(all.series$tasmax,mon.fx=mean)
   tn.mon.series <- monthly_quantity(all.series$tasmin,mon.fx=mean)
   ts.mon.series <- monthly_quantity(all.series$tas,mon.fx=mean)

   ##For each interval compare the months in each season
   ##Calculate the 25th and 75th percentiles of precip
   ##Calculate the 10th and 90th percentiles of temperatures
   ##Count the coincident events of the 4 possible pairs

   comp_fxn <- function(pr,pr.q,tas,tas.q) { sum((pr > pr.q) & (tas < tas.q))/length(pr)}
   wet.cold <- compound_months(pr.mon.series,tn.mon.series,pr.pctl=75,temp.pctl=25,comp_fxn)      
   
   comp_fxn <- function(pr,pr.q,tas,tas.q) { sum((pr < pr.q) & (tas < tas.q))/length(pr)}
   dry.cold <- compound_months(pr.mon.series,tn.mon.series,pr.pctl=25,temp.pctl=25,comp_fxn)      

   comp_fxn <- function(pr,pr.q,tas,tas.q) { sum((pr < pr.q) & (tas > tas.q))/length(pr)}
   dry.warm <- compound_months(pr.mon.series,tx.mon.series,pr.pctl=25,temp.pctl=75,comp_fxn)      

   comp_fxn <- function(pr,pr.q,tas,tas.q) { sum((pr > pr.q) & (tas > tas.q))/length(pr)}
   wet.warm <- compound_months(pr.mon.series,tx.mon.series,pr.pctl=75,temp.pctl=75,comp_fxn)      


   compound.months <- list(dry_cold=dry.cold,dry_warm=dry.warm,
                           wet_warm=wet.warm,wet_cold=wet.cold)

   save(compound.months,file=paste0(save.dir,"compound_months_25_75_",site.name,".trace.RData"))                            
   
}
