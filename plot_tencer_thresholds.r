##PLot the 365 threshold values for each interval and model type

library(scales)
##--------------------------------------------------------------

sites <- c('Fort_Nelson') ##,'Skagway','Prince_George','Prince_Rupert','Golden','Ucluelet')

var.names <- c('pr') ##,'tasmax','tasmin')
intervals <- c('1951-1980','1981-2010','2011-2040','2041-2070','2071-2100')

plot.col <- 'gray22'
line.cols <- c('green','orange','red')
cx <- 0.7
lw <- 2
x <- 1:365

for (site in sites) {
   print(site)
   load.file <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/RData/",
                      site,"/tencer_thresholds_",site,".trace.RData")
   load(load.file)                     
   plot.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/plots/",
                      site,"/")

   for (var.name in var.names) {
     
      plot.file <- paste0(plot.dir,var.name,"_tencer_thresholds_",site,".trace.png")
      tvar <- switch(var.name,pr='pr',tasmax='tx',tasmin='tn')
      y.label <- switch(var.name,pr='Precipitation (mm)',
                                 tasmax='Max. Temperature ()',
                                 tasmin='Min. Temperature ()')
      var.thresh <- thresholds[[tvar]]
   
      ##Make Plot
      png(plot.file,width=4,height=5,units='in',res=300)   
      layout(rbind(c(1,1,2,2),
                c(3,3,4,4),
                c(0,5,5,0)))
      par(mar=c(2,2,1,1))
      par(oma=c(2,2,2,0))
                
      for (i in 1:5) { ##Intervals
         plot.vals <- c(var.thresh$rcm[i,],var.thresh$bccaqv2[i,],var.thresh$mbcn[i,])
         plot.title <- intervals[i]
         y.lim <- range(pretty(plot.vals))
      
         plot(c(),xlim=c(0,366),ylim=y.lim,
             xlab='Day',ylab=y.label,main='',
             cex=cx,cex.axis=cx,cex.lab=cx,axes=FALSE,
             col=plot.col,col.axis=plot.col,col.lab=plot.col)

         lines(x,var.thresh$rcm[i,],col=line.cols[1],lwd=lw)
         lines(x,var.thresh$bccaqv2[i,],col=line.cols[2],lwd=lw)
         lines(x,var.thresh$mbcn[i,],col=line.cols[3],lwd=lw)
         axis(1,at=pretty(1:365),label=pretty(1:365),col.axis=plot.col)
         axis(2,at=pretty(plot.vals),label=pretty(plot.vals),col.axis=plot.col)
         box(which='plot',col=plot.col)

         bccaqv2.rmse <- sqrt(mean((var.thresh$bccaqv2[i,] - var.thresh$rcm[i,])^2))
         mbcn.rmse <- sqrt(mean((var.thresh$mbcn[i,] - var.thresh$rcm[i,])^2))

         text(x=0.2*par('usr')[2],y=0.075*diff(par('usr')[3:4])+par('usr')[3],
               label=paste0('BCCAQv2: ',round(bccaqv2.rmse,2)),col=plot.col,cex=0.8)
         text(x=0.8*par('usr')[2],y=0.075*diff(par('usr')[3:4])+par('usr')[3],
              label=paste0('MBCn: ',round(mbcn.rmse,2)),col=plot.col,cex=0.8)


         text(x=0.13*par('usr')[2],y=0.95*par('usr')[4],label=plot.title,col=plot.col,cex=0.8)
         if (i==5) {
            coord <- par('usr')
            par(xpd=NA)
            legend(x=1.05*coord[2],y=coord[4], col = "gray22", legend=c('CanRCM4','BCCAQv2','MBCn'),
               pch=22, pt.bg = line.cols,bg=alpha('white',0.95),box.col='gray22',text.col='gray22',
               pt.cex=2.0, y.intersp=0.8, xjust=0, cex=1)
            par(xpd=FALSE)
         }
      } ##Intervals
      mtext('Day of the Year',side=1,outer=TRUE,col=plot.col,cex=cx)
      mtext(y.label,side=2,outer=TRUE,col=plot.col,cex=cx)
      mtext(gsub('_',' ',site),side=3,outer=TRUE,col=plot.col,cex=cx)


##      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
##      plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
##      legend("bottomleft", leg=c('CanRCM4','BCCAQv2','MBCn'), xpd = TRUE, inset = c(0,
##            0), bty = "n", pch = c(15, 15), col = line.cols, cex = cx, pt.cex=2,
##            box.col=plot.col,title.col=plot.col,text.col=plot.col)

##      par(xpd=TRUE)
##      legend('bottomleft',leg=c('CanRCM4','BCCAQv2','MBCn'),pch=15,
##             col=line.cols,pt.cex=2,box.col=plot.col,title.col=plot.col,text.col=plot.col)
##      par(xpd=FALSE)       

      dev.off()
browser()
   } ##Varnames
browser()   
}## Sites
