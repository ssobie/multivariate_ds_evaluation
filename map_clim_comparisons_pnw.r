##Plot a projected map of Canada with the locations of all Canadian CWEC files
##identified on the map
##Have an option to represent sites as coloured circles so quantities
##can be added to the map.

source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/multivariate_ds_evaluation/pnw.canrcm4.plot.r',chdir=T)

library(graticule)
library(raster)
library(rgdal)
library(scales)

convert_to_proj_coords <- function(lon,lat,proj.crs="+init=epsg:3005") {

  d <- data.frame(lon=lon,lat=lat)
  coordinates(d) <- c('lon','lat')
  proj4string(d) <- CRS("+init=epsg:4326")
  d.albers <- spTransform(d,CRS(proj.crs))
  rv <- d.albers@coords
  return(rv)
}

##---------------------------------------------------------------------
##*********************************************************************
##Canadian Projection - Lambert Conformal Conic
pnw.crs <- "+init=epsg:3005"

##-------------------------------------------------------------------



make_comparison_plot <- function(plot.file,var.name,main.title,leg.title,
                                 ref.name,past.ref.brick,proj.ref.brick,
                                 sim.name,past.sim.brick,proj.sim.brick,
                                 past.int,proj.int,diff='anomaly',
                                 clim.breaks=NULL,diff.breaks=NULL,
                                 anom.breaks=NULL,anom.diff.breaks=NULL,bp=0) {

   ##----------------------
   aplin.mask <- brick("/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/anusplin.mask.nc")
   mask.proj <- projectRaster(aplin.mask,crs=CRS(pnw.crs))


   png(file=write.file,width=5,height=3.5,units='in',res=250,pointsize=6,bg='white')
   par(mfrow=c(3,3))
   par(mar=c(0,0,1,10))
   par(oma=c(4,4,4,1))

   past.ref.proj <- projectRaster(past.ref.brick,crs=CRS(pnw.crs)) * mask.proj

   past.sim.proj <- projectRaster(past.sim.brick,crs=CRS(pnw.crs)) * mask.proj  

   past.diff.brick <- past.sim.brick - past.ref.brick
   past.diff.proj <- projectRaster(past.diff.brick,crs=CRS(pnw.crs)) * mask.proj  

   proj.ref.proj <- projectRaster(proj.ref.brick,crs=CRS(pnw.crs)) * mask.proj

   proj.sim.proj <- projectRaster(proj.sim.brick,crs=CRS(pnw.crs)) * mask.proj  

   proj.diff.brick <- proj.sim.brick - proj.ref.brick
   proj.diff.proj <- projectRaster(proj.diff.brick,crs=CRS(pnw.crs)) * mask.proj  

   anom.ref.brick <- proj.ref.brick - past.ref.brick
   anom.ref.proj <- projectRaster(anom.ref.brick,crs=CRS(pnw.crs)) * mask.proj  

   anom.sim.brick <- proj.sim.brick - past.sim.brick
   anom.sim.proj <- projectRaster(anom.sim.brick,crs=CRS(pnw.crs)) * mask.proj  

   anom.diff.brick <- anom.sim.brick - anom.ref.brick
   anom.diff.proj <- projectRaster(anom.diff.brick,crs=CRS(pnw.crs)) * mask.proj  

   clim.range <- range(c(cellStats(past.ref.proj,range,na.rm=T),cellStats(past.sim.proj,range,na.rm=T)),
                       c(cellStats(proj.ref.proj,range,na.rm=T),cellStats(proj.sim.proj,range,na.rm=T)))
   diff.range <- range(c(cellStats(past.diff.proj,range,na.rm=T),cellStats(proj.diff.proj,range,na.rm=T)))
   anom.range <- range(c(cellStats(anom.ref.proj,range,na.rm=T),cellStats(anom.sim.proj,range,na.rm=T)))
   anom.diff.range <- cellStats(anom.diff.proj,range,na.rm=T)

   print(paste0('Map Range: ',clim.range))
   if (is.null(clim.breaks)) {
      clim.breaks <- pretty(clim.range) ##get.class.breaks(var.name,type='past',clim.range,manual.breaks='')
   }
   clim.class.breaks.labels <- get.class.break.labels(clim.breaks)

   print(paste0('Diff range: ',diff.range))
   if (is.null(diff.breaks)) {
      diff.breaks <- pretty(diff.range) 
   }
   diff.class.breaks.labels <- get.class.break.labels(diff.breaks)

   print(paste0('Anom diff range: ',anom.range))
   if (is.null(anom.breaks)) {
      anom.breaks <- pretty(anom.range)
   }
   anom.class.breaks.labels <- get.class.break.labels(anom.breaks)

   if (is.null(anom.diff.breaks)) {
      anom.diff.breaks <- pretty(anom.diff.range)
   }
   anom.diff.breaks.labels <- get.class.break.labels(anom.diff.breaks)


   bp <- bp
   colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=clim.range,
                                    my.bp=bp,class.breaks=clim.breaks,
                                    type='past')
   bp <- 0
   diff.ramp <- get.legend.colourbar(var.name=var.name,map.range=diff.range,
                                    my.bp=bp,class.breaks=diff.breaks,
                                    type='anomaly')

   anom.ramp <- get.legend.colourbar(var.name=var.name,map.range=anom.range,
                                    my.bp=bp,class.breaks=anom.breaks,
                                    type='anomaly')
   anom.diff.ramp <- get.legend.colourbar(var.name=var.name,map.range=anom.diff.range,
                                    my.bp=bp,class.breaks=anom.diff.breaks,
                                    type='anomaly')


   ##Past Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(ref.name," (",past.int,")"),
                      morph.proj=past.ref.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,y.axis=T,x.axis=FALSE,inset=c(-0.55,0))
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name," (",past.int,")"),
                      morph.proj=past.sim.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,x.axis=FALSE,inset=c(-0.55,0))

   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name,"-",ref.name," Diff."),
                      morph.proj=past.diff.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=FALSE)

   ##Proj Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(ref.name," (",proj.int,")"),
                      morph.proj=proj.ref.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,y.axis=T,x.axis=FALSE)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name," (",proj.int,")"),
                      morph.proj=proj.sim.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=FALSE)

   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name,"-",ref.name," Diff."),
                      morph.proj=proj.diff.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,x.axis=FALSE)

   ##Diff Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(ref.name," Anomaly"),
                      morph.proj=anom.ref.proj,
                      class.breaks=anom.breaks,colour.ramp=anom.ramp,
                      map.class.breaks.labels=anom.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,y.axis=T,x.axis=TRUE)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name," Anomaly"),
                      morph.proj=anom.sim.proj,
                      class.breaks=anom.breaks,colour.ramp=anom.ramp,
                      map.class.breaks.labels=anom.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,x.axis=TRUE)

   make_pnw_multiplot(var.name,col.name,plot.title="Anomaly Differences",
                      morph.proj=anom.diff.proj,
                      class.breaks=anom.diff.breaks,colour.ramp=anom.diff.ramp,
                      map.class.breaks.labels=anom.diff.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=TRUE)
                                 

   mtext('Longitude (\u00B0E)',side=1,outer=T,line=2.8,col='gray22')
   mtext('Latitude (\u00B0N)',side=2,outer=T,line=2.8,col='gray22')
   mtext(main.title,side=3,outer=T,line=0.75,col='gray22')

   dev.off()



}


