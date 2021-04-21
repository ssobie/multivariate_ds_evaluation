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

make_count_plot <- function(plot.file,var.name,main.title,leg.title,
                                 ref.name,sim.name,mbc.name,
                                 past.ref.brick,past.sim.brick,past.mbc.brick,
                                 interval,diff='anomaly',
                                 clim.breaks=NULL,diff.breaks=NULL,
                                 anom.breaks=NULL,anom.diff.breaks=NULL,bp=0,
                                 c.g=FALSE,c.l=FALSE,d.g=FALSE,d.l=FALSE) {

   ##----------------------
   aplin.mask <- brick("/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/anusplin.mask.nc")
   mask.proj <- projectRaster(aplin.mask,crs=CRS(pnw.crs))

   png(file=write.file,width=5,height=3.5,units='in',res=250,pointsize=6,bg='white')
   par(mfrow=c(3,3))
   par(mar=c(0,0,1,10))
   par(oma=c(4,4,4,1))

   past.ref.proj <- projectRaster(past.ref.brick,crs=CRS(pnw.crs)) * mask.proj
   past.sim.proj <- projectRaster(past.sim.brick,crs=CRS(pnw.crs)) * mask.proj  
   past.mbc.proj <- projectRaster(past.mbc.brick,crs=CRS(pnw.crs)) * mask.proj  

   sim.diff.brick <- past.sim.brick - past.ref.brick
   sim.diff.proj <- projectRaster(sim.diff.brick,crs=CRS(pnw.crs)) * mask.proj  

   mbc.diff.brick <- past.mbc.brick - past.ref.brick
   mbc.diff.proj <- projectRaster(mbc.diff.brick,crs=CRS(pnw.crs)) * mask.proj  

   mbc.sim.brick <- past.mbc.brick - past.sim.brick
   mbc.sim.proj <- projectRaster(mbc.sim.brick,crs=CRS(pnw.crs)) * mask.proj  

   clim.range <- range(c(cellStats(past.ref.proj,range,na.rm=T),cellStats(past.sim.proj,range,na.rm=T),
                       cellStats(past.mbc.proj,range,na.rm=T)))
   diff.range <- range(c(cellStats(sim.diff.proj,range,na.rm=T),cellStats(mbc.diff.proj,range,na.rm=T),
                       cellStats(mbc.sim.proj,range,na.rm=T)))

   print(paste0('Map Range: ',round(clim.range,3)))
   if (is.null(clim.breaks)) {
      clim.breaks <- pretty(clim.range) ##get.class.breaks(var.name,type='past',clim.range,manual.breaks='')
   }
   clim.class.breaks.labels <- get.class.break.labels(clim.breaks,greater.sign=c.g,lesser.sign=c.l)

   print(paste0('Diff range: ',round(diff.range,3)))
   if (is.null(diff.breaks)) {
      diff.breaks <- pretty(diff.range) 
   }
   diff.class.breaks.labels <- get.class.break.labels(diff.breaks,greater.sign=d.g,lesser.sign=d.l)

   bp <- 0
   colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=clim.range,
                                    my.bp=bp,class.breaks=clim.breaks,
                                    type='past')
   bp <- 0
   diff.ramp <- get.legend.colourbar(var.name=var.name,map.range=diff.range,
                                    my.bp=bp,class.breaks=diff.breaks,
                                    type='anomaly')

   ##BCCAQv2 Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(ref.name," (",interval,")"),
                      morph.proj=past.ref.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,y.axis=T,x.axis=FALSE)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name," (",interval,")"),
                      morph.proj=past.sim.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=FALSE)

   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name,"-",ref.name," Diff."),
                      morph.proj=sim.diff.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=FALSE)

   ##MBC Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(ref.name," (",interval,")"),
                      morph.proj=past.ref.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,y.axis=T,x.axis=FALSE)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(mbc.name," (",interval,")"),
                      morph.proj=past.mbc.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,x.axis=FALSE)

   make_pnw_multiplot(var.name,col.name,plot.title=paste0(mbc.name,"-",ref.name," Diff."),
                      morph.proj=mbc.diff.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,x.axis=FALSE)

   ##Diff Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name," (",interval,")"),
                      morph.proj=past.sim.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,y.axis=T,x.axis=TRUE)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(mbc.name," (",interval,")"),
                      morph.proj=past.mbc.proj,
                      class.breaks=clim.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=clim.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=TRUE)

   make_pnw_multiplot(var.name,col.name,plot.title="MBCn-BCCAQv2 Diff.",
                      morph.proj=mbc.sim.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=TRUE)
                                 

   mtext('Longitude (\u00B0E)',side=1,outer=T,line=2.8,col='gray22')
   mtext('Latitude (\u00B0N)',side=2,outer=T,line=2.8,col='gray22')
   mtext(main.title,side=3,outer=T,line=0.75,col='gray22')

   dev.off()



}


