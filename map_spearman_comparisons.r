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

plot.dir <- '/storage/data/projects/rci/data/stat.downscaling/MBCn/CanRCM4/'

aplin.mask <- brick("/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/anusplin.mask.nc")
mask.proj <- projectRaster(aplin.mask,crs=CRS(pnw.crs))



ref.name <- ref <- 'CanRCM4'
ref.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",ref,"_Derived/spearman/")

sim.name <- sim <- 'MBCn'
sim.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",sim,"_Derived/spearman/")

var.pairs <- c('pr_tasmax','pr_tasmin') 

for (var.pair in var.pairs) {



leg.title <- 'Sp. Corr.'
past.int <- '1951-1980'
proj.int <- '2041-2070'

if (var.pair=='pr_tasmax') {
  var.title <- paste0('(Precipitation and Max. Temperature, ',proj.int,')')
} 

if (var.pair=='pr_tasmin') {
  var.title <- paste0('(Precipitation and Min. Temperature, ',proj.int,')')
} 



##--------------------------------------------------------
seasons <- 'summer' ##c('winter','summer','annual') ##'spring','summer','fall','annual')

for (s in seq_along(seasons)) {

   season <- seasons[s]
   print(season)
   cap.seas <- tools::toTitleCase(season)     
   main.title <- paste0(cap.seas,' Spearman Correlation Comparison\n',var.title)

   plot.file <- paste0(season,'.',gsub('_','.',var.pair),'.spearman.seasonal.',ref,'.',sim,'.',proj.int,'.png')
   write.file <- paste0(plot.dir,plot.file)

   ##----------------------

   png(file=write.file,width=5,height=3.5,units='in',res=300,pointsize=6,bg='white')
   par(mfrow=c(3,3))
   par(mar=c(0,0,1,10))
   par(oma=c(4,4,4,1))

   past.ref.file <- paste0(season,'_spearman_corr_',var.pair,"_BC_",ref,"_",past.int,".nc")
   past.ref.brick <- subset(brick(paste0(ref.dir,past.ref.file)),1)
   past.ref.proj <- projectRaster(past.ref.brick,crs=CRS(pnw.crs)) * mask.proj

   past.sim.file <- paste0(season,'_spearman_corr_',var.pair,"_BC_",sim,"_",past.int,".nc")
   past.sim.brick <- subset(brick(paste0(sim.dir,past.sim.file)),1)
   past.sim.proj <- projectRaster(past.sim.brick,crs=CRS(pnw.crs)) * mask.proj  

   past.diff.brick <- past.sim.brick - past.ref.brick
   past.diff.proj <- projectRaster(past.diff.brick,crs=CRS(pnw.crs)) * mask.proj  

   proj.ref.file <- paste0(season,'_spearman_corr_',var.pair,"_BC_",ref,"_",proj.int,".nc")
   proj.ref.brick <- subset(brick(paste0(ref.dir,proj.ref.file)),1)
   proj.ref.proj <- projectRaster(proj.ref.brick,crs=CRS(pnw.crs)) * mask.proj

   proj.sim.file <- paste0(season,'_spearman_corr_',var.pair,"_BC_",sim,"_",proj.int,".nc")
   proj.sim.brick <- subset(brick(paste0(sim.dir,proj.sim.file)),1)
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
   class.breaks <- pretty(clim.range) ##get.class.breaks(var.name,type='past',clim.range,manual.breaks='')
   class.breaks <- round(seq(-0.7,0.8,0.1),1)
   map.class.breaks.labels <- get.class.break.labels(class.breaks)

   print(paste0('Diff range: ',diff.range))
   diff.breaks <- round(seq(-0.25,0.3,0.05),2) ##pretty(diff.range) ##
   diff.class.breaks.labels <- get.class.break.labels(diff.breaks)

   print(paste0('Anom diff range: ',anom.range))
   anom.breaks <- round(seq(-0.25,0.2,0.05),2) ## seq(pretty(anom.range) ##round(seq(-0.2,0.25,0.05),2) ##
   anom.class.breaks.labels <- get.class.break.labels(anom.breaks)

   anom.diff.breaks <- diff.breaks ##pretty(anom.diff.range) ##round(seq(-0.2,0.25,0.05),2) ##
   anom.diff.breaks.labels <- get.class.break.labels(anom.diff.breaks)



   var.name <- 'tasmax'
   col.var <- 'tasmax'
   ##class.breaks <- round(c(seq(-0.6,-0.1,by=0.1),seq(0.1,0.7,by=0.1)),1)

   bp <- 0
   colour.ramp <- get.legend.colourbar(var.name=col.var,map.range=clim.range,
                                    my.bp=bp,class.breaks=class.breaks,
                                    type='past')
   diff.ramp <- get.legend.colourbar(var.name=col.var,map.range=diff.range,
                                    my.bp=bp,class.breaks=diff.breaks,
                                    type='anomaly')
   anom.ramp <- get.legend.colourbar(var.name=col.var,map.range=anom.range,
                                    my.bp=bp,class.breaks=anom.breaks,
                                    type='anomaly')
   anom.diff.ramp <- get.legend.colourbar(var.name=col.var,map.range=anom.diff.range,
                                    my.bp=bp,class.breaks=anom.diff.breaks,
                                    type='anomaly')

 
   ##Past Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(ref.name," (",past.int,")"),
                      morph.proj=past.ref.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,y.axis=T,x.axis=FALSE)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name," (",past.int,")"),
                      morph.proj=past.sim.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
                      leg.title=leg.title,leg.add=TRUE,x.axis=FALSE)

   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name,"-",ref.name," Diff."),
                      morph.proj=past.diff.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,x.axis=FALSE)

   ##Proj Comparison
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(ref.name," (",proj.int,")"),
                      morph.proj=proj.ref.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
                      leg.title=leg.title,leg.add=FALSE,y.axis=T,x.axis=FALSE)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(sim.name," (",proj.int,")"),
                      morph.proj=proj.sim.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
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
}





