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

clim.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/Derived/spearman/"
var.pairs <- c('pr_tasmax','pr_tasmin','pr_tas','tasmax_tasmin')

obs <- 'ANUSPLIN'
var.pair <- 'pr_tasmax'
main.title <- paste0(obs,' Spearman Correlation (Precipitation and Max. Temperature)')


past.int <- '1951-1980'
proj.int <- '1981-2010'


plot.file <- paste0(gsub('_','.',var.pair),'.spearman.seasonal.',obs,'.png')
write.file <- paste0(plot.dir,plot.file)

##----------------------

png(file=write.file,width=5,height=5,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(5,3))
par(mar=c(0,0,0,9))
par(oma=c(4,4,4,1))

seasons <- c('winter','spring','summer','fall','annual')

for (s in seq_along(seasons)) {

   season <- seasons[s]
   print(season)

   past.clim.file <- paste0(season,'_spearman_corr_',var.pair,"_BC_",obs,"_",past.int,".nc")
   proj.clim.file <- paste0(season,'_spearman_corr_',var.pair,"_BC_",obs,"_",proj.int,".nc")


   past.brick <- subset(brick(paste0(clim.dir,past.clim.file)),1)
   proj.brick <- subset(brick(paste0(clim.dir,proj.clim.file)),1)
   diff.brick <- proj.brick - past.brick

   past.proj <- projectRaster(past.brick,crs=CRS(pnw.crs))
   proj.proj <- projectRaster(proj.brick,crs=CRS(pnw.crs))
   diff.proj <- projectRaster(diff.brick,crs=CRS(pnw.crs))

   if (obs=='CanRCM4') {
      past.proj <- past.proj * mask.proj
      proj.proj <- proj.proj * mask.proj
      diff.proj <- diff.proj * mask.proj
   }

   map.range <- range(c(cellStats(past.proj,range,na.rm=T),cellStats(proj.proj,range,na.rm=T)))
   print(paste0('Map Range: ',map.range))
   class.breaks <- round(seq(-0.7,0.8,0.1),1) ## pretty(map.range) ##get.class.breaks(var.name,type='past',map.range,manual.breaks='')
   ###class.breaks <- round(seq(0,1,0.1),1)
   map.class.breaks.labels <- get.class.break.labels(class.breaks)

   diff.range <- cellStats(diff.proj,range,na.rm=T)
   print(paste0('Diff range: ',diff.range))
   diff.breaks <- round(seq(-0.2,0.25,0.05),2) ##pretty(diff.range)
   diff.class.breaks.labels <- get.class.break.labels(diff.breaks)

   var.name <- 'tasmax'
   col.var <- 'tasmax'
   ##class.breaks <- round(c(seq(-0.6,-0.1,by=0.1),seq(0.1,0.7,by=0.1)),1)

   bp <- 0
   colour.ramp <- get.legend.colourbar(var.name=col.var,map.range=map.range,
                                    my.bp=bp,class.breaks=class.breaks,
                                    type='past')
   diff.ramp <- get.legend.colourbar(var.name=col.var,map.range=diff.range,
                                    my.bp=bp,class.breaks=diff.breaks,
                                    type='past')


   if (season=='annual') {
      x.axis=TRUE
   } else {
      x.axis=FALSE
   }
   if (season=='winter') {
      leg.add=TRUE
   } else {
      leg.add=FALSE
   }


   make_pnw_multiplot(var.name,col.name,plot.title=paste0(tools::toTitleCase(season)," (1951-1980)"),
                      morph.proj=past.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
                      leg.title='Sp. Corr.',leg.add=leg.add,y.axis=T,x.axis=x.axis)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(tools::toTitleCase(season)," (1981-2010)"),
                      morph.proj=proj.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
                      leg.title='Sp. Corr.',leg.add=leg.add,x.axis=x.axis)

   make_pnw_multiplot(var.name,col.name,plot.title=paste0(tools::toTitleCase(season)," Diff."),
                      morph.proj=diff.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title='Sp. Corr.',leg.add=leg.add,x.axis=x.axis)

   mtext('Longitude (\u00B0C)',side=1,outer=T,line=2.8,col='gray22')
   mtext('Latitude (\u00B0C)',side=2,outer=T,line=2.8,col='gray22')
   mtext(main.title,side=3,outer=T,line=2,col='gray22')

}

dev.off()




