##Plot a projected map of Canada with the locations of all Canadian CWEC files
##identified on the map
##Have an option to represent sites as coloured circles so quantities
##can be added to the map.

source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)
source('/storage/home/ssobie/code/repos/MBC/R/figures/pnw.canrcm4.plot.r',chdir=T)

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

clim.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/Derived/"
clim.file <- "tasmax_annual_climatology_CanRCM4_1951-1980.nc"

clim.brick <- subset(brick(paste0(clim.dir,clim.file)),1)
clim.proj <- projectRaster(clim.brick,crs=CRS(pnw.crs))

map.range <- cellStats(clim.proj,range,na.rm=T)
class.breaks <- pretty(map.range) ##get.class.breaks(var.name,type='past',map.range,manual.breaks='')
map.class.breaks.labels <- get.class.break.labels(class.breaks)

var.name <- 'tasmax'
col.var <- 'tasmax'
##class.breaks <- round(c(seq(-0.6,-0.1,by=0.1),seq(0.1,0.7,by=0.1)),1)

bp <- 0
colour.ramp <- get.legend.colourbar(var.name=col.var,map.range=map.range,
                                    my.bp=bp,class.breaks=class.breaks,
                                    type='past')



##-------------------------------------------------------------------

plot.title <- 'TASMAX ANUSPLIN TEST'
plot.dir <- '/storage/data/projects/rci/data/stat.downscaling/MBCn/CanRCM4/'
plot.file <- 'tx.can.test.png'

make_pnw_plot(var.name,col.name,plot.dir,plot.file,plot.title,
                 morph.proj=clim.proj,
                 class.breaks=class.breaks,colour.ramp=colour.ramp,
                 map.class.breaks.labels=map.class.breaks.labels,
                 leg.title='degC')







