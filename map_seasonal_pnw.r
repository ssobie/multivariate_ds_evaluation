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

mbcn.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/MBCn_Derived/seasonal/"
can.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/CanRCM4_Derived/seasonal/"

main.title <- 'CanRCM4 MBCn Comparison (Precipitation)'

proj.bnds <- c(1981,2010)
proj.int <- paste0(proj.bnds,collapse='-')


plot.file <- paste0('seasonal.pr.comparison.',proj.int,'.png')
write.file <- paste0(plot.dir,plot.file)

##----------------------

png(file=write.file,width=5,height=5,units='in',res=600,pointsize=6,bg='white')
par(mfrow=c(5,3))
par(mar=c(0,0,0,9))
par(oma=c(4,4,4,1))

seasons <- c('winter','spring','summer','fall','annual')

for (s in seq_along(seasons)) {

   seas.grep <- switch(season,
                    winter='(12|01|02)', ##Winter
                    spring='(03|04|05)', ##Spring
                    summer='(06|07|08)', ##Summer
                    fall='(09|10|11)', ##Fall
                    annual='*') ##Annual
   

   season <- seasons[s]
   print(season)

   mbcn.file <- paste0(var.name,'_MBCn_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_1950-2100.nc')
   mbcn.brick <- brick(paste0(clim.dir,mbcn.file))
   mbcn.dates <- mbcn.brick@z$Date
   mbcn.date.ix <- which(mbcn.dates >= as.Date(paste0(proj.bnds[1],'-01-01')) & 
                    mbcn.dates <= as.Date(paste0(proj.bnds[2],'-12-31')))
   mbcn.date.sub <- subset(mbcn.brick,mbcn.date.ix)
browser()
   mbcn.seas.ix <- grep(seas.grep,format(mbcn.dates[mbcn.date.ix],'%m'))


                                  
   canrcm4.file <- paste0(var.name,'_seasonal_total_CanRCM4_historical+rcp85_r1i1p1_1950-2100.nc')


   mbcn.brick <- subset(),s)
   can4.brick <- subset(brick(paste0(clim.dir,can4.file)),s)
   diff.brick <- mbcn.brick - can4.brick

   mbcn.proj <- projectRaster(mbcn.brick,crs=CRS(pnw.crs))
   can4.proj <- projectRaster(can4.brick,crs=CRS(pnw.crs))
   diff.proj <- projectRaster(diff.brick,crs=CRS(pnw.crs))

   if (obs=='CanRCM4' | obs=='MBCn') {
      mbcn.proj <- mbcn.proj * mask.proj
      can4.proj <- can4.proj * mask.proj
      diff.proj <- diff.proj * mask.proj
   }

   map.range <- range(c(cellStats(mbcn.proj,range,na.rm=T),cellStats(can4.proj,range,na.rm=T)))
   print(paste0('Map Range: ',map.range))
   class.breaks <- pretty(map.range) ##get.class.breaks(var.name,type='past',map.range,manual.breaks='')

   map.class.breaks.labels <- get.class.break.labels(class.breaks)

   diff.range <- cellStats(diff.proj,range,na.rm=T)
   print(paste0('Diff range: ',diff.range))
   diff.breaks <- pretty(diff.range)
   diff.class.breaks.labels <- get.class.break.labels(diff.breaks)

   bp <- 0
   colour.ramp <- get.legend.colourbar(var.name=var.name,map.range=map.range,
                                    my.bp=bp,class.breaks=class.breaks,
                                    type='past')
   diff.ramp <- get.legend.colourbar(var.name=var.name,map.range=diff.range,
                                    my.bp=bp,class.breaks=diff.breaks,
                                    type='anomaly')

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


   make_pnw_multiplot(var.name,col.name,plot.title=paste0(tools::toTitleCase(season)," (",proj.int,")"),
                      morph.proj=past.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
                      leg.title=leg.units,leg.add=leg.add,y.axis=T,x.axis=x.axis)
   make_pnw_multiplot(var.name,col.name,plot.title=paste0(tools::toTitleCase(season)," (",proj.int,")"),
                      morph.proj=proj.proj,
                      class.breaks=class.breaks,colour.ramp=colour.ramp,
                      map.class.breaks.labels=map.class.breaks.labels,
                      leg.title=leg.units,leg.add=leg.add,x.axis=x.axis)

   make_pnw_multiplot(var.name,col.name,plot.title=paste0(tools::toTitleCase(season)," Diff."),
                      morph.proj=diff.proj,
                      class.breaks=diff.breaks,colour.ramp=diff.ramp,
                      map.class.breaks.labels=diff.class.breaks.labels,
                      leg.title=leg.units,leg.add=leg.add,x.axis=x.axis)

   mtext('Longitude (\u00B0C)',side=1,outer=T,line=2.8,col='gray22')
   mtext('Latitude (\u00B0C)',side=2,outer=T,line=2.8,col='gray22')
   mtext(main.title,side=3,outer=T,line=2,col='gray22')

}

dev.off()




