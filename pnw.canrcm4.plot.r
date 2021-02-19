##Plot a projected map of Canada with the locations of all Canadian CWEC files
##identified on the map
##Have an option to represent sites as coloured circles so quantities
##can be added to the map.

source('/storage/data/projects/rci/assessments/code/resource.region.map.support.r',chdir=T)

library(graticule)
library(raster)
library(rgdal)
library(scales)

add_graticules <- function(lons,lats,crs) {

  xl <-  range(lons)
  yl <- range(lats)
  grat <- graticule(lons, lats, proj = CRS(crs),xlim=c(-180,0),ylim=c(30,89))
  rv <- grat
  return(rv)
}

##X-Axis Ticks
get.proj.xaxis <- function(lons,crs,plot.window.ylim) {

  y <- seq(0,80,0.1)
  xm <- sapply(lons,rep,length(y))
  S <- apply(xm,2,function(x) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'y']-plot.window.ylim[1]))})
  xticks <- mapply(FUN=function(s,i){s@coords[,'x'][i]},S2,indices)
  return(xticks)
}


 ##Y-Axis Ticks
get.proj.yaxis <- function(lats,crs,plot.window.xlim) {

  x <- seq(-180,-80,0.1)
  ym <- sapply(lats,rep,length(x))
  S <- apply(ym,2,function(y) {SpatialPoints(cbind(x,y), proj4string = CRS("+proj=longlat +datum=WGS84"))})
  S2<- lapply(S,spTransform, crs)
  indices <- lapply(S2,function(x){which.min(abs(x@coords[,'x']-plot.window.xlim[1]))})
  yticks <- mapply(FUN=function(s,i){s@coords[,'y'][i]},S2,indices)
  return(yticks)
}

convert_to_proj_coords <- function(lon,lat,proj.crs="+init=epsg:3005") {

  d <- data.frame(lon=lon,lat=lat)
  coordinates(d) <- c('lon','lat')
  proj4string(d) <- CRS("+init=epsg:4326")
  d.albers <- spTransform(d,CRS(proj.crs))
  rv <- d.albers@coords
  return(rv)
}

##-------------------------------------------------------------------

make_pnw_plot <- function(morph.name,var.name,plot.dir,plot.file,plot.title,
                          morph.proj=NULL,epw.proj.coords=NULL,
                          class.breaks=NULL,colour.ramp=NULL,
                          map.class.breaks.labels=NULL,leg.title='') {
 
   pnw.crs <- "+init=epsg:3005"

   plot.window.xlim <- c(200000,2100000)
   plot.window.ylim <- c(300000,2000000)

   lons <- seq(-180,0,by=5)
   lats <- seq(30,90,by=2.5)
   grats <- add_graticules(lons,lats,pnw.crs)
  
   xtks <- get.proj.xaxis(lons,pnw.crs,plot.window.ylim)
   ytks <- get.proj.yaxis(lats,pnw.crs,plot.window.xlim)

   shape.dir <- '/storage/data/gis/basedata/base_layers'

   provinces.region <- 'canada_provinces'
   can.shp <- readOGR(shape.dir, provinces.region, stringsAsFactors=F, verbose=F)

   na.region <- 'north_america_state_provincial_boundaries'
   na.shp <- readOGR(shape.dir, na.region, stringsAsFactors=F, verbose=F)
   na.proj.shp <- spTransform(na.shp,CRS(pnw.crs))
   ocean.shp <- readOGR(shape.dir, 'ocean_sym_diff_continents', stringsAsFactors=F, verbose=F)
   us.shp <- readOGR(shape.dir, 'united_states', stringsAsFactors=F, verbose=F)
   lakes.shp <- readOGR(shape.dir, 'lakes', stringsAsFactors=F, verbose=F)
   rivers.shp <- readOGR(shape.dir, 'rivers', stringsAsFactors=F, verbose=F)

   ocean.transformed <- spTransform(ocean.shp, CRS(pnw.crs))
   us.transformed <- spTransform(us.shp, CRS(pnw.crs))
   can.transformed <- spTransform(can.shp, CRS(pnw.crs))
   lakes.transformed <- spTransform(lakes.shp, CRS(pnw.crs))
   rivers.transformed <- spTransform(rivers.shp, CRS(pnw.crs))

   ###plot.dir <- '/storage/data/projects/rci/weather_files/plots/'
   write.file <- paste0(plot.dir,plot.file)
  
   cx <- 1.5

   png(file=write.file,width=5,height=5,units='in',res=600,pointsize=6,bg='white')
   par(mar=c(4,4.3,4.1,2.1))
   plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',col.lab='gray22',col.main='gray22', # 'gray94',
   xlab='Longitude (\u00B0E)',ylab='Latitude (\u00B0N)',main=plot.title,
   cex.axis=cx,cex.lab=cx,cex.main=1.5,mgp=c(2.5,1.75,0),axes=F)
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')

   if (!is.null(morph.proj)) {  
      image(morph.proj,col=alpha(colour.ramp,0.9),breaks=class.breaks,axes=FALSE,
            xlim=plot.window.xlim,ylim=plot.window.ylim,bg='white',
            main='',xlab='',ylab='',add=T)
   }

   lw <- 0.3
##   plot(us.transformed,col='lightgray',add=TRUE,lwd=lw)
##   plot(lakes.transformed,col='aliceblue',border='aliceblue',add=TRUE,lwd=lw)
   plot(can.transformed,add=TRUE,lwd=lw)

##   plot(ocean.transformed,col='aliceblue',add=TRUE,lwd=lw)
   plot(grats,add=TRUE,lty=5,col='gray',lwd=lw)

   if (!is.null(epw.proj.coords)) {
      points(epw.proj.coords$lon,epw.proj.coords$lat,
             pch=21,cex=1.5,col='gray5',lwd=0.5,
             bg=epw.proj.coords$col)
##      points(epw.proj.coords,pch=18,cex=0.5,col='green')
   }
   
   axis(1,at=xtks,label=rep('',length(lons)),cex.axis=cx-1,col='gray22')
   axis(1,at=xtks,label=lons,cex.axis=cx,col='gray22',col.axis='gray22')
   axis(2,at=ytks,label=rep('',length(lats)),cex.axis=cx-1,col='gray22')
   axis(2,at=ytks,label=lats,cex.axis=cx,col='gray22',col.axis='gray22')

   legend('topright',col = "gray22", legend=rev(map.class.breaks.labels), 
           pch=22, pt.bg = rev(colour.ramp),bg=alpha('white',0.95),box.col='gray22',text.col='gray22',
           pt.cex=3.0, y.intersp=0.8, title.adj=0.2, title=leg.title, xjust=0, cex=1.25)
   box(which='plot',lwd=2,col='gray22')

   dev.off()
}


##---------------------------------------------------------------------

make_pnw_multiplot <- function(morph.name,var.name,plot.title,
                          morph.proj=NULL,epw.proj.coords=NULL,
                          class.breaks=NULL,colour.ramp=NULL,
                          map.class.breaks.labels=NULL,leg.title='',
                          x.axis=FALSE,y.axis=FALSE,leg.add=FALSE) {
 
   pnw.crs <- "+init=epsg:3005"

   plot.window.xlim <- c(200000,2100000)
   plot.window.ylim <- c(300000,2000000)

   lons <- seq(-180,0,by=5)
   lats <- seq(30,90,by=2.5)
   grats <- add_graticules(lons,lats,pnw.crs)
  
   xtks <- get.proj.xaxis(lons,pnw.crs,plot.window.ylim)
   ytks <- get.proj.yaxis(lats,pnw.crs,plot.window.xlim)

   shape.dir <- '/storage/data/gis/basedata/base_layers'

   provinces.region <- 'canada_provinces'
   can.shp <- readOGR(shape.dir, provinces.region, stringsAsFactors=F, verbose=F)

   na.region <- 'north_america_state_provincial_boundaries'
   na.shp <- readOGR(shape.dir, na.region, stringsAsFactors=F, verbose=F)
   na.proj.shp <- spTransform(na.shp,CRS(pnw.crs))
   ocean.shp <- readOGR(shape.dir, 'ocean_sym_diff_continents', stringsAsFactors=F, verbose=F)
   us.shp <- readOGR(shape.dir, 'united_states', stringsAsFactors=F, verbose=F)
   lakes.shp <- readOGR(shape.dir, 'lakes', stringsAsFactors=F, verbose=F)
   rivers.shp <- readOGR(shape.dir, 'rivers', stringsAsFactors=F, verbose=F)

   ocean.transformed <- spTransform(ocean.shp, CRS(pnw.crs))
   us.transformed <- spTransform(us.shp, CRS(pnw.crs))
   can.transformed <- spTransform(can.shp, CRS(pnw.crs))
   lakes.transformed <- spTransform(lakes.shp, CRS(pnw.crs))
   rivers.transformed <- spTransform(rivers.shp, CRS(pnw.crs))
 
   cx <- 1.5

   plot(c(),xlim=plot.window.xlim,ylim=plot.window.ylim,xaxs='i',yaxs='i',
   bg='white',col.lab='gray22',col.main='gray22', # 'gray94',
   xlab='',ylab='',main='',
   cex.axis=cx,cex.lab=cx,cex.main=1.5,mgp=c(2.5,1.75,0),axes=F)
   rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col='lightgray')

   if (!is.null(morph.proj)) {  
      image(morph.proj,col=alpha(colour.ramp,0.9),breaks=class.breaks,axes=FALSE,
            xlim=plot.window.xlim,ylim=plot.window.ylim,bg='white',
            main='',xlab='',ylab='',add=T)
   }

   lw <- 0.3
##   plot(us.transformed,col='lightgray',add=TRUE,lwd=lw)
##   plot(lakes.transformed,col='aliceblue',border='aliceblue',add=TRUE,lwd=lw)
   plot(can.transformed,add=TRUE,lwd=lw)

##   plot(ocean.transformed,col='aliceblue',add=TRUE,lwd=lw)
   plot(grats,add=TRUE,lty=5,col='gray',lwd=lw)


   rect(0.3*par("usr")[2], 0.92*par("usr")[4], 0.8*par("usr")[2], par("usr")[4], col='white',border='gray22')
   text(0.55*par("usr")[2], 0.96*par("usr")[4],plot.title,col='gray22')

   if (!is.null(epw.proj.coords)) {
      points(epw.proj.coords$lon,epw.proj.coords$lat,
             pch=21,cex=1.5,col='gray5',lwd=0.5,
             bg=epw.proj.coords$col)
##      points(epw.proj.coords,pch=18,cex=0.5,col='green')
   }

   if (x.axis) {   
      axis(1,at=xtks,label=rep('',length(lons)),cex.axis=cx-1,col='gray22')
      axis(1,at=xtks,label=lons,cex.axis=cx,col='gray22',col.axis='gray22')
   }
   if (y.axis) {
      axis(2,at=ytks,label=rep('',length(lats)),cex.axis=cx-1,col='gray22')
      axis(2,at=ytks,label=lats,cex.axis=cx,col='gray22',col.axis='gray22')
   }
   if (leg.add) {
      par(xpd=NA)
      legend('topright', col = "gray22", legend=rev(map.class.breaks.labels), inset=c(-0.55,0),
           pch=22, pt.bg = rev(colour.ramp),bg=alpha('white',0.95),box.col='gray22',text.col='gray22',
           pt.cex=3.0, y.intersp=0.8, title.adj=0.2, title=leg.title, xjust=0, cex=1.25)
      par(xpd=FALSE)
   }
   box(which='plot',lwd=1,col='gray22')
}



##---------------------------------------------------------------------

##*********************************************************************
