##Test using Bivariate Conditional Bias Correction with one cell of data

library(scales)
library(ncdf4)
library(PCICt)
source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----------------------------------------------------
##Extract a time series of BCCI 


extract_series <- function(nc,var.name,lonc,latc,save.file) {

   lon <- ncvar_get(nc,'lon')   
   lat <- ncvar_get(nc,'lat')   

   lon.ix <- which.min(abs(lonc-lon))
   lat.ix <- which.min(abs(latc-lat))
print(lon.ix)
print(lat.ix)

   cell <- ncvar_get(nc,var.name,start=c(lon.ix,lat.ix,1),count=c(1,1,-1))

   if (var.name=='pr' & grepl('_CanRCM4_BC_',save.file)) {
      cell <- cell * 86400
   }
   time.series <- as.Date(as.character(netcdf.calendar(nc)))

   rv <- list(time=time.series,data=cell)
   return(rv)
}

##----------------------------------------------------

save_series <- function(var.name,file.name,read.dir,
                        save.file,
                        lonc,latc,site.name) {

   nc <- nc_open(paste0(read.dir,file.name))
   var.series <- extract_series(nc,var.name,lonc,latc,save.file)   
   nc_close(nc)
   save(var.series,file=save.file)

   return(var.series)   
}

##----------------------------------------------------
make_cell_series <- function(read.suffix,read.dir,
                             save.suffix,save.dir,
                             lonc,latc,site.name) {
                                                   
   can.series <- vector(mode='list',length=4)
   names(can.series) <- c('pr','tasmax','tasmin','tas')
   var.list <- c('pr','tasmax','tasmin')
   for (var.name in var.list) {
      rcm.file <- paste0(var.name,read.suffix)
      save.file <- paste0(save.dir,var.name,save.suffix)

      can.series[[var.name]] <- save_series(var.name,file.name=rcm.file,read.dir=read.dir,
                                            save.file=save.file,
                                            lonc,latc,site.name)
   }
   ##Add average temperature
   tas.save.file <- paste0(save.dir,"tas",save.suffix)
   can.tas <- (can.series$tasmax$data + can.series$tasmin$data) / 2
   var.series <- list(time=can.series$tasmax$time,data=can.tas)
   save(var.series,file=tas.save.file)         

}


##****************************************************

##site.name <- 'Kelowna'
##lonc <- -119.4
##latc <- 49.9

sites <- list(list(site.name="Fort_Nelson",lonc=-122.71,latc=58.79),
              list(site.name="Skagway",lonc=-135.0,latc=59.25),
              list(site.name="Prince_George",lonc=-122.75,latc=53.92), 
              list(site.name="Prince_Rupert",lonc=-130.3,latc=54.3),
              list(site.name="Golden",lonc=-116.99,latc=51.3),
              list(site.name="Ucluelet",lonc=-125.5,latc=48.9),
              list(site.name="Prince_George",lonc=-122.75,latc=53.92),
              list(site.name="Kelowna",lonc=-119.5,latc=49.9),
              list(site.name="Victoria",lonc=-123.4,latc=48.6))


for (site in sites) {

site.name <- site$site.name
print(site.name)
lonc <- site$lonc
latc <- site$latc

save.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/save_dir/",site.name,"/")
if (!file.exists(save.dir)) {
  dir.create(save.dir,recursive=TRUE)
}

bcci.dir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/mvbc/BCCI/"
bccaq.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"

##------------------------------------------------------
##CanRCM4 Acting as observations
if (1==0) {
rcm.suffix <- paste0("_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
save.suffix <- paste0("_day_CanRCM4_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")

make_cell_series(read.suffix=rcm.suffix,read.dir=bccaq.dir,
                 save.suffix=save.suffix,save.dir,
                 lonc,latc,site.name)
}

##------------------------------------------------------
##BCCI Acting as GCM input
if (1==0) {
bcci.suffix <- "_BCCI_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
save.suffix <- paste0("_day_BCCI_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")

make_cell_series(read.suffix=bcci.suffix,read.dir=bcci.dir,
                 save.suffix=save.suffix,save.dir,
                 lonc,latc,site.name)
}
##------------------------------------------------------
##BCCAQv2 for Comparison
if (1==0) {
bccaqv2.suffix <- "_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"

save.suffix <- paste0("_day_BCCAQv2_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")

make_cell_series(read.suffix=bccaqv2.suffix,read.dir=bccaq.dir,
                 save.suffix=save.suffix,save.dir,
                 lonc,latc,site.name)
}

##------------------------------------------------------
##MBCn from completed downscaling 

mbcn.dir <- "/storage/data/climate/downscale/BCCAQ2+PRISM/mvbc/MBCn/"
mbcn.suffix <- "_MBCn_iterated_trace_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc"
save.suffix <- paste0("_day_MBCn_iterated_trace_BC_RCM-Grid_",site.name,"_cell_1950-2100.RData")

make_cell_series(read.suffix=mbcn.suffix,read.dir=mbcn.dir,
                 save.suffix=save.suffix,save.dir,
                 lonc,latc,site.name)



}
