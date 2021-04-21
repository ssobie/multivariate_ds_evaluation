##Run BCCAQv2 using ClimDown

library(ClimDown)


##---------------------------------
testing <- FALSE

if (testing) {

   tmpdir <- '/local_temp/ssobie/bccaqc2/'
   varname <- 'tasmax'
} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
   }
}



tmpdir <- paste0(tmpdir,'/bccaqv2_CanRCM4_',varname,'_tmp')
if (!file.exists(tmpdir)) {
   dir.create(tmpdir,recursive=TRUE)
}

var.name <- varname
obs.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"
gcm.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/gcm_grid/"
bccaq.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"


obs.file <- paste0(obs.dir,"pr_BC_RCM-Grid_CanRCM4_day_19500101-21001231.nc")
###              var.name,"_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
gcm.file <- paste0(gcm.dir,"pr_BC_GCM-Grid_CanRCM4_day_19500101-21001231.nc")
###              var.name,"_BC_GCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")

bccaq.file <- gsub(paste0(var.name,"_BC_GCM-Grid"),paste0(var.name,'_BCCAQv2_BC_RCM-Grid'),gcm.file) 

gcm.tmp <- paste0(tmpdir,'/',basename(gcm.file))
obs.tmp <- paste0(tmpdir,'/',basename(obs.file))
bccaq.tmp <- paste0(tmpdir,'/',basename(bccaq.file))

file.copy(from=obs.file,to=tmpdir)
file.copy(from=gcm.file,to=tmpdir)
print('Copying Complete')

options(max.GB=1)
options(calibration.start=as.POSIXct('1951-01-01', tz='GMT'))
options(calibration.end=as.POSIXct('1980-12-31', tz='GMT'))

bccaq.netcdf.wrapper(gcm.file=gcm.tmp,obs.file=obs.tmp,out.file=bccaq.tmp,varname=var.name)
##ca.netcdf.wrapper(gcm.file=gcm.tmp,obs.file=obs.tmp,varname=var.name)
print('BCCAQv2 Complete')

file.copy(from=bccaq.tmp,to=bccaq.dir)


