
library(ncdf4)

var.name <- "pr"

input.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"
input.file <- paste0(var.name,"_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")

base.file <- paste0(var.name,"_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")

bnc <- nc_open(paste0(input.dir,base.file))
nc <- nc_open(paste0(input.dir,input.file),write=TRUE)

global.atts <- ncatt_get(bnc,varid=0)

global.names <- names(global.atts)
for (g in 1:length(global.names)) {
  ncatt_put(nc,varid=0,attname=global.names[g],attval=global.atts[[g]])
}

nc_close(bnc)
nc_close(nc)