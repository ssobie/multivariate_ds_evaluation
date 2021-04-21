##Calculates precipitation as snow matching how Plan2Adapt calculates that variable

source('/storage/data/projects/rci/stat.downscaling/bccaq2/code/new.netcdf.calendar.R',chdir=T)

##----
library(foreach)
library(doParallel)
registerDoParallel(cores=2) 

library(ncdf4)
library(PCICt)


##--------------------------------------------------------------

make_pas_file <- function(gcm,scenario,run,
                          gcm.file,tmp.dir) {

   gcm.nc <- nc_open(paste0(tmp.dir,gcm.file))
   gcm.series <- netcdf.calendar(gcm.nc)
   gcm.time <- ncvar_get(gcm.nc,'time')
   gcm.cal <- attr(gcm.series,'cal')
   nc_close(gcm.nc)
   
   pas.file <- gsub('pr_','pas_',gcm.file)

   file.copy(from=paste0(tmp.dir,gcm.file),to=paste0(paste0(tmp.dir,pas.file)))
   Sys.sleep(1)

   work <- paste0('ncrename -v pr,pas ',tmp.dir,pas.file)
   system(work)
   Sys.sleep(1)

   print('Created PAS file')
   nc <- nc_open(paste0(tmp.dir,pas.file),write=TRUE)

   ncatt_put(nc,'pas','standard_name','Precipitation-as-snow')
   ncatt_put(nc,'pas','long_name','Precipitation-as-snow')
   ncatt_put(nc,0,'history','')
   ncatt_put(nc,0,'method','PAS = (Pr | TAS <= 0), same as Plan2Adapt')

   nc_close(nc)
   return(pas.file)

}

##--------------------------------------------------------------
calculate_pas <- function(tas,prec) {
   temp <- tas > 0
   pas <- prec
   pas[temp] <- 0   

   return(pas)
}

##--------------------------------------------------------------

pas_for_model <- function(pas.ncs,lat.ix,n.lon,flag,
                          pr.list,tas.list) {
   ##Variables
   flen <- sum(!flag)
   n.col <- length(pr.list[[1]])
   pas.matrix <- matrix(NA,nrow=n.lon,ncol=n.col)

   if (flen!=0) { ##Some Real Values
      pr.sub.list <- pr.list[!flag]
      rm(pr.list)
      tas.sub.list <- tas.list[!flag]
      rm(tas.list)

      sub.values <- foreach(
                        tas=tas.sub.list,
                        prec=pr.sub.list,
                        .export='calculate_pas'
                        ) %dopar% {
                              pas.values <- calculate_pas(tas,prec)
                        }

      sub.matrix <- matrix(unlist(sub.values),nrow=flen,ncol=n.col,byrow=TRUE)

      rm(tas.sub.list)
      rm(pr.sub.list)

      rm(sub.values)
      pas.matrix[!flag,] <- sub.matrix
      rm(sub.matrix)    
     
   } else {
      print('All NA values')
   }

   ncvar_put(pas.ncs,varid='pas',vals=pas.matrix,
             start=c(1,lat.ix,1),count=c(-1,1,-1))

   rm(pas.matrix)
   gc()
}



##--------------------------------------------------------------

##****************************************************************
testing <- TRUE

res <- NULL
if (testing) {
###   tmpdir <- '/local_temp/ssobie'
###    gcm <- 'ACCESS1-0'
###   scenario <- 'rcp45'
###   run <- 'r1i1p1'
   res <- NULL

##   args <- commandArgs(trailingOnly=TRUE)
##   for(i in 1:length(args)){
##       eval(parse(text=args[[i]]))
##   }


   base.dir <- gcm.dir <- "/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/"
   type <- 'MBCn'
   
   file.end <- switch(type,
                      CanRCM4="BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc",
                      BCCAQv2="BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc",
                      MBCn="MBCn_iterated_trace_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_19500101-21001231.nc")
   pr.file <- paste0("pr_",file.end)
   tasmax.file <- paste0("tasmax_",file.end)
   tasmin.file <- paste0("tasmin_",file.end)
   write.dir <- paste0(base.dir,type,"_Derived/pas/")
   tmpdir <- '/local_temp/ssobie'
   if (!file.exists(write.dir)) {
      dir.create(write.dir,recursive=TRUE)
   }

} else {
   args <- commandArgs(trailingOnly=TRUE)
   for(i in 1:length(args)){
       eval(parse(text=args[[i]]))
   }

   pr.file <- prfile
   tasmax.file <- txfile
   tasmin.file <- tnfile
   gcm.dir <- gcmdir
   write.dir <- writedir
}

##-----------------------------------------------

tmp.dir <- paste0(tmpdir,'/',type,'_pas/')
if (!file.exists(tmp.dir)) {
  dir.create(tmp.dir,recursive=TRUE)
}

print('Beginning to copy files')
file.copy(from=paste0(gcm.dir,pr.file),to=tmp.dir)
print('Copied precipitation')
file.copy(from=paste0(gcm.dir,tasmax.file),to=tmp.dir)
print('Copied maximum temperature')
file.copy(from=paste0(gcm.dir,tasmin.file),to=tmp.dir)
print('Copied minimum temperature')
print('Done copying daily files')



##---------------------------------------------------------------------------
##PAS Files for writing

pas.file <-  make_pas_file(gcm,scenario,run,
                           pr.file,tmp.dir)

pas.ncs <- nc_open(paste0(tmp.dir,pas.file),write=TRUE)


##---------------------------------------------------------------------------
##---------------------------------------------------------------------------


print('Data opening')
pr.nc <- nc_open(paste0(tmp.dir,pr.file),write=FALSE)
pr.dates <- netcdf.calendar(pr.nc)

tx.nc <- nc_open(paste0(tmp.dir,tasmax.file),write=FALSE)
tx.dates <- netcdf.calendar(tx.nc)
tn.nc <- nc_open(paste0(tmp.dir,tasmin.file),write=FALSE)
tn.dates <- netcdf.calendar(tn.nc)

if (length(pr.dates) != length(tx.dates) | length(pr.dates) != length(tn.dates)) {
###browser()
###   stop('Input data is not all the same length')

}

scale <- 1
if (type=='CanRCM4') {scale <- 86400}

lon <- ncvar_get(pr.nc,'lon')
lat <- ncvar_get(pr.nc,'lat')
n.lon <- length(lon)
n.lat <- length(lat)

for (j in 1:n.lat) { 
   print(paste0('Latitude: ',j,' of ',n.lat))
   lat.ix <- j
   tasmax.subset <- ncvar_get(tx.nc,'tasmax',start=c(1,j,1),count=c(-1,1,-1))
   tasmin.subset <- ncvar_get(tn.nc,'tasmin',start=c(1,j,1),count=c(-1,1,-1))
   tas.subset <- (tasmax.subset + tasmin.subset) / 2    
   rm(tasmax.subset)
   rm(tasmin.subset)

   flag <- is.na(tas.subset[,1])
   tas.list <- lapply(seq_len(nrow(tas.subset)), function(k) tas.subset[k,])
   rm(tas.subset)

   pr.subset <- ncvar_get(pr.nc,'pr',start=c(1,j,1),count=c(-1,1,-1))*scale
   pr.list <- lapply(seq_len(nrow(pr.subset)), function(k) pr.subset[k,])
   rm(pr.subset)

   ##Calculate PAS
   pas_for_model(pas.ncs,lat.ix,n.lon,flag,
                  pr.list,tas.list)

   rm(pr.list)
   rm(tas.list)      

}##Latitude Loop

nc_close(pr.nc)
nc_close(tx.nc)
nc_close(tn.nc)
nc_close(pas.ncs)

file.copy(from=paste0(tmp.dir,pas.file),to=write.dir,overwrite=TRUE)
Sys.sleep(1)

file.remove(paste0(tmp.dir,pr.file))
file.remove(paste0(tmp.dir,tasmax.file))
file.remove(paste0(tmp.dir,tasmin.file))
file.remove(paste0(tmp.dir,pas.file))
