##Make plots comparing a reference and simulated dataset for seasonal climatologies 

source('/storage/home/ssobie/code/repos/multivariate_ds_evaluation/map_count_comparisons_pnw.r',chdir=T)

##---------------------------------------------------------------------
##*********************************************************************

##-------------------------------------------------------------------

plot.dir <- '/storage/data/projects/rci/data/stat.downscaling/MBCn/CanRCM4/'
run <- "r1i1p1"
scenario <- "rcp85"
sim <- 'BCCAQv2'

aplin.mask <- brick("/storage/data/climate/downscale/RCM/CanRCM4/daily/NAM_44/rcm_grid/anusplin.mask.nc")

intervals <- c('1951-1980','1981-2010','2011-2040','2041-2070','2071-2100')

var.name <- 'tasmax'
pctl <- '90'

mean.matrix <- matrix(NA,nrow=3,ncol=length(intervals))
bias.matrix <- matrix(NA,nrow=2,ncol=length(intervals))
rmse.matrix <- matrix(NA,nrow=2,ncol=length(intervals))

for (j in seq_along(intervals)) {

   interval <- intervals[j]
   print(interval)
var.title <- switch(var.name,
                     tasmax=paste0(tools::toTitleCase(var.name),' (',toupper(pctl),') Daily'),
                     tasmin=paste0(tools::toTitleCase(var.name),' (',toupper(pctl),') Daily'))
                    
leg.title <- switch(var.name,
                    pr='mm',
                    tasmax='\u00B0C',
                    tasmin='\u00B0C')

main.title <- var.title



##--------------------------------
##Reference Data

ref.name <- ref <- 'CanRCM4'
ref.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",ref,"_Derived/percentiles/")
past.ref.file <- paste0(var.name,"_quantile_",pctl,"_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_",interval,".nc")
past.ref.brick <- subset(brick(paste0(ref.dir,past.ref.file))*aplin.mask,1)

##--------------------------------
##BCCAQv2 Data

sim.name <- 'BCCAQv2'
sim.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",sim,"_Derived/percentiles/")
past.sim.file <- paste0(var.name,"_quantile_",pctl,"_BCCAQv2_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_",interval,".nc")
past.sim.brick <- subset(brick(paste0(sim.dir,past.sim.file))*aplin.mask,1)

##--------------------------------
##MBCn Data

mbc.name <- 'MBCn'
mbc.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/MBCn_Derived/percentiles/")
past.mbc.file <- paste0(var.name,"_quantile_",pctl,"_MBCn_iterated_BC_RCM-Grid_CCCma-CanESM2_historical-r1_r1i1p1_CCCma-CanRCM4_day_",interval,".nc")
past.mbc.brick <- subset(brick(paste0(mbc.dir,past.mbc.file))*aplin.mask,1)

##-------------------------------
##Save Stats

ref.mean <- cellStats(past.ref.brick,mean)
sim.mean <- cellStats(past.sim.brick,mean)
mbc.mean <- cellStats(past.mbc.brick,mean)

sim.diff <- past.sim.brick - past.ref.brick
sim.bias <- cellStats(sim.diff,mean)
sim.rmse <- sqrt(cellStats(sim.diff^2,mean))

mbc.diff <- past.mbc.brick - past.ref.brick
mbc.bias <- cellStats(mbc.diff,mean)
mbc.rmse <- sqrt(cellStats(mbc.diff^2,mean))

mean.matrix[1,j] <- ref.mean
mean.matrix[2,j] <- sim.mean
mean.matrix[3,j] <- mbc.mean

bias.matrix[1,j] <- sim.bias
bias.matrix[2,j] <- mbc.bias

rmse.matrix[1,j] <- sim.rmse
rmse.matrix[2,j] <- mbc.rmse



plot.file <- paste0(var.name,'.',pctl,'.daily.percentile.CanRCM4.BCCAQv2.MBCn.',interval,'.png')
write.file <- paste0(plot.dir,plot.file)

if (var.name=='pr') {


}

if (var.name=='tasmax') {
   clim.breaks <- seq(0,7,0.5)
   diff.breaks <- seq(-3,3.5,0.5)
}

if (var.name=='tasmin') {
    clim.breaks <- seq(-0.25,1.5,0.25)
    diff.breaks <- seq(-0.5,1,0.25)
}


if (1==1) {
make_count_plot(plot.file,var.name,main.title,leg.title,
                 ref.name,sim.name,mbc.name,
                 past.ref.brick,past.sim.brick,past.mbc.brick,
                 interval,bp=0)
##                 clim.breaks=clim.breaks,
##                 diff.breaks=diff.breaks,
}


} 

