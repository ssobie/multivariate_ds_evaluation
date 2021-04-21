##Make plots comparing a reference and simulated dataset for seasonal climatologies 

source('/storage/home/ssobie/code/repos/multivariate_ds_evaluation/map_clim_comparisons_pnw.r',chdir=T)

##---------------------------------------------------------------------
##*********************************************************************

Rprof('plotting.out')

##-------------------------------------------------------------------

plot.dir <- '/storage/data/projects/rci/data/stat.downscaling/MBCn/CanRCM4/'
var.name <- 'pr'

past.int <- '1951-1980'
proj.intervals <- c('1981-2010','2011-2040','2041-2070')
sim <- 'BCCAQv2'

main.title <- switch(var.name,
                     tasmax=paste0('Annual Maximum Temperature Comparison (CanRCM4, ',sim,')'),
                     tasmin=paste0('Annual Minimum Temperature Comparison (CanRCM4, ',sim,')'),            
                     pr=paste0('Annual Maximum Precipitation Comparison (CanRCM4, ',sim,')'))

leg.title <- switch(var.name,
                    pr='mm',    
                    tasmax='\u00B0C',
                    tasmin='\u00B0C')


run <- "r1i1p1"
scenario <- "rcp85"


var.type <- switch(var.name,
                   pr='maximum',
                   tasmax='maximum',
                   tasmin='minimum')

for (proj.int in proj.intervals) {

##--------------------------------
##Reference Data

ref.name <- ref <- 'CanRCM4'
ref.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",ref,"_Derived/annual_extremes/climatologies/")
past.ref.file <- paste0(var.name,"_annual_",var.type,"_climatology_",ref,"_historical+",scenario,"_",run,"_",past.int,".nc")
past.ref.brick <- subset(brick(paste0(ref.dir,past.ref.file)),1)

proj.ref.file <- paste0(var.name,"_annual_",var.type,"_climatology_",ref,"_historical+",scenario,"_",run,"_",proj.int,".nc")
proj.ref.brick <- subset(brick(paste0(ref.dir,proj.ref.file)),1)

##--------------------------------
##Simulation Data

sim.name <- sim 
sim.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",sim,"_Derived/annual_extremes/climatologies/")
past.sim.file <- paste0(var.name,"_annual_",var.type,"_climatology_",sim,"_historical+",scenario,"_",run,"_",past.int,".nc")
past.sim.brick <- subset(brick(paste0(sim.dir,past.sim.file)),1)

proj.sim.file <- paste0(var.name,"_annual_",var.type,"_climatology_",sim,"_historical+",scenario,"_",run,"_",proj.int,".nc")
proj.sim.brick <- subset(brick(paste0(sim.dir,proj.sim.file)),1)


plot.file <- paste0('annual_extremes.',var.name,'.',var.type,'.clim.comparisons.',ref,'.',sim,'.',past.int,'.',proj.int,'.png')
write.file <- paste0(plot.dir,plot.file)

if (var.name=='pr') {
   clim.breaks <- seq(0,100,10)
   diff.breaks <- seq(-12.5,10,2.5)
   anom.breaks <- seq(-7.5,20,2.5)
   anom.diff.breaks <- diff.breaks ##seq(-12.5,15,2.5)
}

if (var.name=='tasmax') {
   clim.breaks <- seq(0,40,5)
   diff.breaks <- seq(-3,4,0.5)
   anom.breaks <- seq(-1,8,1)
   anom.diff.breaks <- diff.breaks ##seq(-3,3,0.5)
}

if (var.name=='tasmin') {
    clim.breaks <- seq(-50,10,5)
    diff.breaks <- seq(-6,3,1)
    anom.breaks <- seq(-2,15,2)
    anom.diff.breaks <- diff.breaks ##seq(-6,3,1)
}


make_comparison_plot(plot.file,var.name,main.title,leg.title,
                     ref.name,past.ref.brick,proj.ref.brick,
                     sim.name,past.sim.brick,proj.sim.brick,
                     past.int,proj.int,
                     clim.breaks=clim.breaks,
                     diff.breaks=diff.breaks,
                     anom.breaks=anom.breaks,
                     anom.diff.breaks=anom.diff.breaks)


}

Rprof(NULL)