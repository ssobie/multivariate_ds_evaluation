##Make plots comparing a reference and simulated dataset for seasonal climatologies 

source('/storage/home/ssobie/code/repos/multivariate_ds_evaluation/map_clim_comparisons_pnw.r',chdir=T)

##---------------------------------------------------------------------
##*********************************************************************

##-------------------------------------------------------------------

plot.dir <- '/storage/data/projects/rci/data/stat.downscaling/MBCn/CanRCM4/'
var.name <- 'pr'
run <- "r1i1p1"
scenario <- "rcp85"
sim <- 'BCCAQv2'

season <- 'summer'
seas.ix <- 3


var.type <- switch(var.name,
                   pr='total',
                   tasmax='average',
                   tasmin='average')

var.title <- switch(var.name,
                     tasmax=paste0(' ',tools::toTitleCase(var.type),' Maximum Temperature Comparison (CanRCM4, ',sim,')'),
                     tasmin=paste0(' ',tools::toTitleCase(var.type),' Minimum Temperature Comparison (CanRCM4, ',sim,')'),
                     pr=paste0(' Total Precipitation Comparison (CanRCM4, ',sim,')'))

leg.title <- switch(var.name,
                    pr='mm',
                    tasmax='\u00B0C',
                    tasmin='\u00B0C')

main.title <- paste0(tools::toTitleCase(season),var.title)


past.int <- '1951-1980'
proj.intervals <- '2041-2070' ##c('1981-2010','2011-2040','2041-2070')

for (proj.int in proj.intervals) {

##--------------------------------
##Reference Data

ref.name <- ref <- 'CanRCM4'
ref.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",ref,"_Derived/seasonal/climatologies/")
past.ref.file <- paste0(var.name,"_seasonal_",var.type,"_climatology_",ref,"_historical+",scenario,"_",run,"_",past.int,".nc")
past.ref.brick <- subset(brick(paste0(ref.dir,past.ref.file)),seas.ix)

proj.ref.file <- paste0(var.name,"_seasonal_",var.type,"_climatology_",ref,"_historical+",scenario,"_",run,"_",proj.int,".nc")
proj.ref.brick <- subset(brick(paste0(ref.dir,proj.ref.file)),seas.ix)

##--------------------------------
##Simulation Data

sim.name <- sim 
sim.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",sim,"_Derived/seasonal/climatologies/")
past.sim.file <- paste0(var.name,"_seasonal_",var.type,"_climatology_",sim,"_historical+",scenario,"_",run,"_",past.int,".nc")
past.sim.brick <- subset(brick(paste0(sim.dir,past.sim.file)),seas.ix)

proj.sim.file <- paste0(var.name,"_seasonal_",var.type,"_climatology_",sim,"_historical+",scenario,"_",run,"_",proj.int,".nc")
proj.sim.brick <- subset(brick(paste0(sim.dir,proj.sim.file)),seas.ix)


plot.file <- paste0(season,'.',var.name,'.seasonal.clim.comparisons.',ref,'.',sim,'.',past.int,'.',proj.int,'.png')
write.file <- paste0(plot.dir,plot.file)

if (var.name=='pr') {

   if (season=='winter') {
      clim.breaks <- seq(0,1200,100)
      diff.breaks <- seq(-40,80,10)
      anom.breaks <- seq(-40,200,20)
      anom.diff.breaks <- diff.breaks ##seq(-12.5,15,2.5)
   }
   if (season=='summer') {
      clim.breaks <- seq(0,600,50)
      diff.breaks <- seq(-40,60,10)
      anom.breaks <- seq(-80,80,20)
      anom.diff.breaks <- diff.breaks ##seq(-12.5,15,2.5)
   }


}

if (var.name=='tasmax') {
   clim.breaks <- seq(0,35,5)
   diff.breaks <- seq(-3,4,0.5)
   anom.breaks <- seq(0,8,1)
   anom.diff.breaks <- diff.breaks ##seq(-3,3,0.5)
}

if (var.name=='tasmin') {
    clim.breaks <- seq(-30,10,5)
    diff.breaks <- seq(-2,1,0.5)
    anom.breaks <- seq(0,8,1)
    anom.diff.breaks <- diff.breaks
}


make_comparison_plot(plot.file,var.name,main.title,leg.title,
                     ref.name,past.ref.brick,proj.ref.brick,
                     sim.name,past.sim.brick,proj.sim.brick,
                     past.int,proj.int,
                     clim.breaks=clim.breaks,
                     diff.breaks=diff.breaks,
                     anom.breaks=anom.breaks,
                     anom.diff.breaks=anom.diff.breaks)

browser()

}