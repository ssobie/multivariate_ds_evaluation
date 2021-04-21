##Make plots comparing a reference and simulated dataset for seasonal climatologies 

source('/storage/home/ssobie/code/repos/multivariate_ds_evaluation/map_clim_comparisons_pnw.r',chdir=T)

##---------------------------------------------------------------------
##*********************************************************************

##-------------------------------------------------------------------

plot.dir <- '/storage/data/projects/rci/data/stat.downscaling/MBCn/CanRCM4/'
var.name <- 'tasmin'
run <- "r1i1p1"
scenario <- "rcp85"
sim <- 'BCCAQv2'

leg.title <- 'Corr.'

seasons <- c('winter','spring','summer','fall','annual')
var.title <- switch(var.name,
                    pr='Precipitation',
                    tasmax='Max. Temperature',
                    tasmin='Min. Temperature')

for (season in seasons) {
print(season)

main.title <- paste0(tools::toTitleCase(season),' ',var.title,' Lag-1 Autocorrelation Comparison (CanRCM4, ',sim,')')

past.int <- '1951-1980'
proj.intervals <- '2041-2070' ##c('1981-2010','2011-2040','2041-2070')

for (proj.int in proj.intervals) {
 

 ##--------------------------------
##Reference Data

ref.name <- ref <- 'CanRCM4'
ref.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",ref,"_Derived/lag_autocorrelation/")
past.ref.file <- paste0(season,'_lag_corr_',var.name,"_BC_",ref,"_",past.int,".nc")
past.ref.brick <- subset(brick(paste0(ref.dir,past.ref.file)),1)

proj.ref.file <- paste0(season,'_lag_corr_',var.name,"_BC_",ref,"_",proj.int,".nc")
proj.ref.brick <- subset(brick(paste0(ref.dir,proj.ref.file)),1)

##--------------------------------
##Simulation Data

sim.name <- sim 
sim.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",sim,"_Derived/lag_autocorrelation/")
past.sim.file <- paste0(season,'_lag_corr_',var.name,"_BC_",sim,"_",past.int,".nc")
past.sim.brick <- subset(brick(paste0(sim.dir,past.sim.file)),1)

proj.sim.file <- paste0(season,'_lag_corr_',var.name,"_BC_",sim,"_",proj.int,".nc")
proj.sim.brick <- subset(brick(paste0(sim.dir,proj.sim.file)),1)


plot.file <- paste0(season,'.',var.name,'.lag1.autocorrelation.comparisons.',ref,'.',sim,'.',past.int,'.',proj.int,'.png')
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
                     past.int,proj.int)

#                     clim.breaks=clim.breaks,
#                     diff.breaks=diff.breaks,
#                     anom.breaks=anom.breaks,
#                     anom.diff.breaks=anom.diff.breaks)

##browser()

}
}