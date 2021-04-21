##Make plots comparing a reference and simulated dataset for seasonal climatologies 

source('/storage/home/ssobie/code/repos/multivariate_ds_evaluation/map_clim_comparisons_pnw.r',chdir=T)

##---------------------------------------------------------------------
##*********************************************************************

##-------------------------------------------------------------------

plot.dir <- '/storage/data/projects/rci/data/stat.downscaling/MBCn/CanRCM4/'
run <- "r1i1p1"
scenario <- "rcp85"

sim <- 'BCCAQv2'
var.name <- 'tasmax'
pctl <- 'q90'
season <- 'summer'
seas.ix <- 3

var.title <- switch(var.name,
                     tasmax=paste0(' Pr (Q75), ',tools::toTitleCase(var.name),'(',toupper(pctl),') Tencer Compound Events Ratios (CanRCM4, ',sim,')'),
                     tasmin=paste0(' Pr,(Q75), ',tools::toTitleCase(var.name),'(',toupper(pctl),') Tencer Compound Events Ratios (CanRCM4, ',sim,')'))
                    
leg.title <- switch(var.name,
                    pr='Ratio',
                    tasmax='Ratio',
                    tasmin='Ratio')

main.title <- paste0(tools::toTitleCase(season),var.title)


past.int <- '1951-1980'
proj.intervals <- '1981-2010' ##c('1981-2010','2011-2040','2041-2070')

for (proj.int in proj.intervals) {

##--------------------------------
##Reference Data

ref.name <- ref <- 'CanRCM4'
ref.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",ref,"_Derived/tencer_compound/")
##past.ref.file <- paste0("pr_",var.name,"_",pctl,"_ratio_BC_RCM-Grid_",ref,"_day_",past.int,".nc")
past.ref.file <- paste0("pr_",var.name,"_",pctl,"_seasonal_average_compound_days_over_q90_BC_RCM-Grid_",ref,"_day_",past.int,".nc")
past.ref.brick <- subset(brick(paste0(ref.dir,past.ref.file)),seas.ix)

##proj.ref.file <- paste0("pr_",var.name,"_",pctl,"_ratio_BC_RCM-Grid_",ref,"_day_",proj.int,".nc") 
proj.ref.file <- paste0("pr_",var.name,"_",pctl,"_seasonal_average_compound_days_over_q90_BC_RCM-Grid_",ref,"_day_",proj.int,".nc")
proj.ref.brick <- subset(brick(paste0(ref.dir,proj.ref.file)),seas.ix)

##--------------------------------
##Simulation Data

sim.name <- sim 
sim.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",sim,"_Derived/tencer_compound/")
##past.sim.file <- paste0("pr_",var.name,"_",pctl,"_ratio_BC_RCM-Grid_",sim,"_day_",past.int,".nc")
past.sim.file <- paste0("pr_",var.name,"_",pctl,"_seasonal_average_compound_days_over_q90_BC_RCM-Grid_",sim,"_day_",past.int,".nc")
past.sim.brick <- subset(brick(paste0(sim.dir,past.sim.file)),seas.ix)

##proj.sim.file <- paste0("pr_",var.name,"_",pctl,"_ratio_BC_RCM-Grid_",sim,"_day_",proj.int,".nc")
proj.sim.file <- paste0("pr_",var.name,"_",pctl,"_seasonal_average_compound_days_over_q90_BC_RCM-Grid_",sim,"_day_",proj.int,".nc")
proj.sim.brick <- subset(brick(paste0(sim.dir,proj.sim.file)),seas.ix)


plot.file <- paste0(season,'.pr.',var.name,'.',pctl,'.tencer.compound.counts.',ref,'.',sim,'.',past.int,'.',proj.int,'.png')
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
   clim.breaks <- seq(0,7,0.5)
   diff.breaks <- seq(-1.5,2,0.5)
   anom.breaks <- seq(-1,3,0.5)
   anom.diff.breaks <- diff.breaks ##seq(-3,3,0.5)
}

if (var.name=='tasmin') {
    clim.breaks <- seq(-0.25,1.5,0.25)
    diff.breaks <- seq(-0.5,1,0.25)
    anom.breaks <- seq(-1.5,0.5,0.25)
    anom.diff.breaks <- diff.breaks
}


make_comparison_plot(plot.file,var.name,main.title,leg.title,
                     ref.name,past.ref.brick,proj.ref.brick,
                     sim.name,past.sim.brick,proj.sim.brick,
                     past.int,proj.int,bp=1,
                     clim.breaks=clim.breaks,
                     diff.breaks=diff.breaks,
                     anom.breaks=anom.breaks,
                     anom.diff.breaks=anom.diff.breaks)
 

}