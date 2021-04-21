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

intervals <- c('1951-1980','1981-2010','2011-2040','2041-2070','2070-2099')
##interval <- '2041-2070'

var.name <- 'tasmin'
pctl <- 'q10'
compare <- 'under'
seasons <- c('winter','spring','summer','fall')
##seas.ix <- 4

mean.list <- vector(mode='list',length=4)
names(mean.list) <- seasons

rmse.list <- bias.list <- mean.list


for (i in 1:4) {
   seas.ix <- i
   season <- seasons[i]
   print(season)
   mean.matrix <- matrix(NA,nrow=3,ncol=length(intervals))
   bias.matrix <- matrix(NA,nrow=2,ncol=length(intervals))
   rmse.matrix <- matrix(NA,nrow=2,ncol=length(intervals))

for (j in seq_along(intervals)) {

   interval <- intervals[j]
   print(interval)
var.title <- switch(var.name,
                     tasmax=paste0(' Pr (Q75), ',tools::toTitleCase(var.name),'(',toupper(pctl),') Tencer Compound Events Counts'),
                     tasmin=paste0(' Pr,(Q75), ',tools::toTitleCase(var.name),'(',toupper(pctl),') Tencer Compound Events Counts'))
                    
leg.title <- switch(var.name,
                    pr='Count',
                    tasmax='Count',
                    tasmin='Count')

main.title <- paste0(tools::toTitleCase(season),var.title)



##--------------------------------
##Reference Data

ref.name <- ref <- 'CanRCM4'
ref.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",ref,"_Derived/tencer_compound/")
past.ref.file <- paste0("pr_",var.name,"_seasonal_total_compound_days_",compare,"_",pctl,"_BC_RCM-Grid_",ref,"_day_",interval,".nc")
past.ref.brick <- subset(brick(paste0(ref.dir,past.ref.file))*aplin.mask,seas.ix)

##--------------------------------
##BCCAQv2 Data

sim.name <- 'BCCAQv2'
sim.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/",sim,"_Derived/tencer_compound/")
past.sim.file <- paste0("pr_",var.name,"_seasonal_total_compound_days_",compare,"_",pctl,"_BC_RCM-Grid_",sim,"_day_",interval,".nc")
past.sim.brick <- subset(brick(paste0(sim.dir,past.sim.file))*aplin.mask,seas.ix)

##--------------------------------
##MBCn Data

mbc.name <- 'MBCn'
mbc.dir <- paste0("/storage/data/climate/downscale/RCM/CanRCM4/multivar_evaluation/rcm_grid/MBCn_Derived/tencer_compound/")
past.mbc.file <- paste0("pr_",var.name,"_seasonal_total_compound_days_",compare,"_",pctl,"_BC_RCM-Grid_MBCn_day_",interval,".nc")
past.mbc.brick <- subset(brick(paste0(mbc.dir,past.mbc.file))*aplin.mask,seas.ix)

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



plot.file <- paste0(season,'.pr.',var.name,'.',pctl,'.tencer.compound.totals.CanRCM4.BCCAQv2.MBCn.',interval,'.png')
write.file <- paste0(plot.dir,plot.file)

c.g <- c.l <- d.g <- d.l <- FALSE

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

##Averages
if (var.name=='tasmax') {
   clim.breaks <- seq(0,7,0.5)
   diff.breaks <- seq(-3,3.5,0.5)
   anom.breaks <- seq(-1,3,0.5)
   anom.diff.breaks <- diff.breaks ##seq(-3,3,0.5)
}

##Totals
if (var.name=='tasmax') {
   clim.breaks <-seq(0,200,20)
   c.l <- FALSE
   diff.breaks <- c(-100,seq(-30,30,5),100)
   d.g <- d.l <- TRUE
}

if (var.name=='tasmin') {
    clim.breaks <- seq(0,100,10)
    diff.breaks <- c(-100,seq(-20,20,5),100)
    d.g <- d.l <- TRUE
}


if (1==1) {
make_count_plot(plot.file,var.name,main.title,leg.title,
                 ref.name,sim.name,mbc.name,
                 past.ref.brick,past.sim.brick,past.mbc.brick,
                 interval,bp=0,
                 clim.breaks=clim.breaks,
                 diff.breaks=diff.breaks,
                 c.g=c.g,c.l=c.l,d.g=d.g,d.l=d.l)
##                 anom.breaks=anom.breaks,
##                 anom.diff.breaks=anom.diff.breaks)
}


} 


mean.list[[i]] <- round(mean.matrix,2)
bias.list[[i]] <- round(bias.matrix,2)
rmse.list[[i]] <- round(rmse.matrix,2)

##Save the stats for use in markdown table


}