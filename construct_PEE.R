library(lhs)
library(ncdf4)
library(readr)
library(abind)

setwd("~/Documents/beijing/emissions/MEIC2013_new/")
######################################
####This section generates PPE design####

#load predefined params and ranges
par <- read_csv("PPE_optimised_param.csv")

#draw a latin hypercube sample from uniform distributions as PPE design
unit.design <- maximinLHS(n = nrow(par)*10, k = nrow(par), dup = 5)
summary(unit.design)

#transform to actual parameter ranges
design <- matrix(0, nrow = nrow(unit.design), ncol = ncol(unit.design))
#loop over each col(par) and rescale
for (i in 1:ncol(design)){
  design[,i] <- qunif(unit.design[,i], min = par$min[i], max = par$max[i])
}
colnames(design) <- par$param
summary(design)

#can check that each range are evenly covered
hist(design[,11], xlab=par$param[11], breaks = seq(0.05, 1.5, 0.05))

#2d visualisation of parameter settings
pairs(design, pch=17, cex = 0.5, upper.panel=NULL)

#write.csv(design, "PPE_optimised_design.csv", row.names = F)
##############################################
####This section calculates perturbed tran dirunal profile for each PPE member####

#load PPE design
design <- read.csv("PPE_optimised_design.csv")

#load tran diurnal profile: same for all pollutants (see plot_emissions_full.R)
diurnal <- read.csv("tran_diurnal.csv")

#create data frames to store perturbed profiles
CO <- matrix(0, nrow = nrow(diurnal), ncol = nrow(design))
colnames(CO) <- paste("Run", seq(1, 140, 1), sep = "")

NOX <- matrix(0, nrow = nrow(diurnal), ncol = nrow(design))
colnames(NOX) <- paste("Run", seq(1, 140, 1), sep = "")

#perturb the total weight of 0-5 am (incl.) emissions -> new weight for each hour
for (i in 1:nrow(design)){
  
  CO[c(1:6), i] <- diurnal$CO[c(1:6)] * design$CO.tran.night[i] / sum(diurnal$CO[c(1:6)])
  CO[7:24, i] <- diurnal$CO[7:24] * (1-design$CO.tran.night[i]) / sum(diurnal$CO[7:24])
  
  NOX[c(1:6), i] <- diurnal$NOX[c(1:6)] * design$NOX.tran.night[i] / sum(diurnal$NOX[c(1:6)])
  NOX[7:24, i] <- diurnal$NOX[7:24] * (1-design$NOX.tran.night[i]) / sum(diurnal$NOX[7:24])
  
}

#add Run0
CO <- data.frame(hour = diurnal$hour, Run0 = diurnal$CO, CO)
#write.csv(CO, file = "tran_CO_diurnal_perturbed_optimised.csv", row.names = F)
NOX <- data.frame(hour = diurnal$hour, Run0 = diurnal$NOX, NOX)
#write.csv(NOX, file = "tran_NOX_diurnal_perturbed_optimised.csv", row.names = F)

##############################################################
####This section generates perturbed emissions for each PPE member####

#load raw emissions
ind <- nc_open("meic2013_new_industry_plus1hr.nc")
pow <- nc_open("meic2013_new_power_plus1hr.nc")
res <- nc_open("meic2013_new_residential_plus1hr.nc")
tran <- nc_open("meic2013_new_transport_plus1hr.nc")
all <- nc_open("meic2013_new_all_sectors_plus1hr.nc")

sectors <- c("ind", "pow", "res", "tran")
for (i in sectors){
  x <- get(i)
  y <- ncvar_get(x, "NO2")
  z <- ncvar_get(x, "NOX")
  b <- ncvar_get(x, "CO")
  assign(paste(i, "no2", "raw", sep = "."), y)
  assign(paste(i, "nox", "raw", sep = "."), z)
  assign(paste(i, "co", "raw", sep = "."), b)
  rm(x, y, z, b)
}

#leave VOC unscaled
voc <- ncvar_get(all, "VOC")

#free up memory (ind used later)
rm(pow, res, tran, i)

#load PPE design
design <- read.csv("PPE_optimised_design.csv")

#load perturbed diurnal profiles
diurnal.co <- read.csv("tran_CO_diurnal_perturbed_optimised.csv")
diurnal.nox <- read.csv("tran_NOX_diurnal_perturbed_optimised.csv")

for (i in seq(1,140,1)) {
  
  ind.co.scaled <- abind(ind.co.raw[,,1,]*design$CO.ind.ground[i], 
                         ind.co.raw[,,2:7,]*design$CO.ind.upper[i], along = 3)
  pow.co.scaled <- abind(pow.co.raw[,,1:3,]*design$CO.pow.below.200[i],
                         pow.co.raw[,,4:7,]*design$CO.pow.above.200[i], along = 3)
  res.co.scaled <- res.co.raw*design$CO.res[i]
  tran.co.scaled <- tran.co.raw*design$CO.tran[i]
  
  #create a new array to store diurnal profile modified tran
  tran.co.scaled2 <- tran.co.scaled
  #loop over the time dim, select 24h at a time and apply new weights
  for (j in seq(1, 8761, 24)){
    k <- j+23
    day <- tran.co.scaled[,,1,j:k] 
    #new weight/old weight: percent change 
    day.scaled <- sweep(day, MARGIN = 3, unlist(diurnal.co[i+2]/diurnal.co[2], use.names = F), FUN = "*")
    #24h total remains the same
    #all.equal(sum(day), sum(day.scaled))
    #append 
    tran.co.scaled2[,,1,j:k] <- day.scaled
  }
  rm(j, k, day, day.scaled)
  #all.equal(sum(tran.co.scaled), sum(tran.co.scaled2))
  
  #add up sector-specific emissions  
  co.scaled <- ind.co.scaled + pow.co.scaled + res.co.scaled + tran.co.scaled2
  
  
  #repeat the steps above for NOX
  ind.nox.scaled <- abind(ind.nox.raw[,,1,]*design$NOX.ind.ground[i], 
                          ind.nox.raw[,,2:7,]*design$NOX.ind.upper[i], along = 3)
  pow.nox.scaled <- abind(pow.nox.raw[,,1:3,]*design$NOX.pow.below.200[i], 
                          pow.nox.raw[,,4:7,]*design$NOX.pow.above.200[i], along = 3)
  res.nox.scaled <- res.nox.raw*design$NOX.res[i]
  tran.nox.scaled <- tran.nox.raw*design$NOX.tran[i]
  
  tran.nox.scaled2 <- tran.nox.scaled
  for (j in seq(1, 8761, 24)){
    k <- j+23
    day <- tran.nox.scaled[,,1,j:k] 
    day.scaled <- sweep(day, MARGIN = 3, unlist(diurnal.nox[i+2]/diurnal.nox[2], use.names = F), FUN = "*")
    tran.nox.scaled2[,,1,j:k] <- day.scaled
  }
  rm(j, k, day, day.scaled)
  #all.equal(sum(tran.nox.scaled), sum(tran.nox.scaled2))
  nox.scaled <- ind.nox.scaled + pow.nox.scaled + res.nox.scaled + tran.nox.scaled2
  
  
  #repeat the steps above for NO2: use same scalings for NOX (fNO2 remains the same)
  ind.no2.scaled <- abind(ind.no2.raw[,,1,]*design$NOX.ind.ground[i], 
                          ind.no2.raw[,,2:7,]*design$NOX.ind.upper[i], along = 3)
  pow.no2.scaled <- abind(pow.no2.raw[,,1:3,]*design$NOX.pow.below.200[i],
                          pow.no2.raw[,,4:7,]*design$NOX.pow.above.200[i], along = 3)
  res.no2.scaled <- res.no2.raw*design$NOX.res[i]
  tran.no2.scaled <- tran.no2.raw*design$NOX.tran[i]
  
  tran.no2.scaled2 <- tran.no2.scaled
  for (j in seq(1, 8761, 24)){
    k <- j+23
    day <- tran.no2.scaled[,,1,j:k] 
    day.scaled <- sweep(day, MARGIN = 3, unlist(diurnal.nox[i+2]/diurnal.nox[2], use.names = F), FUN = "*")
    tran.no2.scaled2[,,1,j:k] <- day.scaled
  }
  rm(j, k, day, day.scaled)
  no2.scaled <- ind.no2.scaled + pow.no2.scaled + res.no2.scaled + tran.no2.scaled2
  
  
  #define dimension vars for a new .nc file using those in ind
  x <- ncdim_def("x", units = ind$dim$x$units, vals = ind$dim$x$vals)
  y <- ncdim_def("y", units = ind$dim$y$units, vals = ind$dim$y$vals)
  z <- ncdim_def("z", units = ind$dim$z$units, vals = ind$dim$z$vals)
  time <- ncdim_def("time", units = ind$dim$time$units, vals = ind$dim$time$vals, unlim = T)
  
  #define 4D vars 
  CO <- ncvar_def("CO", units = ind$var$CO$units, dim = list(x, y, z, time),
                  longname = "gridded emissions of CO", prec = "float")
  NOX <- ncvar_def("NOX", units = ind$var$NOX$units, dim = list(x, y, z, time),
                   longname = "gridded emissions of NOX", prec = "float")
  NO2 <- ncvar_def("NO2", units = ind$var$NO2$units, dim = list(x, y, z, time),
                   longname = "gridded emissions of NO2", prec = "float")
  VOC <- ncvar_def("VOC", units = all$var$VOC$units, dim = list(x, y, z, time),
                   longname = "gridded emissions of VOC", prec = "float")
  
  #create the new .nc file and assign values 
  new <- nc_create("new.nc", list(CO, NOX, NO2, VOC))
  ncvar_put(new, varid = "CO", vals = co.scaled)
  ncvar_put(new, varid = "NOX", vals = nox.scaled)
  ncvar_put(new, varid = "NO2", vals = no2.scaled)
  ncvar_put(new, varid = "VOC", vals = voc)
  nc_close(new)
  
  #run bash script to edit meta data (need to allow to execute first)
  system("./PPE_edit_meta.sh")
  
  #rename file
  filename <- paste("meic2013_new_all_sectors_PPE_optimised_Run", i, ".nc", sep = "")
  system2("mv", args = c("new.nc", filename))
  
  #free up memory
  rm(list = ls(pattern = "scaled"))
}



