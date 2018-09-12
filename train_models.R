## This script runs 20 parameterization routines for each model in the
## phenor package and saves optimal parameter sets and RMSE for each run

library(phenor)
library(jsonlite)
library(foreach)
library(iterators)
library(doParallel)
library(landsat)

comb <- function(x, ...) {
  lapply(seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

#Register the parallel backend
registerDoParallel(16)

args = commandArgs(trailingOnly=T)
M = args[1]
#M = 'M1_Jan1_5C'

# Load phenocam transition dates and associated Daymet time series
#load("~/Code/GitHub/phenor/data/phenocam_DB.rda") #old version from KH phenor package
load("/projectnb/modislc/users/emelaas/MSB/PhenoCam/phenocam_DB.RData")

for (i in 1:138){
  phenocam_DB[[i]]$site <- names(phenocam_DB)[i]
}

# modify silas little 2017 DOY
phenocam_DB[[105]]$transition_dates[5] <- 134
data <- phenocam_DB

# concat locations data into a matrix with the first row
# being the latitude and the second longitude
location = do.call("cbind",lapply(data,function(x){
  if(!is.null(x)){
    matrix(rep(x$location, ncol(x$Ti)), 2, ncol(x$Ti))
  }
}))

# concat sitenames into a vector using a do.call()
site = as.character(do.call("c",lapply(data, function(x){
  if(!is.null(x)){
    rep(x$site, ncol(x$Ti))
  }
})))

site_ind <- match(site,names(phenocam_DB))

# Find all sites located west of 100 deg W and eliminate them
w <- which(location[2,] < -100)
site_w <- as.numeric(names(table(site_ind[w])))
phenocam_DB[site_w] <- NULL


# # Replace existing SOS 25% with SOS 50% dates
# setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/PhenoCam/')
# local_database = jsonlite::fromJSON("site_information.json")
# w = which(local_database$veg_type=='DB')
# local_database = local_database[w,]
# 
# #Read in PhenoCam site-metrics
# setwd('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/PhenoCam/ORNL/PhenoCam_V1_1511/PhenoCam_V1_1511/data/')
# phenocam_files = list.files(path=getwd(),pattern=glob2rx("*DB*3day_tr*csv"),full.names=T,include.dirs=T,recursive=TRUE)
# 
# for (i in 1:length(phenocam_files)){
#   if (i == 1){
#     phen = read.table(phenocam_files[i],header=TRUE,
#       sep = ",",
#       skip = 16)
#   } else {
#     tmp = read.table(phenocam_files[i],header=TRUE,
#       sep = ",",
#       skip = 16)
#     phen = rbind(phen,tmp)
#   }
# }
# 
# # only retain the gcc_90 values (the other time series are non-standard)
# phen = phen[which(phen$gcc_value == "gcc_90"),]
# phen$year = format(as.Date(phen$transition_10),"%Y")
# phen$doy_transition_10 = as.numeric(format(as.Date(phen$transition_10),"%j"))
# phen$doy_transition_25 = as.numeric(format(as.Date(phen$transition_25),"%j"))
# phen$doy_transition_50 = as.numeric(format(as.Date(phen$transition_50),"%j"))
# w = which(phen$direction=='rising')
# phen = phen[w,]
# phen = data.frame(phen$sitename,phen$roi_id,phen$year,phen$doy_transition_10,phen$doy_transition_25,phen$doy_transition_50)
# colnames(phen) = c('sites','roi-id','year','s10','s25','s50')
# 
# nam = names(phenocam_data)
# for (i in 1:length(nam)){
#   w = which(phen$sites == nam[i])
#   yrs = phenocam_data[[i]]$year
#   w2 = which(phen$year[w] %in% yrs)
#   phenocam_data[[i]]$transition_dates = phen$s50[w[w2]]
# }

# Loop through each model and run 20 model calibration routines
mod_cal <- foreach(i = 1:20, .combine='comb', .multicombine=TRUE,
  .init=list(list(), list(), list(), list())) %dopar% {

    print(i)

    cal = model_calibration(model = M,
      data = phenocam_DB,
      control = list(max.call = 1e5, temperature = 10000),
      par_ranges = sprintf("%s/extdata/parameter_ranges.csv",path.package("phenor")),
      random_seed = round(runif(1,1,1000)),
      plot = FALSE)

    list(cal$par,cal$rmse,cal$rmse_null,cal$aic$AIC)
  }

# Save model calibration results
save(mod_cal,file=paste('/projectnb/modislc/users/emelaas/scratch15/NASA_TE/Daymet/mod_cal_',M,sep = ''))
