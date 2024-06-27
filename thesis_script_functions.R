# Activating required packages
library(terra)
library(earlywarnings)
library(ggplot2)
library(ggpubr)
## Setting working directory 
#Replace ''your_file_location" with a path to your data folder. Using ctrl+f & replace is recommended for the rest of the script.
setwd("your_file_location")

### RESILLIENCE INDICATOR SCRIPT ###
# Developed for MSc Thesis ''Relations between temporal resilience indicators and trend breakpoints in a dryland high-mountain catchment"
# Developed by Willem Grootoonk (2024)

#References 

#Vermeer, A. (2021). Ecological stability in the face of climatic disturbances. MSc Thesis. Utrecht University.

# Dakos, V., Carpenter, S. R., Brock, W. A., Ellison, A. M., Guttal, V., Ives, A. R., Kéfi, S., Livina, 
#V., Seekell, D. A., van Nes, E. H., & Scheffer, M. (2012). Methods for detecting early warnings of critical 
#transitions in time series illustrated using simulated ecological data. PLoS ONE, 7(7). https://doi.org/10.1371/journal.pone.0041010

###################################
## DATA SELECTION               ###
###################################

## INPUT
#selecting the entire NDVI raster brick, catchment outline, channel outlines and plantation outlines from Vermeer (2021)
catchment_outline <- terra::vect("your_file_location/Ounila_catchment_reprojected.shp")
channels <- terra::rast("your_file_location/Rchannels_buffer100m.grd")
plantations <- terra::rast("your_file_location/Rplantations.grd")
ras_brick_full <- terra::rast("your_file_location/NDVI_final_OLS.grd")

##PARAMETERS
starttime <- 1 #this represents the start of the selected timeseries (in image number).
endtime <- 570 #the end of the selected timeseries (in image number). 
# Image 800 represents the last image for 2002 (end drought), 570 is last image of 1998 (start drought). See Vermeer (2021) for more info
# Increasing the length of the used timeseries increases RAM/memory consumption! Take your memory availability into account. 

#selecting desired layers corresponding with the starttime and endtime 
ras_brick_selection <- subset(ras_brick_full, starttime:endtime) 

#Masking out irrigated lands in channels and plantations 
irrigated <- mosaic(channels, plantations)
irr_mask <- ifel(is.na(irrigated), 1, NA)
irr_mask_proj <- project(irr_mask, ras_brick_selection) #matching projection of the mask to the ras_brick_selection
NDVI_data <- mask(ras_brick_selection, irr_mask_proj)

###################################
## DEFINING FUNCTIONS           ###
###################################
#These functions are adapted versions of the "surrogates_EWS" function from the earlywarnings package (Dakos et al., 2012). 
#They have been adjusted for application using the ''app'' function from the terra package. 
#Some simplifications were made to data input and output. However, the function should still support most original operations.

compute_variance <- function (i, indicator = "sd", winsize = 50, detrending = "loess", 
                         bandwidth = NULL, span = NULL, degree = NULL, boots = 100, 
                         logtransform = FALSE, interpolate = TRUE, maxNA = 114) #MaxNA represents the max number of NA's allowed. 
{
  skewness <- moments::skewness
  kurtosis <- moments::kurtosis
  Y <- i
  if ((sum(is.na(Y)) > maxNA)){
    output <- NA
  }
  else  {
  
  timeindex <- 1:length(Y)

  if (interpolate) {
    YY <- approx(timeindex, Y, n = length(Y), method = "linear")
    Y <- YY$y
  }
  else {
    Y <- Y
  }
  if (logtransform) {
    Y <- log(Y + 1)
  }
  detrending <- match.arg(detrending)
  if (detrending == "gaussian") {
    if (is.null(bandwidth)) {
      bw <- round(bw.nrd0(timeindex))
    }
    else {
      bw <- round(length(Y) * bandwidth)/100
    }
    smYY <- ksmooth(timeindex, Y, kernel = c("normal"), 
                    bandwidth = bw, range.x = range(timeindex), n.points = length(timeindex))
    nsmY <- Y - smYY$y
    smY <- smYY$y
  }
  else if (detrending == "loess") {
    if (is.null(span)) {
      span <- 25/100
    }
    else {
      span <- span/100
    }
    if (is.null(degree)) {
      degree <- 2
    }
    else {
      degree <- degree
    }
    smYY <- loess(Y ~ timeindex, span = span, degree = degree, 
                  normalize = FALSE, family = "gaussian")
    smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
    nsmY <- Y - smY
  }
  else if (detrending == "linear") {
    nsmY <- resid(lm(Y ~ timeindex))
    smY <- fitted(lm(Y ~ timeindex))
  }
  else if (detrending == "first-diff") {
    nsmY <- diff(Y)
    timeindexdiff <- timeindex[1:(length(timeindex) - 1)]
  }
  else if (detrending == "no") {
    smY <- Y
    nsmY <- Y
  }
  mw <- round(length(Y) * winsize)/100
  omw <- length(nsmY) - mw + 1
  low <- 6
  high <- omw
  nMR <- matrix(data = NA, nrow = mw, ncol = omw)
  for (i in 1:omw) {
    Ytw <- nsmY[i:(i + mw - 1)]
    nMR[, i] <- Ytw
  }
  indicator = match.arg(indicator)
  if (indicator == "ar1") {
    indic <- apply(nMR, 2, function(x) {
      nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, 
                     intercept = FALSE)
      nAR1$ar
    })
  }
  else if (indicator == "sd") {
    indic <- apply(nMR, 2, sd)
  }
  else if (indicator == "sk") {
    indic <- apply(nMR, 2, skewness)
  }
  else if (indicator == "kurt") {
    indic <- apply(nMR, 2, kurtosis)
  }
  else if (indicator == "acf1") {
    indic <- apply(nMR, 2, function(x) {
      nACF <- acf(x, lag.max = 1, type = c("correlation"), 
                  plot = FALSE)
      nACF$acf[2]
    })
  }
  else if (indicator == "returnrate") {
    indic <- apply(nMR, 2, function(x) {
      nACF <- acf(x, lag.max = 1, type = c("correlation"), 
                  plot = FALSE)
      1 - nACF$acf[2]
    })
  }
  else if (indicator == "cv") {
    indic <- apply(nMR, 2, function(x) {
      sd(x)/mean(x)
    })
  }
  else if (indicator == "densratio") {
    indic <- apply(nMR, 2, function(x) {
      spectfft <- spec.ar(x, n.freq = omw, plot = FALSE, 
                          order = 1)
      spectfft$spec
      spectfft$spec[low]/spectfft$spec[high]
    })
  }
  timevec <- seq(1, length(indic))
  Kt <- cor.test(timevec, indic, alternative = c("two.sided"), 
                 method = c("kendall"), conf.level = 0.95)
  Ktauestindorig <- Kt$estimate
  output <- as.numeric(Ktauestindorig)
  
  }
  return(output)
}

compute_AC1 <- function (i, indicator = "acf1", winsize = 50, 
                              detrending = "loess", 
                              bandwidth = NULL, span = NULL, degree = NULL, boots = 100, 
                              logtransform = FALSE, interpolate = TRUE, maxNA = 114) 
{
  skewness <- moments::skewness
  kurtosis <- moments::kurtosis
  Y <- i
  if ((sum(is.na(Y)) > maxNA)){
    output <- NA
  }
  else  {
    
    timeindex <- 1:length(Y)
    
    if (interpolate) {
      YY <- approx(timeindex, Y, n = length(Y), method = "linear")
      Y <- YY$y
    }
    else {
      Y <- Y
    }
    if (logtransform) {
      Y <- log(Y + 1)
    }
    detrending <- match.arg(detrending)
    if (detrending == "gaussian") {
      if (is.null(bandwidth)) {
        bw <- round(bw.nrd0(timeindex))
      }
      else {
        bw <- round(length(Y) * bandwidth)/100
      }
      smYY <- ksmooth(timeindex, Y, kernel = c("normal"), 
                      bandwidth = bw, range.x = range(timeindex), n.points = length(timeindex))
      nsmY <- Y - smYY$y
      smY <- smYY$y
    }
    else if (detrending == "loess") {
      if (is.null(span)) {
        span <- 25/100
      }
      else {
        span <- span/100
      }
      if (is.null(degree)) {
        degree <- 2
      }
      else {
        degree <- degree
      }
      smYY <- loess(Y ~ timeindex, span = span, degree = degree, 
                    normalize = FALSE, family = "gaussian")
      smY <- predict(smYY, data.frame(x = timeindex), se = FALSE)
      nsmY <- Y - smY
    }
    else if (detrending == "linear") {
      nsmY <- resid(lm(Y ~ timeindex))
      smY <- fitted(lm(Y ~ timeindex))
    }
    else if (detrending == "first-diff") {
      nsmY <- diff(Y)
      timeindexdiff <- timeindex[1:(length(timeindex) - 1)]
    }
    else if (detrending == "no") {
      smY <- Y
      nsmY <- Y
    }
    mw <- round(length(Y) * winsize)/100
    omw <- length(nsmY) - mw + 1
    low <- 6
    high <- omw
    nMR <- matrix(data = NA, nrow = mw, ncol = omw)
    for (i in 1:omw) {
      Ytw <- nsmY[i:(i + mw - 1)]
      nMR[, i] <- Ytw
    }
    indicator = match.arg(indicator)
    if (indicator == "ar1") {
      indic <- apply(nMR, 2, function(x) {
        nAR1 <- ar.ols(x, aic = FALSE, order.max = 1, dmean = FALSE, 
                       intercept = FALSE)
        nAR1$ar
      })
    }
    else if (indicator == "sd") {
      indic <- apply(nMR, 2, sd)
    }
    else if (indicator == "sk") {
      indic <- apply(nMR, 2, skewness)
    }
    else if (indicator == "kurt") {
      indic <- apply(nMR, 2, kurtosis)
    }
    else if (indicator == "acf1") {
      indic <- apply(nMR, 2, function(x) {
        nACF <- acf(x, lag.max = 1, type = c("correlation"), 
                    plot = FALSE)
        nACF$acf[2]
      })
    }
    else if (indicator == "returnrate") {
      indic <- apply(nMR, 2, function(x) {
        nACF <- acf(x, lag.max = 1, type = c("correlation"), 
                    plot = FALSE)
        1 - nACF$acf[2]
      })
    }
    else if (indicator == "cv") {
      indic <- apply(nMR, 2, function(x) {
        sd(x)/mean(x)
      })
    }
    else if (indicator == "densratio") {
      indic <- apply(nMR, 2, function(x) {
        spectfft <- spec.ar(x, n.freq = omw, plot = FALSE, 
                            order = 1)
        spectfft$spec
        spectfft$spec[low]/spectfft$spec[high]
      })
    }
    timevec <- seq(1, length(indic))
    Kt <- cor.test(timevec, indic, alternative = c("two.sided"), 
                   method = c("kendall"), conf.level = 0.95)
    Ktauestindorig <- Kt$estimate
    output <- as.numeric(Ktauestindorig)
  }
  return(output)
}

#A simple function to calculate Kendall's tau for greening/browning trends in the area 

compute_greening <- function (i, maxNA = 114) 
{
  Y <- i
  if ((sum(is.na(Y)) > maxNA)){
    output <- NA
  }
  else  {
  timeindex <- 1:length(Y)
  kt <- cor.test(timeindex, Y, alternative = "two.sided", method = "kendall")
  Kttest <- kt$estimate
  output <- as.numeric(Kttest)
  }
  return(output)
}

###################################
## APPLYING FUNCTIONS & sAVING  ###
###################################
#Note that these are set up for parallel processing across 7 cores. Adjust this for the available processing power in your system. 
#Parallel processing can be disabled by just removing the argument but processing will take a long time. 

greening_output <- app(NDVI_data, compute_greening, cores=7)
f <- file.path(getwd(), "greening_output.grd")
writeRaster(greening_output, f)

######################################
## CHECKING CLIMATE VARIABILITY    ###
######################################

#temporarily changing the wd for import
setwd("your_file_location/precipitation data")
rastlist <- list.files(path = "your_file_location/precipitation data", pattern='.tif$', 
                       all.files= T, full.names= T)
precipitation <- terra::rast(list)

#changing back the wd, applying the compute_variance function and saving for different WS. 
setwd("your_file_location")
precip_tau_output <- app(precipitation, compute_variance, cores=7)
f <- file.path(getwd(), "precip_tau_output_WS50.grd")
writeRaster(precip_tau_output, f)

######################################
## PLOTTING MAPS FOR INITIAL CHECK ###
###################################### 
#This can be used to quickly check the current output for variance and AC1. For further figures and maps, see data quality and figures scripts. 

mapcolours <- hcl.colors(100, palette = "Temps")

plot(variance_tau_output,col= mapcolours, main = "Kendall's Tau Variance", range=c(-1,1))
lines(catchment_outline)

plot(AC1_tau_output,col= mapcolours, main = "Kendall's Tau AC1", range=c(-1,1))
lines(catchment_outline)




