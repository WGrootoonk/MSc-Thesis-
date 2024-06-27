# Activating required packages
library(terra)
library(earlywarnings)
library(ggplot2)
library(ggpubr)
library(reshape2)
## Setting working directory. 
#Replace ''your_file_location" with your data folder. Using ctrl+f & replace is recommended for the rest of the script. 
setwd("your_file_location")

### DATA QUALITY CHECK SCRIPT ###
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
#selecting the entire NDVI raster brick, catchment outline and channel outlines from Vermeer (2021)
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
## DATA QUALITY CHECK          ###
###################################

##checking for amount of NA values 
# A quick function to check the amount of NA's for each pixel and turn this into a raster
NAcounter <- function(i){
  if(sum(is.na(i))==length(i)){ #This is done so 100% NA pixels (i.e. out of bounds pixels) are excluded
    out <- NA }
  else{
    out <-(sum(is.na(i))/(length(i))*100)}
  return(out)
} 

NAcount <- app(ras_brick_selection, NAcounter)
plot(NAcount, col= (hcl.colors(100, palette="Temps")), main = "Percentage missing values for each pixel", legend = TRUE)
lines(catchment_outline)


# A function to check if the amount of NA's per pixel changes through time and turns this into a raster

pxNAtrend <- function(i){
  if(sum(is.na(i))==length(i)){ #This is done so 100% NA pixels (i.e. out of bounds pixels) are excluded
    out <- NA }
  else{
    timeindex <- 1:length(i)
    num_na <- as.numeric(is.na(i))
    cortest <- cor.test(timeindex, num_na, alternative = "two.sided", method= "kendall")
    out <- as.numeric(cortest$estimate)
  }
  return(out)
}

pxNAtrend_output <- app(NDVI_data, pxNAtrend, cores=7)
mapcolours <- hcl.colors(100, palette = "Temps")
plot(pxNAtrend_output,col=mapcolours, main = "Kendall's Tau NA values")

#Alternatively the trend in NA frequency can also be computed on a layer by layer basis. 
#Note that this will included NA values making up out of bounds areas, which are excluded by earlier functions. 
NAtrend <- data.frame(layers = (1:nlyr(NDVI_data)), percentage = NA)

for(i in 1:nlyr(NDVI_data)){
  layer <- NDVI_data[[i]]
  NApercent <- (sum(values(allNA(layer)))/ncell(layer))
  NAtrend$percentage[i] <- (NApercent*100)
  print(i)
}

ggscatter(NAtrend, x = "layers", y = "percentage", size = 1, 
          add = "reg.line", cor.coef = TRUE, cor.method = "pearson", 
          title = "Pearson correlation NA values", xlab = "Layer number", ylab = "NA values (in % of total cells)")

############################################
## NULL MODEL & CORRECTION RASTERs       ###
############################################

#Null Script. This uses the same function as in the main resilience indicator script, except that the non-NA values of the timeseries are replaced with 
# a random value between 0 and 1. If NA's have no effect, a noisy normal distribution across the catchment is expected (given the input). Any deviations
# represent bias introduced by missing values. 


tau_null <- function (i, indicator = "acf1", winsize = 50, detrending = "loess", 
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
    Y[!is.na(Y)] <- runif((sum(!is.na(Y))), min=0, max=1)
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

#Running the above several times for variance and AC1 at different window sizes and saving. 
#Note that these are set up for parallel processing across 7 cores. Adjust this for the available processing power in your system. 
#Parallel processing can be disabled by just removing the argument but processing will take a long time. 

tau_null <- app(NDVI_data, tau_null, cores=7)
f <- file.path(getwd(), "outputname.grd")
writeRaster(tau_null, f)

####################################
## ANALYZING CORRECTION RASTERs  ###
####################################

#Importing files produced using null model as described above. 
var_null_WS20 <- terra::rast("your_file_location/var_null_WS20.grd")
var_null_WS30 <- terra::rast("your_file_location/var_null_WS30.grd")
var_null_WS40 <- terra::rast("your_file_location/var_null_WS40.grd")
var_null_WS50 <- terra::rast("your_file_location/var_null_WS50.grd")
var_null_mean <- mean(c(var_null_WS20, var_null_WS30, var_null_WS40, var_null_WS50))

#Saving mean result
f <- file.path(getwd(), "var_null_mean.grd")
writeRaster(var_null_mean, f)

AC1_null_WS20 <- terra::rast("your_file_location/AC1_null_WS20.grd")
AC1_null_WS30 <- terra::rast("your_file_location/AC1_null_WS30.grd")
AC1_null_WS40 <- terra::rast("your_file_location/AC1_null_WS40.grd")
AC1_null_WS50 <- terra::rast("your_file_location/AC1_null_WS50.grd")
AC1_null_mean <- mean(c(AC1_null_WS20, AC1_null_WS30, AC1_null_WS40, AC1_null_WS50))

#Saving mean result
f <- file.path(getwd(), "AC1_null_mean.grd")
writeRaster(AC1_null_mean, f)

var_null_mean <- terra::rast("your_file_location/var_null_mean.grd")
AC1_null_mean <- terra::rast("your_file_location/AC1_null_mean.grd")

#The restulting data is very noisy. Here we downscale the raster by a factor 5 by taking the mean and then resize it to the original
# resolution, leading to smoother values 
var_null_mean_agg <- terra::aggregate(var_null_mean, fact=5, fun='mean')
var_correction_raster <- project((terra::disagg(var_null_mean_agg, fact=5, method="bilinear")), var_null_mean)

AC1_null_mean_agg <- terra::aggregate(AC1_null_mean, fact=5, fun='mean')
AC1_correction_raster <- project((terra::disagg(AC1_null_mean_agg, fact=5, method="bilinear")), var_null_mean)


#Saving the correction rasters for later use
f <- file.path(getwd(), "var_correction_raster.grd")
writeRaster(var_correction_raster, f, overwrite=TRUE)

f <- file.path(getwd(), "AC1_correction_raster.grd")
writeRaster(AC1_correction_raster, f, overwrite=TRUE)


#Plotting the correction rasters to check what they look like
mapcolours <- hcl.colors(100, palette = "Temps")
plot(var_correction_raster, col= mapcolours, main = "Correction raster variance", legend = TRUE, range= c(-1,1))
lines(catchment_outline)
plot(AC1_correction_raster, col= mapcolours, main = "Correction raster AC1", legend = TRUE, range= c(-1,1))
lines(catchment_outline)


# Making some additional figures for the correction rasters 
var_cor_values <- as.numeric(values(var_correction_raster))
AC1_cor_values <- as.numeric(values(AC1_correction_raster))

data = data.frame(Variance = var_cor_values, AC1 = AC1_cor_values)

ggplot(data, aes(x=Variance))+
  xlim(-1,1)+
  theme(axis.text.y = element_blank())+
  geom_boxplot(color = "black", fill="darkgrey")+
  labs(x="Kendall's Tau (variance corr.)")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)

ggplot(data, aes(x=AC1))+
  xlim(-1,1)+
  theme(axis.text.y = element_blank())+
  geom_boxplot(color = "black", fill="darkgrey")+
  labs(x="Kendall's Tau (AC1 corr.)")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)

ggarrange(var_cor_box, AC1_cor_box, legend=FALSE)


# Plotting overlaying density graphs with a mean line to check for value distribution
x <- data.frame(Variance=var_null_values, AC1=AC1_null_values)
data <- melt(x)
ggplot(data,aes(x=value, fill=variable)) + geom_density(alpha=0.25, na.rm = TRUE) + 
  labs(title = "Value density null models", x="Kendall's Tau", y="Density")+ 
  geom_vline(aes(xintercept=mean(var_null_values, na.rm = TRUE)), color="red", linetype="dashed", linewidth=1)+
  geom_vline(aes(xintercept=mean(AC1_null_values, na.rm=TRUE)), color="blue", linetype="dashed", linewidth=1)










