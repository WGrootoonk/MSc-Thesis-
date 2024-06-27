library(terra)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggridges)
library(viridis)
library(ggcorrplot)

## Setting working directory. 
#Replace ''your_file_location" with your data folder. Using ctrl+f & replace is recommended for the rest of the script. 
setwd("your_file_location") 

### FIGURES SCRIPT ###
# Developed for MSc Thesis ''Relations between temporal resilience indicators and trend breakpoints in a dryland high-mountain catchment"
# Developed by Willem Grootoonk (2024)

#References 

#Vermeer, A. (2021). Ecological stability in the face of climatic disturbances. MSc Thesis. Utrecht University.

# Dakos, V., Carpenter, S. R., Brock, W. A., Ellison, A. M., Guttal, V., Ives, A. R., Kéfi, S., Livina, 
#V., Seekell, D. A., van Nes, E. H., & Scheffer, M. (2012). Methods for detecting early warnings of critical 
#transitions in time series illustrated using simulated ecological data. PLoS ONE, 7(7). https://doi.org/10.1371/journal.pone.0041010


#########################################
## MAPPING RESILLIENCE INDICATORS     ###
#########################################
##Importing correction for NA values from null model generated in the data quality script 
var_correction_raster <- terra::rast("your_file_location/var_correction_raster.grd")
AC1_correction_raster <- terra::rast("your_file_location/AC1_correction_raster.grd")

#Importing outputs for various windowsizes, applying correction. 
variance_output_WS20 <- (terra::rast("your_file_location/variance_output_WS20.grd")-
                                    var_correction_raster)
variance_output_WS30 <- (terra::rast("your_file_location/variance_output_WS30.grd")-
                                    var_correction_raster)
variance_output_WS40 <- (terra::rast("your_file_location/variance_output_WS40.grd")-
                                    var_correction_raster)
variance_output_WS50 <- (terra::rast("your_file_location/variance_output_WS50.grd")-
                                    var_correction_raster)
variance_output_mean <- mean(c(variance_output_WS20,variance_output_WS30,variance_output_WS40,variance_output_WS50))


AC1_output_WS20 <- (terra::rast("your_file_location/AC1_output_WS20.grd")- 
                               AC1_correction_raster)
AC1_output_WS30 <- (terra::rast("your_file_location/AC1_output_WS30.grd")- 
                               AC1_correction_raster)
AC1_output_WS40 <- (terra::rast("your_file_location/AC1_output_WS40.grd")- 
                               AC1_correction_raster)
AC1_output_WS50 <- (terra::rast("your_file_location/AC1_output_WS50.grd")- 
                               AC1_correction_raster)
AC1_output_mean <- mean(c(AC1_output_WS20,AC1_output_WS30,AC1_output_WS40,AC1_output_WS50))

#Producing maps for each WS and the means
mapcolours <- hcl.colors(100, palette = "Temps")
catchment_outline <- terra::vect("your_file_location/Ounila_catchment_reprojected.shp")

plot(variance_output_WS20, col= mapcolours, main = "WS=20", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(variance_output_WS30, col= mapcolours, main = "WS=30", range =c(-1,1),legend = FALSE)
lines(catchment_outline)
plot(variance_output_WS40, col= mapcolours, main = "WS=40", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(variance_output_WS50, col= mapcolours, main = "WS=50", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(variance_output_mean, col= mapcolours, main = "Mean Kendall's Tau Variance", range =c(-1,1), legend = TRUE)
lines(catchment_outline)

plot(AC1_output_WS20, col= mapcolours, main = "WS=20", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(AC1_output_WS30, col= mapcolours, main = "WS=30", range =c(-1,1),legend = FALSE)
lines(catchment_outline)
plot(AC1_output_WS40, col= mapcolours, main = "WS=40", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(AC1_output_WS50, col= mapcolours, main = "WS=50", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(AC1_output_mean, col= mapcolours, main = "Mean Kendall's Tau AC1", range =c(-1,1), legend = TRUE)
lines(catchment_outline)


#Importing rasters for different cut off dates & applying correction raster
variance_output_1998_01_04 <- (terra::rast("your_file_location/variance_output_19980104.grd")
-var_correction_raster)                                     
variance_output_1997_01_10 <- (terra::rast("your_file_location/variance_output_19970110.grd")
-var_correction_raster) 
variance_output_1996_01_08 <- (terra::rast("your_file_location/variance_output_19960108.grd")
-var_correction_raster) 

AC1_output_1998_01_04 <- (terra::rast("your_file_location/AC1_output_19980104.grd")
-AC1_correction_raster)                                     
AC1_output_1997_01_10 <- (terra::rast("your_file_location/AC1_output_19970110.grd")
-AC1_correction_raster) 
AC1_output_1996_01_08 <- (terra::rast("your_file_location/AC1_output_19960108.grd")
-AC1_correction_raster) 

#plotting these figures
plot(variance_output_1998_01_04, col= mapcolours, main = "Until 04/01/1998", range =c(-1,1), legend = T)
lines(catchment_outline)
plot(variance_output_1997_01_10, col= mapcolours, main = "Until 10/01/1997", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(variance_output_1996_01_08, col= mapcolours, main = "Until 08/01/1996", range =c(-1,1), legend = FALSE)
lines(catchment_outline)

plot(AC1_output_1998_01_04, col= mapcolours, main = "Until 04/01/1998", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(AC1_output_1997_01_10, col= mapcolours, main = "Until 10/01/1997", range =c(-1,1), legend = FALSE)
lines(catchment_outline)
plot(AC1_output_1996_01_08, col= mapcolours, main = "Until 08/01/1996", range =c(-1,1), legend = FALSE)
lines(catchment_outline)

#########################################
## COMPARING BETWEEN BREAKPOINT TYPES ###
#########################################
#importing breakpoint type raster from Vermeer (2021)
DBPType_raster <- terra::rast("your_file_location/DBPType.grd")

##creating masks for dominant responses & 
#matching the projections of the rasters so this doesn't cause issues later. The variance raster is used but AC1 is identical. 

# DBPtype = 0 No significant trend (NT)
# DBPtype = 1 interrupted increase (II)
# DBPtype = 2 increase after no trend (IN)
# DBPtype = 3 No trend after decrease (ND)
# DBPtype = 4 positive reversal (PR)
# DBPtype = 5 interrupted decrease (ID)
# DBPtype = 6 decrease after no trend (DN)
# DBPtype = 7 No trend after in increase (NI)
# DBPtype = 8 Negative reversal (NR)

#Defining names and colours corresponding with Vermeer (2021)
BPcolours <- c("grey85",rev(paste0(viridis(n=5,option="D"))[1:4]),
               rev(paste0(viridis(n=6,option="plasma"))[3:6]))
BPnames <- c("No significant trend (NT)",
             "Interrupted increase (II) ᵖ","Increase after no trend (IN) ᵖ","No trend after decrease (ND) ᵖ","Positive reversal (PR) ᵖ",
             "Interrupted decrease (ID) ⁿ","Decrease after no trend (DN) ⁿ","No trend after increase (NI) ⁿ","Negative reversal (NR) ⁿ")

#Creating masks for each breakpoint type, matching projection with output. 
NT_mask <- project((ifel(DBPType_raster == 0, 1, NA)), variance_output_mean)
II_mask <- project((ifel(DBPType_raster == 1, 1, NA)), variance_output_mean)
IN_mask <- project((ifel(DBPType_raster == 2, 1, NA)), variance_output_mean)
ND_mask <- project((ifel(DBPType_raster == 3, 1, NA)), variance_output_mean)
PR_mask <- project((ifel(DBPType_raster == 4, 1, NA)), variance_output_mean)
ID_mask <- project((ifel(DBPType_raster == 5, 1, NA)), variance_output_mean)
DN_mask <- project((ifel(DBPType_raster == 6, 1, NA)), variance_output_mean)
NI_mask <- project((ifel(DBPType_raster == 7, 1, NA)), variance_output_mean)
NR_mask <- project((ifel(DBPType_raster == 8, 1, NA)), variance_output_mean)
PosBP_mask <- merge(II_mask, IN_mask, PR_mask, ND_mask)
NegBP_mask <- merge(ID_mask, DN_mask, NI_mask, NR_mask)

#Variance masking
var_NT_mask <- mask(variance_output_mean, NT_mask)
var_II_mask <- mask(variance_output_mean, II_mask)
var_IN_mask <- mask(variance_output_mean, IN_mask)
var_ND_mask <- mask(variance_output_mean, ND_mask)
var_PR_mask <- mask(variance_output_mean, PR_mask)
var_ID_mask <- mask(variance_output_mean, ID_mask)
var_DN_mask <- mask(variance_output_mean, DN_mask)
var_NI_mask <- mask(variance_output_mean, NI_mask)
var_NR_mask <- mask(variance_output_mean, NR_mask)
var_PosBP_mask <- mask(variance_output_mean, PosBP_mask)
var_NegBP_mask <- mask(variance_output_mean, NegBP_mask)

##Autocorrelation masking
AC1_NT_mask <- mask(AC1_output_mean, NT_mask)
AC1_II_mask <- mask(AC1_output_mean, II_mask)
AC1_IN_mask <- mask(AC1_output_mean, IN_mask)
AC1_ND_mask <- mask(AC1_output_mean, ND_mask)
AC1_PR_mask <- mask(AC1_output_mean, PR_mask)
AC1_ID_mask <- mask(AC1_output_mean, ID_mask)
AC1_DN_mask <- mask(AC1_output_mean, DN_mask)
AC1_NI_mask <- mask(AC1_output_mean, NI_mask)
AC1_NR_mask <- mask(AC1_output_mean, NR_mask)
AC1_PosBP_mask <- mask(AC1_output_mean, PosBP_mask)
AC1_NegBP_mask <- mask(AC1_output_mean, NegBP_mask)

# Extracting values as a numeric for use in further figures
var_NT <- as.numeric(values(var_NT_mask))
var_II <- as.numeric(values(var_II_mask))
var_IN <- as.numeric(values(var_IN_mask))
var_ND <- as.numeric(values(var_ND_mask))
var_PR <- as.numeric(values(var_PR_mask))
var_ID <- as.numeric(values(var_ID_mask))
var_DN <- as.numeric(values(var_DN_mask))
var_NI <- as.numeric(values(var_NI_mask))
var_NR <- as.numeric(values(var_NR_mask))
var_PosBP <- as.numeric(values(var_PosBP_mask))
var_NegBP <-as.numeric(values(var_NegBP_mask))

AC1_NT <- as.numeric(values(AC1_NT_mask))
AC1_II <- as.numeric(values(AC1_II_mask))
AC1_IN <- as.numeric(values(AC1_IN_mask))
AC1_ND <- as.numeric(values(AC1_ND_mask))
AC1_PR <- as.numeric(values(AC1_PR_mask))
AC1_ID <- as.numeric(values(AC1_ID_mask))
AC1_DN <- as.numeric(values(AC1_DN_mask))
AC1_NI <- as.numeric(values(AC1_NI_mask))
AC1_NR <- as.numeric(values(AC1_NR_mask))
AC1_PosBP <- as.numeric(values(AC1_PosBP_mask))
AC1_NegBP <-as.numeric(values(AC1_NegBP_mask))


########################################
## Making plots for breakpoint types ###
########################################

##Creating density plots for postive and negative breakpoints
#collecting negative and postive values in a dataframe. 
vardata <- data.frame(P = var_PosBP, N = var_NegBP)
vardata_m <- melt(vardata, variable.name = 'Breakpoint', na.rm=TRUE)
AC1data <- data.frame(P = AC1_PosBP, N = AC1_NegBP)
AC1data_m <- melt(AC1data, variable.name = 'Breakpoint', na.rm=TRUE)
data <- data.frame(Breakpoint = vardata_m$Breakpoint, Variance = vardata_m$value, AC1 = AC1data_m$value)
#I am sure there is a better way to do the above but this works

vardens <- ggplot(data,aes(x=Variance, fill=Breakpoint)) + geom_density(alpha=0.25, na.rm = TRUE) + 
  labs(title = "Value density variance", x="Kendall's Tau", y="Density")+ xlim(-1,1)+
  #geom_vline(aes(xintercept=mean(var_PosBP, na.rm =TRUE)), color= "blue", linetype="dashed", linewidth=1)+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)+
  #geom_vline(aes(xintercept=mean(var_NegBP, na.rm=TRUE)), color= "red", linetype="dashed", linewidth=1)+
  scale_fill_manual(values = c("blue","red"), labels= c("Postive Breakpoints", "Negative Breakpoints"),name=NULL)

AC1dens <- ggplot(data,aes(x=AC1, fill=Breakpoint)) + geom_density(alpha=0.25, na.rm = TRUE) + 
  labs(title = "Value density AC1", x="Kendall's Tau", y="Density")+ xlim(-1,1)+
  #geom_vline(aes(xintercept=mean(AC1_PosBP, na.rm = TRUE)), color= "blue", linetype="dashed", linewidth=1)+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)+
  #geom_vline(aes(xintercept=mean(AC1_NegBP, na.rm=TRUE)), color= "red", linetype="dashed", linewidth=1)+
  scale_fill_manual(values = c("blue","red"), labels= c("Postive Breakpoints", "Negative Breakpoints"), name=NULL)

varbox <- ggplot(data, aes(x=Variance, y=Breakpoint, fill=Breakpoint))+
  xlim(-1,1)+
  geom_boxplot(alpha=0.25, na.rm=TRUE)+
  labs(title = "Boxplot variance", x="Kendall's Tau", y="Breakpoint")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)+
  scale_fill_manual(values = c("blue","red"), labels= c("Postive Breakpoints", "Negative Breakpoints"),name=NULL)

AC1box <- ggplot(data, aes(x=AC1, fill = Breakpoint, y=Breakpoint))+
  xlim(-1,1)+
  geom_boxplot(alpha=0.25, na.rm=TRUE)+
  labs(title = "Boxplot AC1", x="Kendall's Tau", y="Breakpoint")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)+
  scale_fill_manual(values = c("blue","red"), labels= c("Postive Breakpoints", "Negative Breakpoints"),name=NULL)

ggarrange(vardens, AC1dens, varbox, AC1box, ncol=2, nrow=2, common.legend = TRUE, labels = c("A","B","C","D"), legend = 'bottom')

## Making density plots for all breakpoint types
#collecting  values in a dataframe. 
vardata <- data.frame(NT=var_NT, II=var_II, IN=var_IN, ND=var_ND, PR=var_PR, ID=var_ID, DN=var_DN, NI=var_NI, NR=var_NR)
AC1data <- data.frame(NT=AC1_NT, II=AC1_II, IN=AC1_IN, ND =AC1_ND, PR=AC1_PR, ID=AC1_ID, DN=AC1_DN, NI=AC1_NI, NR=AC1_NR)
vardata_m <- melt(vardata, variable.name = 'Breakpoint', na.rm=TRUE)
AC1data_m <- melt(AC1data, variable.name = 'Breakpoint', na.rm=TRUE)
data <- data.frame(Breakpoint = vardata_m$Breakpoint, Variance = vardata_m$value, AC1 = AC1data_m$value)

vardens_all <- ggplot(data,aes(x=Variance, y=Breakpoint, fill = Breakpoint)) + geom_density_ridges()+
  xlim(-1,1)+
  scale_fill_manual(values = BPcolours, labels=BPnames, name="Breakpoint type")+
  labs(title = "Value density variance", x="Kendall's Tau", y="Breakpoint")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)

AC1dens_all <- ggplot(data,aes(x=AC1, y=Breakpoint, fill = Breakpoint)) + geom_density_ridges()+
  xlim(-1,1)+
  scale_fill_manual(values = BPcolours, labels=BPnames, name="Breakpoint type")+
  labs(title = "Value density AC1", x="Kendall's Tau", y="Breakpoint")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)

var_box <- ggplot(data, aes(x=Variance, y=Breakpoint, fill = Breakpoint))+
xlim(-1,1)+
scale_fill_manual(values = BPcolours, labels=BPnames, name="Breakpoint type")+
geom_boxplot()+
labs(title = "Boxplot variance", x="Kendall's Tau", y="Breakpoints")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)

AC1_box <- ggplot(data, aes(x=AC1, y=Breakpoint, fill = Breakpoint))+
  xlim(-1,1)+
  scale_fill_manual(values = BPcolours, labels=BPnames, name="Breakpoint type")+
  geom_boxplot()+
labs(title = "Boxplot AC1", x="Kendall's Tau", y="Breakpoints")+
  geom_vline(xintercept=0.0, color="black", linetype="dashed", linewidth=0.5)

ggarrange(vardens_all, AC1dens_all, var_box, AC1_box, ncol=2, nrow=2, common.legend = TRUE, 
          labels = c("A","B","C","D"), legend = 'right')


#########################################
## COMPARING BETWEEN OTHER FACTORS    ###
#########################################
# Importing precipitation variability output from the functions script. 
precip_output_WS20 <- terra::rast("your_file_location/precip_tau_output_WS20.grd")
precip_output_WS30 <- terra::rast("your_file_location/precip_tau_output_WS30.grd")
precip_output_WS40 <- terra::rast("your_file_location/precip_tau_output_WS40.grd")
precip_output_WS50 <- terra::rast("your_file_location/precip_tau_output_WS50.grd")
precip_output_mean <- mean(c(precip_output_WS20,precip_output_WS30,precip_output_WS40,precip_output_WS50))
precip_output_mean_res <- resample(precip_output_mean, AC1_output_mean)

#Importing data for slope, elevation, salinity, greening, 1998 NDVI and the earlier precip. variability. 
slope <- terra::rast("your_file_location/slope.grd")
DEM_res<- resample((terra::rast("your_file_location/DEM.grd")), AC1_output_mean)
DEM <- project(DEM_res, AC1_output_mean) #For some reason the DEM has a different resolution and CRS as the rest, corrected here.
salinity <- project(terra::rast("your_file_location/NSImap.tif"), AC1_output_mean)
greening <- terra::rast("your_file_location/greening_output.grd")
NDVI_1998 <- terra::rast("your_file_location/NDVI_1998_mean.grd")
precip <- project(crop(precip_output_mean_res, catchment_outline), AC1_output_mean)

#This can be used to select the elevation range you are interested in. Default = 1200 to 3600 (entire catchment)
# This range can roughly be divided in 4 quartiles, see: quantile((as.numeric(values(DEM))), probs = c(0,0.25,0.5,0.75,1), na.rm=TRUE) 
# Q1 (Low)  = 1200 to 1750
# Q2 (Low-intermediate)  = 1750 to 1950
# Q3 (Intermediate) = 1950 to 2200
# Q4 (High)  = 2200 to 3600 
minele <- ifel(DEM > 1200, 1, NA) 
maxele <- ifel(DEM < 3600, 1, NA)
elemask <- maxele*minele

#Masking data to desired elevation range
slope_masked <- mask(slope, elemask)
salinity_masked <- mask(salinity, elemask)
DEM_masked <- mask(DEM, elemask)
greening_masked <- mask(greening, elemask)
NDVI_1998_masked <- mask(NDVI_1998, elemask)
precip_masked <- mask(precip, elemask)
#Turning into numeric values
slopevalues <- as.numeric(values(slope_masked))
DEMvalues <- as.numeric(values(DEM_masked))
sal_values <- as.numeric(values(salinity_masked))
greeningvalues <- as.numeric(values(greening_masked))
NDVI_1998values <- as.numeric(values(NDVI_1998_masked))
precip_values <- as.numeric(values(precip_masked))
#Repeating variance and AC1 values from earlier for convenience
variance_tau_values <- as.numeric(values(variance_output_mean))
AC1_tau_values <- as.numeric(values(AC1_output_mean))

#Creating a correlation matrix
data <- data.frame(Salinity = sal_values, Slope = slopevalues, Mean_NDVI_1998 = NDVI_1998values, "Greening trend" = greeningvalues, 
                   Elevation = DEMvalues, precip_var_trend = precip_values, "Variance Kendall's tau" = variance_tau_values, "AC1 Kendall's tau" = AC1_tau_values) 
                   
colnames(data) <- c("Salinity","Slope", "Mean NDVI in 1998","Greening trend","Elevation","Precip. variance trend",
                    "Variance Kendall's tau", "AC1 Kendall's tau")

M <- cor(data, use= "complete.obs", method='pearson')
p.mat <- cor_pmat(data, alternative = 'two.sided')
ggcorrplot(M, show.legend = FALSE, show.diag = FALSE, type= "upper", sig.level = 0.05,
           title = "Correlation Matrix Kendall's tau AC1 and variance", p.mat=p.mat, insig= "blank", lab = TRUE)

#Mapping some of the supplementary figures 

plot(slope, col= mapcolours, main = "Slope (in degrees)", legend = TRUE)
lines(catchment_outline)

plot(DEM, col= mapcolours, main = "Elevation (in meters)", legend = TRUE)
lines(catchment_outline)

plot(salinity, col= (hcl.colors(100, palette="YLOrbr")) , main = "Soil salinity (NSI) ", legend = TRUE)
lines(catchment_outline)

plot(greening, col= (hcl.colors(100, palette="BrBG")) , main = "Greening trend (Kendall's tau)", legend = TRUE, range = c(-1,1))
lines(catchment_outline)

plot(NDVI_1998, col= (hcl.colors(100, palette="Greens", rev=TRUE)) , main = "Mean NDVI in 1998", legend = TRUE)
lines(catchment_outline)

plot(precip, col= mapcolours, main = "Precipitation variance trend (Kendall's tau)", legend = TRUE, range = c(-1,1))
lines(catchment_outline)



