# MSc-Thesis-
Scripts for MSc thesis "Relations between temporal resilience indicators and trend breakpoints in a dryland high-mountain catchment" by Willem Grootoonk

Utrecht University (01/07/2024)

##DATA SOURCES###

All of the used data is provided here: https://zenodo.org/records/12759684?token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImIxOTVhNDZmLTM5ZTQtNDZiMi04OTYzLTBkYTAwNDFmYWE1OCIsImRhdGEiOnt9LCJyYW5kb20iOiIwZDJiNzc3ZTk5ZjlhMWRmMjU5ZGQ5MTk4OWQxNzE2ZiJ9.Dr60es82tjUCwx2SKwZ7dYLKRgk_6Q_ayePh9kEdl-V-Kbc1EySzP_21lNAfMkZLVQHvoEVVXJH50y1YsZrcfQ 
and can be used with proper referencing 

##HOW TO USE THE SCRIPTS###

Three main R scripts were used for this thesis, and one script for Google Earth Engine.
The operation of these scripts should be clear from comments in the scripts themselves. 
Feel free to use these scripts for further research (with proper referencing) 

## thesis_script_functions: Main script that includes the functions to compute resillience indicators at pixel level.
#Input: 
- A clean and prepared temporal vegetation raster brick (i.e. a stack of rasters forming a timeseries for each pixel). I use one produduced by Vermeer, A. L. (2021). Ecological stability in the face of climatic disturbances: a case study of a dryland ecosystem in the Moroccan High Atlas Mountains (Master's thesis). Any environmental dataset of sufficient length will work however.
- A catchment (study area) outline , channels mask, and plantation mask produced by Vermeer (2021). Optional, will work without.
- A CHIRP dataset of precipitation for the area created using thesis_script_earthengine. Optional, will work without.

#Output:
- Resillience indicator rasters with per-pixel results for Kendall's tau corresponding with the chosen indicator and parameters. By defaults, these are standard deviation (variance) and lag-one autocorrelation.
- A raster with Kendall's tau for vegetation trend corresponding with vegetation trends in the same time period
- A raster with Kendall's tau for precipitation variability trend in the same time period

## thesis_script_dataquality: Script that computes figures to analyize data quality and that contains a null model to produce correction rasters to account for missing values 
#Input: 
- A clean and prepared temporal vegetation raster brick (i.e. a stack of rasters forming a timeseries for each pixel).
- A catchment outline, channels mask, and plantation mask produced by Vermeer (2021). Optional, will work without.

#Output:
- Various figures related to data quality
- Correction rasters for chosen indicators and parameters that can be used to correction for the effect of missing values on the output of the main function.
- Various figures describing the correction rasters

## thesis_script_figures: Script that produces figures, maps and does some comparative analysis. 
#Input: 
- All outputs of the above scripts
- A Breakpoint typology raster created by Vermeer (2021). This is specific to this research, can be skipped.
- A NIR-SWIR salinity index raster for the study area produced using thesis_script_earthengine This is specific to this research, can be skipped.
- Data for slope, elevation, and the cathcment outline from Vermeer (2021). This is specific to this research, can be skipped.

#Output:
- Various figures
  
## thesis_script_earthengine: Script that is used to download CHRIPS data and to produce a salinity map 
#Input: 
- Shapefile for the study area
- Landsat 5 data from EE
- CHIRPS data from EE

#Output:
- A NIR-SWIR salinity index raster for the study area
- A series of CHIRPS images for the study area. Honestly, this is not a good way to import these, would be way better to use some dedicated export function.

