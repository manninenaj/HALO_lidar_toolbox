# HALO_lidar_toolbox
HALO Streamline Photonics Doppler lidar data postprocessing toolbox

This toolbox has been succesfully tested at several sites: Finnish Dopper lidar network site, ARM Dopper lidar network 
and campaign site, Jülich (Germany), Granada (Spain). However, please note that the toolbox is still in its 
development phase and errors mights occur. Let me know if that happens. You are also more than welcome to suggest 
changes or additions in to the toolbox! 

Antti Manninen
University of Helsinki, Finland


8 November 2018,
Antti J Manninen
University of Helsinki, Finland
antti.j.manninen(at)helsinki.fi

----
### Sequence for processing HALO Doppler lidar data by using this toolbox:
1) halo_config.txt
2) calibrateHALO
3) calculateHALOwindvadProduct **AND/OR** calculateHALOwinddbsProduct
4) calculateHALOwindShearProduct
5) calculateHALOwStatsProduct
6) calculateHALOverticalTKEproduct
7) calculateHALOcloudProduct
8) calculateHALOatmBLclassificationProduct
9) calculateHALObetaVeloCovarianceProduct

----
### DETAILED INSTRUCTIONS

More information of the matlab functions can be found by typing >>>  help nameofthefunction  <<<

### 1) Fill in the halo_config.txt file for your site.

Check the paths!

If instrument parameters are changed during deployment, first specify the date from when
the changes are to be taken into account. Add at the very end of your site's parameter list
a new field specifying that date:

parameters_valid_from_including = YYYYMMDD

Then, add the all of the parameters which have changed right below the above mentioned field.
For example, if HALO unit was replaced and it has a different focus range:

halo_unit_ID = XX

focus_stare_co = 1500

### 2) Calibrate HALO data by using the calibrateHALO.m function:

calibrateHALO('site',[YYYMMDD YYYYMMDD])

The function reads all of the data, which is available for the site and date as specified in the halo_config.txt
file. It then corrects the background (with ripple correction if *background*.txt files are available),
corrects focus (currently at specified sites only), and writes a netcdf file per day and per measurement
mode into their respective specified paths.

### 3) Calculate winds with uncertainties

  #### 3.1) If VAD/PPI wind scans are available, calculate them by using calculateHALOwindvadProduct.m function:

  calculateHALOwindvadProduct('site',[YYYMMDD YYYYMMDD],'NN')

  where 'NN' is the elevation angle of the scan; if 75 degrees then type '75' (string input).

  The function reads ppi files and calculates (u,v,w) wind components, wind speed, wind direction, respective errors 
  due to random instumental noise, and overall errors using VAD tehcnique (Päschke et al. 2015; Newsom et al. 2017), 
  and writes the retrieved winds into a netcdf file.
   
  #### 3.2) If DBS winds are available, calculate them by using calculateHALOwinddbsProduct.m function:

  calculateHALOwinddbsProduct('site',[YYYMMDD YYYYMMDD],'noofbeams')

  where 'noofbeams' is the number of beams in the DBS scans; if 3 beams then type '3beams' (string input).

  The function reads dbs files and calculates (u,v,w) wind components, wind speed, wind direction and writes the 
  retrieved winds into a netcdf file.

### 4) calculateHALOwindShearProduct

calculateHALOwindShearProduct('site',[YYYMMDD YYYYMMDD],'windproduct','typeof')

The vector wind shear can be calculated with using either 1) vad or 2) dbs (depending on what is available):
1) 'windproduct' and 'typeof' are 'windvad' and 'eleNN', respectively with NN being the elevation angle in degrees (e.g. 'ele75' or 'ele09')
2) 'windproduct' and 'typeof' are 'windbs' and 'Nbeams', respectively with N specifying the number of dbs beams (e.g. '3beams')

The function calculates vector wind shear ()in temporal resolution(s) based on what is/are available in the vertical 
velocity statistics product, and writes the results into a netcdf file.

### 5) Calculate vertical velocity statistics

calculateHALOwStatsProduct('site',[YYYMMDD YYYYMMDD])

By default, the function calculates the following quantitites from vertically pointing measurements at 3 and 60 min resolutions:
- radial velocity mean, std. deviation, variance, skewness and kurtosis with respective standard errors
- attenuated backscatter coefficient mean and variance with respective standard errors
- signal (SNR+1) mean and variance with respective standard errors
- instrumental precision mean and variance

### 6) calculateHALOverticalTKEproduct

calculateHALOverticalTKEproduct('site',[YYYMMDD YYYYMMDD],'windproduct','typeof')

The dissipation rate of turbulent kinetic energy (TKE) can be calculated with using either 1) vad or 2) dbs wind product (depending on what is available):
1) 'windproduct' and 'typeof' are 'windvad' and 'eleNN', respectively with NN being the elevation angle in degrees (e.g. 'ele75' or 'ele09')
2) 'windproduct' and 'typeof' are 'windbs' and 'Nbeams', respectively with N specifying the number of dbs beams (e.g. '3beams')

The function calculates the dissipation rate of turbulent kinetic energy directly from vertical velocity variance 
(O'Connor et al 2010) with temporal resolution the vertical velocity statistics product is provided, and writes the 
results into a netcdf file.

### 7) calculateHALOcloudProduct

calculateHALOcloudProduct('site',[YYYMMDD YYYYMMDD])

The function calculates cloud base height, cloud base velocity, and provides cloud-precipitation mask in temporal
resolution(s) based on what is/are available in the vertical velocity statistics product.

### 8) calculateHALOatmBLclassificationProduct

calculateHALOatmBLclassificationProduct('site',[YYYMMDD YYYYMMDD])

The function generates boundary layer classification product from calculated Doppler lidar quantities, and writes the 
product into a netcdf file. See detailed description of the method in Manninen et al. (2018), doi: 
10.1029/2017JD028169.

### 9) calculateHALObetaVeloCovarianceProduct

calculateHALObetaVeloCovarianceProduct('site',[YYYMMDD YYYYMMDD])

The function calculates covariance between the attenuated backscatter coefficient and vertical velocity, which are 
read from the vertical velocity statistics product by using a default window size of 90 min and 6 range bins, and outputs the covariance with standard errors and confidence interval into a netcdf file.
