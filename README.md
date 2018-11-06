# HALO_lidar_toolbox
HALO Streamline Photonics lidar data postprocessing toolbox

----
### Recommended, and partly compulsory sequence for processing HALO Doppler lidar data by using this toolbox:
1) halo_config.txt
2) calibrateHALO
3) calculateHALOwindvadProduct **AND/OR** calculateHALOwinddbsProduct
4) calculateHALOwStatsProduct
5) calculateHALOverticalTKEproduct
6) calculateHALOwindShearProduct
7) calculateHALOcloudProduct
8) calculateHALOatmBLclassificationProduct
9) calculateHALObetaVeloCovarianceProduct

----
### DETAILED INSTRUCTIONS

More information of the matlab functions can be found by typing >>>  help nameofthefunction  <<<


### 1) Fill in the halo_config.txt file for your site.

Check the paths!

If instrument parameters are changed during deployment, first specify the date from when
the changes are to be taken into account. Add at the very end of our site's parameter list
a new field specifying that date:

parameters_valid_from_including = YYYYMMDD

Then, add the all of the parameters which have changed right below the above mentioned field.
For example, if HALO unit was replaced and it has a different focus range:

halo_unit_ID = XX
focus_stare_co = 1500


### 2) Calibrate HALO data by using the calibrateHALO.m function:

calibrateHALO('site',[YYYMMDD YYYYMMDD])

See help calibrateHALO. The functions reads all of the data that is specified in the halo_config.txt
file, corrects the background (with ripple correction if *background*.txt files are available),
correct focus (currently at specified site only) and writes a netcdf file per day per measurement
mode into their respective specified paths.


### 3) Calculate winds with uncertainties

  #### 3.1) If VAD/PPI wind scans are available, calculate them by using calculateHALOwindvadProduct.m function:

  calculateHALOwindvadProduct('site',[YYYMMDD YYYYMMDD],'NN')

  where 'NN' is the elevation angle of the scan; if 75 degrees then type '75' (string input).

  #### 3.2) If DBS winds are available, calculate them by using calculateHALOwinddbsProduct.m function:

  calculateHALOwinddbsProduct('site',[YYYMMDD YYYYMMDD],'noofbeams')

  where 'noofbeams' is the number of beams in the DBS scans; if 3 beams then type '3beams' (string input).


### 4) Calculate vertical velocity statistics

calculateHALOwStatsProduct('site',[YYYMMDD YYYYMMDD])

By default, the function calculates the following quantitites from vertically pointing measurements at 3 and 60 min resolutions:
- radial velocity mean, std. deviation, variance, skewness and kurtosis with respective standard errors
- attenuated backscatter coefficient mean and variance with respective standard errors
- signal (SNR+1) mean and variance with respective standard errors
- instrumental precision mean and variance

### 5) calculateHALOverticalTKEproduct

calculateHALOverticalTKEproduct('site',[YYYMMDD YYYYMMDD],'windproduct','typeof')

The dissipation rate of turbulent kinetic energy (TKE) can be calculated with using either 1) vad or 2) dbs:
1) 'windproduct' and 'typeof' are 'windvad' and 'eleNN', respectively with NN being the elevation angle in degrees (e.g. 'ele75' or 'ele09')
2) 'windproduct' and 'typeof' are 'windbs' and 'Nbeams', respectively with N specifying the number of dbs beams (e.g. '3beams')

### 6) calculateHALOwindShearProduct


28 October 2018,
Antti J Manninen
