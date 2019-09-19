# HALO lidar toolbox
HALO Streamline Photonics Doppler lidar data postprocessing toolbox

This software is released under DRAPL license, for now. See the CRAPL-LICENSE.txt for more details.

This toolbox has been succesfully tested at several sites: FMI Dopper lidar network sites, ARM sites where Doppler 
lidar was deployed, Jülich (Germany), and Granada (Spain). However, note that the toolbox is still in development 
phase and errors might occur. Please let me know if that happens. You are also more than welcome to suggest changes or 
additions in to the toolbox! 

Please note that the toolbox requires that the original hpl files are converted into netcdf format!


Height - range - height_agl - altitude - altitude_instrument notation:
-----------------------------------------------------------------------
                 
         - - - - - - - - - - - - - - - - - - 
          \     |               |   |
      range\    |range          |   |
            \   |               |   |
             \  |     height_agl|   |height
              \ |               |   |  
               ____             |   |
              _|__|_ - - - - - -|- -|- - -
        lidar |    |            |   |   |
           ___|____|___         |   |   |
           | __ __ __ |         |   |   |
           | || || || |         |   |   |
  building | __ __ __ |         |   |   |
           | || || || |         |   |   |altitude_instrument 
           | __ __ __ |         |   |   |
           | || || || |         |   |   |
      _____|__________|_________|_  |   | 
      ground level          |       |   |                           
                    altitude|       |   |  
                            |       |   |  
                           ~~~~~~~~~~~~~~~~~~~ 
                           mean sea level
________________________________ 
Created 2019-11-08, last modified 2019-09-17  
Antti J Manninen  
Finnish Meteorological Institute  
antti.manninen(at)fmi.fi  

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
file. It then corrects the background artefacts (if present) in the signal with a method described by Manninen et al., 
(2016, doi:10.5194/amt-9-817-2016) and if *background*.txt files are available, the functions does a ripple removal as 
described by Vakkari et al. (2018, doi:10.5194/amt-2018-323), and finally corrects focus (currently at specified sites 
only). The instrumental precision of radial velocities are (re-)estimated with the method given by Rye and Hardesty 
(1997, doi:10.1364/AO.36.009425), and Pearson et al. (2009, doi:10.1175/2008JTECHA1128.1), and attenuated backscatter 
coefficients with uncertainties are also re-calculated. Calibrated data are written into a netcdf file per day and per 
measurement mode into their respective specified paths.

### 3) Calculate winds with uncertainties

  #### 3.1) If VAD/PPI wind scans are available, calculate them by using calculateHALOwindvadProduct.m function:

  calculateHALOwindvadProduct('site',[YYYMMDD YYYYMMDD],'NN')

  where 'NN' is the elevation angle of the scan; if 75 degrees then type '75' (string input).

  The function reads ppi files and calculates (u,v,w) wind components, wind speed, wind direction, respective errors 
  due to random instumental noise, and overall errors using VAD technique with formal methods given by Päschke et al. 
  (2015, doi:10.5194/amt-8-2251-2015) and Newsom et al. (2017, doi:10.5194/amt-10-1229-2017). The product is written 
  into a netcdf file.
   
  #### 3.2) If DBS winds are available, calculate them by using calculateHALOwinddbsProduct.m function:

  calculateHALOwinddbsProduct('site',[YYYMMDD YYYYMMDD],'noofbeams')

  where 'noofbeams' is the number of beams in the DBS scans; if 3 beams then type '3beams' (string input).

  The function reads dbs files and calculates (u,v,w) wind components, wind speed, wind direction and writes the 
  retrieved winds into a netcdf file.

### 4) calculateHALOwindShearProduct

calculateHALOwindShearProduct('site',[YYYMMDD YYYYMMDD],'windproduct','typeof')

The vector wind shear can be calculated with using either 1) vad or 2) dbs (depending on what is available):
1) 'windproduct' and 'typeof' are 'windvad' and 'eleNN', respectively with NN being the elevation angle in degrees 
(e.g. 'ele75' or 'ele09')
2) 'windproduct' and 'typeof' are 'windbs' and 'Nbeams', respectively with N specifying the number of dbs beams (e.g. 
'3beams')

The function calculates vector wind shear (ICAO, 2005; https://www.skybrary.aero/bookshelf/books/2194.pdf) in temporal 
resolution(s) based on what is/are available in the vertical velocity statistics product, and writes the results into 
a netcdf file.

### 5) Calculate vertical velocity statistics

calculateHALOwStatsProduct('site',[YYYMMDD YYYYMMDD])

By default, the function calculates the following quantitites from vertically pointing measurements at 3, 30, and 60 
min resolutions. The statistics which are unbiased by random noise and sample size are calculated as given by 
Rimoldini, (2014, doi:10.1016/j.ascom.2014.02.001), and standard errors are estimated with a bootstrap method 
described by Kleiner et al. (2014, doi:10.1111/rssb.12050).
- unbiased mean, std. dev., variance, skewness, and kurtosis of radial velocity with their respective standard errors
- unbiased mean and variance of attenuated backscatter coefficient with their respective standard errors 
- unbiased mean and variance of signal (SNR+1) with their respective standard errors 
- mean and variance of radial velocity and signal instrumental precision

### 6) calculateHALOverticalTKEproduct

calculateHALOverticalTKEproduct('site',[YYYMMDD YYYYMMDD],'windproduct','typeof')

The dissipation rate of turbulent kinetic energy (TKE) can be calculated with using either 1) vad or 2) dbs wind 
product (depending on what is available):
1) 'windproduct' and 'typeof' are 'windvad' and 'eleNN', respectively, with NN being the elevation angle in degrees 
(e.g. 'ele75' or 'ele09')
2) 'windproduct' and 'typeof' are 'windbs' and 'Nbeams', respectively with N specifying the number of dbs beams (e.g. 
'3beams')

The function calculates the dissipation rate of turbulent kinetic energy directly from vertical velocity variance 
(O'Connor et al., 2010, doi:10.1175/2010JTECHA1455.1) with temporal resolution the vertical velocity statistics 
product is provided, and writes the results into a netcdf file.

### 7) calculateHALOcloudProduct

calculateHALOcloudProduct('site',[YYYMMDD YYYYMMDD])

The function calculates cloud base height, cloud base velocity, and provides cloud mask in temporal resolution(s) 
based on what is/are available in the vertical velocity statistics product.

### 8) calculateHALOatmBLclassificationProduct

calculateHALOatmBLclassificationProduct('site',[YYYMMDD YYYYMMDD])

The function generates boundary layer classification product from calculated Doppler lidar quantities, and writes the 
product into a netcdf file. Plese refer to Manninen et al. (2018, doi:10.1029/2017JD028169) for an in-detail 
description of the method and product.

### 9) calculateHALObetaVeloCovarianceProduct

calculateHALObetaVeloCovarianceProduct('site',[YYYMMDD YYYYMMDD])

The function calculates covariance between the attenuated backscatter coefficient and vertical velocity, which are 
read from the vertical velocity statistics product. A window of 90 min and 6 range bins is used to estimate the 
covariance, standard errors, and confidence intervals. The function uses the method given by Engelmann et al., (2008, 
doi:10.1175/2007JTECHA967.1). Standard errors and confidence intervals are estimated with a bootstrap method described 
by Kleiner et al. (2014, doi:10.1111/rssb.12050). Product is written into a netcdf file.
