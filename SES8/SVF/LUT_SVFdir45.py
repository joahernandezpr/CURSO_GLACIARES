#### Calculating directed Sky View Factor
#### 22.03.22, last update 02.03.2023
#### Author: Franziska Temme
########################################################################################################################

#### import
import numpy as np
import pandas as pd
import xarray as xr
import math

# Import user defined routines and packages
# loading input data
from relshadFunct_SD import *

#### read in necessary input
ds = xr.open_dataset('/home/christian/curso-glacio-apuandes/SES6/dom/Bell_dom.nc')
DEM = np.array(ds.HGT)
MASK = np.array(ds.MASK)
ASP = np.array(ds.ASPECT)+180
SLO = np.array(ds.SLOPE)
lats = np.array(ds.lat)
lons = np.array(ds.lon)


########################################################################################################################

slo = np.radians(SLO)
asp = np.radians(ASP)
directions = np.arange(0,360,5)

#### empty array to fill
LUT_SVFdir = np.zeros((len(directions),len(lats),len(lons)))
i = 0


############ START

for DIR in directions:    # go trough all wind directions in 5 deg steps
    print(DIR)
    res = DEM*0
    count = 0

    # azimuth with winddireciton +/- 45 °
    AZI =  np.array([DIR-20, DIR-15, DIR-10, DIR-5, DIR, DIR+5, DIR+10, DIR+15, DIR+20])

    # check that AZI is between 0 and 360
    AZI[AZI < 0] = AZI[AZI < 0] + 360
    AZI[AZI >= 360] = AZI[AZI >= 360] - 360

    EL = np.arange(2,90,10)

    for azi in AZI:
        for el in EL:
            illu, Hmax = relshad(DEM,MASK,lats,lons,el,azi)
            a = ((math.cos(np.radians(el)) * np.sin(slo) * np.cos(asp - np.radians(azi))) + (np.sin(np.radians(el)) * np.cos(slo)))
            a[a < 0] = 0
            a[a > 0] = 1
            a[illu == 0] = 0
            res = res + a
            count = count + 1

    vsky = DEM*0
    vsky[:,:] = np.nan
    vsky[MASK == 1] = res[MASK == 1]/(len(AZI)*len(EL))

    LUT_SVFdir[int(DIR/5),:,:] = vsky
    #LUT_SVFdir[i,:,:] = vsky
    #i = i+1



LUT_SVFdir = xr.DataArray(LUT_SVFdir, dims = ['count','lat', 'lon'])



### writing out
ds_SVF = xr.Dataset()
ds_SVF.coords['lon'] = ds.lon.values
ds_SVF.lon.attrs['standard_name'] = 'lon'
ds_SVF.lon.attrs['long_name'] = 'longitude'
ds_SVF.lon.attrs['units'] = 'degrees_east'

ds_SVF.coords['lat'] = ds.lat.values
ds_SVF.lat.attrs['standard_name'] = 'lat'
ds_SVF.lat.attrs['long_name'] = 'latitude'
ds_SVF.lat.attrs['units'] = 'degrees_north'

ds_SVF.coords['count'] = directions
ds_SVF.lat.attrs['standard_name'] = 'count'

ds_SVF['HGT'] = ds.HGT
ds_SVF['MASK'] = ds.MASK
ds_SVF['SLOPE'] = ds.SLOPE
ds_SVF['ASPECT'] = ds.ASPECT+180
ds_SVF['SVFdir'] = LUT_SVFdir


ds_SVF.to_netcdf('./LUT_SVFdir45.nc')
