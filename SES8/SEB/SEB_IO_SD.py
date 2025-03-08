#!/usr/bin/python
# -*- coding: utf-8 -*-

## Simplified energy balance model
## Author:      Franziska Temme, Johannes Fürst
## Last update: 05.04.2023
########################################################################################################################


########################################################################################################################
####  IMPORT EXTERNAL PACKAGES 

## Import general Python packages
import numpy as np
import pandas as pd
import xarray as xr

#### READ IN DATA   
########################################################################################################################

## Input - Output path
INPUT_NAME = '../data/ERA5_G_input_bell.nc'
OUTPUT_NAME = '../output/Bell_SMB_out_SEB.nc'

## read data
ds = xr.open_dataset(INPUT_NAME)      # input netcdf file
T2 = np.array(ds.T2)-273.15           # air temperature (converted from °C to K)
G = np.array(ds.G)                    # radiation (potential or global, adjust in SEB_param)
Prec = np.array(ds.RRR)/1000          # precipitation (converted from mm to m)
MASK = np.array(ds.MASK)              # glacier mask
lats = np.array(ds.lat)        
lons = np.array(ds.lon)
time = np.array(ds.time)
HGT = ds.HGT                          # DEM

#### read in SVF
dsSVF = xr.open_dataset('../SVF/LUT_SVFdir45.nc')    # path to the look-up table of directed sky-view factors
LUT_SVF = dsSVF.SVFdir                # look-up table of directed sky-view factors (can be created with the script LUT_SVFdir45.py)

#### read in Wind data
Wind = pd.read_csv('../data/wind_data.csv', sep='\t')

DIR = np.array(Wind.DIR)              # timeseries of wind direction (°)
WS = np.array(Wind.WS)                # timeseries of wind speed (m/s)

