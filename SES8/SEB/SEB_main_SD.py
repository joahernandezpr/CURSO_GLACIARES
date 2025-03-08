#!/usr/bin/python
# -*- coding: utf-8 -*-

## Simplified energy balance model
## Author:      Franziska Temme, Johannes Fürst
## Last update: 05.04.2023
################################################################################################


########################################################################################################################
# IMPORT EXTERNAL PACKAGES
########################################################################################################################

# Import general Python packages

import os
import numpy as np
import importlib
from datetime import datetime

cwd = os.getcwd()
print('Your current working directory is :')
print('   ', cwd)


# Import user defined routines and packages
# loading input data
from SEB_IO_SD import *

# model parameter setttings
from SEB_param_SD import *


########################################################################################################################
#  PREPROCESSING
########################################################################################################################


##### SOLID PRECIPITATION

snow = 1.0*Prec
snow[T2>temp_thresh] = 0.0         # where temperature is below threshold, preciptation is solid

snow[:,MASK == 0] = np.nan         # only interested over the glacier


#### SNOWDRIFT
Snowcorr = np.zeros(np.shape(snow))+np.nan   # empty array

Emin = np.min(HGT)                           # minimum altitude
Emax = np.max(HGT)                           # maximum altitude
E = ((HGT - Emin) / (Emax - Emin))           # altitude scaling
E2 = np.array(E.copy())
E2[E2<0] = 0.0

print('START SNOWDRIFT')
for i in np.arange(0,len(snow)):             # loop through all time steps

    dir = 5 * int(DIR[i]/5)                  # round directions to 5 deg sectors
    SVF = np.array(LUT_SVF.sel(count = dir)) # directed sky-view factor for the current wind direction

    Cwind = (WS[i] / 4.56) * E2 * (Dmax * (1 - SVF) - 1) + 0.0  # correction field

    Snowcorr[i,:,:] = (snow[i,:,:] + Cwind * snow[i,:,:])       # snowfall amount is corrected accordingly

Snowcorr[Snowcorr < 0.0] = 0.0                # avoid too much  snow been blown away (more than fallen)


########################################################################################################################
#  SMB MODEL
########################################################################################################################

tmax = len(time)

## create arrays
albedo     = np.zeros(np.shape(T2))+albedo_ice  # ice albedo
Qm         = np.zeros(np.shape(T2))             # surface energy balance (W m**-2)
snow_depth = np.zeros(np.shape(T2))             # snow depth (m w.e.)
smb_cum    = np.zeros(np.shape(T2))             # cumulative surface mass balance (m w.e.)
smb        = np.zeros(np.shape(T2))             # surface mass balance (m w.e.)
abl        = np.zeros(np.shape(T2))             # surface abaltion (m w.e.)
acc        = np.zeros(np.shape(T2))             # surface accumulation (m w.e.)

## only needed over the glacier mask
albedo[:,MASK == 0] = np.nan
Qm[:,MASK == 0] = np.nan
snow_depth[:,MASK == 0] = np.nan
smb[:,MASK == 0] = np.nan
smb_cum[:,MASK == 0] = np.nan
abl[:,MASK == 0] = np.nan
acc[:,MASK == 0] = np.nan

## create arrays
a = np.arange(0,st_p_day)
Tmax_day = np.zeros((int(len(T2)/st_p_day),len(lats),len(lons)))   # maximum temperature
snow_day = np.zeros(np.shape(Tmax_day)) + np.nan                   # daily accumulated snow
tacc_day = np.zeros(np.shape(Tmax_day))                            # daily accumulated temperature

print('START ALBEDO PARAMETRIZATION')

########################################################################################################################
#  ALBEDO PARAMETRIZATION
########################################################################################################################

## daily snow sum and T max
for dd in np.arange(0,len(Tmax_day)):
    Tmax_day[dd,:,:] = T2[a+(dd*st_p_day),:,:].max(axis = 0)
    snow_day[dd,:,:] = Snowcorr[a+(dd*st_p_day),:,:].sum(axis = 0)

tacc_day[np.isnan(snow_day)] = np.nan

## accumulated Tmax since last snowfall
for ii in np.arange(0,len(Tmax_day)):
    tacc_day[ii,snow_day[ii,:,:] > 0.0] = 0.0
    tacc_day[ii,(snow_day[ii,:,:]==0) & (Tmax_day[ii,:,:] > 0.0)] = tacc_day[ii-1,(snow_day[ii,:,:]==0) & (Tmax_day[ii,:,:] > 0.0)] + Tmax_day[ii,(snow_day[ii,:,:]==0) & (Tmax_day[ii,:,:] > 0.0)]
    tacc_day[ii,(snow_day[ii] == 0) & (Tmax_day[ii] <= 0.0)] = tacc_day[ii-1,(snow_day[ii] == 0) & (Tmax_day[ii] <= 0.0)]

tacc_day = tacc_day + 1.0

## compute snow albedo
# daily
alb_s = 0.9 - p2 * np.log10(tacc_day)

# snow albedo to model time step
alb_s_3h = np.repeat(alb_s,st_p_day, axis=0)

#############
# TIME LOOP START
# of the mass balance model
# time index : tt
# time step  : dt
# end time   : tmax

print('START SEB SIMULATION')

for tt in np.arange(0,tmax-1):

    # get solid accumulation / snowfall
    acc[tt,:,:] = Snowcorr[tt,:,:]

    # adjust albedo: if snowfree --> albedo ice, if snow --> albedo snow
    albedo[tt, snow_depth[tt,:,:] == 0] = albedo_ice
    albedo[tt, snow_depth[tt,:,:] > 0.0] = alb_s_3h[tt, snow_depth[tt,:,:] > 0.0]

    # calculate surface energy balance
    Qm[tt,:,:] = (1 - albedo[tt,:,:]) * tau * G[tt,:,:] + c1 * T2[tt,:,:] + c0

    abl[tt,Qm[tt,:,:] > 0.0] = -Qm[tt,Qm[tt,:,:] > 0.0] * sec_per_hr * dt / lm / rho_water

    abl[tt,Qm[tt,:,:] < 0.0] = 0.0

    smb[tt,:,:] = acc[tt,:,:] + abl[tt,:,:]


    #### TIME STEP UPDATE : tt + 1

    # cumulative mass balance
    smb_cum[tt + 1,:,:] = smb_cum[tt,:,:] + smb[tt,:,:]

    # calculate snowdepth
    snow_depth[tt+1,:,:] = snow_depth[tt,:,:] + (acc[tt,:,:] + abl[tt,:,:])

    # the snowdepth cannot be smaller than zero
    snow_depth[tt+1,snow_depth[tt+1,:,:] < 0.0] = 0


########################################################################################################################
#  WRITE OUTPUT            #
########################################################################################################################
print('WRITING OUTPUT TO FILE')

ALB = xr.DataArray(albedo, dims=['time','lat','lon'], attrs=dict(long_name='Albedo'))
ds['alpha'] = ALB

QM = xr.DataArray(Qm, dims=['time','lat','lon'], attrs=dict(units='W/m2', long_name='Melt energy'))
ds['Qm'] = QM

ACC = xr.DataArray(acc, dims=['time','lat','lon'], attrs=dict(units='m w.e.', long_name='Accumulation'))
ds['acc'] = ACC

ABL = xr.DataArray(abl, dims=['time','lat','lon'], attrs=dict(units='m w.e.', long_name='Ablation'))
ds['abl'] = ABL

SMB = xr.DataArray(smb, dims=['time','lat','lon'], attrs=dict(units='m w.e.', long_name='Surface mass balance'))
ds['smb'] = SMB

ds.to_netcdf(OUTPUT_NAME)

