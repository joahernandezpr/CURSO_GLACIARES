#!/usr/bin/python
# -*- coding: utf-8 -*-

## Simplified energy balance model
## Author:      Franziska Temme, Johannes Fürst
## Last update: 05.04.2023
########################################################################################################################

## Parameter defintions

## Accessed by:
##    SEB_main_3D.py


## constants
sec_per_hr = 3600.0        # seconds per hour    
lm         = 334000.       # latent heat of melting (J/kg)
rho_water  = 1000.         # density of fresh water (kg/m3)
rho_ice    = 917.          # density of ice (kg/m3)


## general model settings
dt        = 6              # model time step (hours)
st_p_day = int(24/dt)          # how many steps per day do we have?

## input data parameters
temp_thresh =  1.8         # temperature threshold for solid prec.  (˚C )


## SNOWDRIFT model parameters
Dmax = 8                   # maximum deposition


## SEB model parameters

potential_radiation  = False   # for using the model with potential radiation set to True, with global radiation set to False
transmissivity       = 0.38    # atmospheric transmissivity

if potential_radiation == False:
    tau = 1.0
else:
    tau = transmissivity

c1                   =  30     # temperature melt factor (W m**-2 K**-1) 
c0                   = -40     # empirical factor (W m**-2)

albedo_snow          =  0.9    # albedo of snow
albedo_ice           =  0.4    # albedo of ice
p2                   = 0.155   # empirical parameter (Gabbi et al. 2014, Journal of Glaciology)
