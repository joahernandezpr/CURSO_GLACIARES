#### Positive Degree-Day Model
#### 02.09.2021
#### last update 01.03.2023
#### Author: Franziska Temme
########################################################################################################################

#### import required packages
import numpy as np
import xarray as xr

#### set file directories
infile = '../SES7/data/ERA5_input_bell.nc'   # input path
outfile = './output/Bell_SMB_out.nc'               # output path


#### specify parameters
dt = 6                   # time step in hours

M_thresh =  0.0           # Temperature threshold to distinguish between melt and no melt
S_thresh =  1.0           # Temperature threshold to distinguish between rain and snowfall
DDFice = 7.0             # melt factor for ice in mm/d/°C
DDFsnow = 5.0            # melt factor for snow in mm/d/°C

Temp_in_K = False        # Is the input temperature given in Kelvin or degrees Celsius?

#### read in the input data
ds = xr.open_dataset(infile)

if Temp_in_K == False:
    Temp = np.array(ds.T2 - 273.15)
else:
    Temp = np.array(ds.T2)
Prec = np.array(ds.RRR)
MASK = np.array(ds.MASK)
HGT = ds.HGT
lats = ds.lat
lons = ds.lon
time = ds.time

Temp[:,MASK == 0.0] = np.nan   # cut temperature input to glacier


#### create arrays that will be filled in the simulations
SNOWD = (Temp.copy() * 0.0)**2                 # snow depth
ACC = np.zeros(np.shape(Temp))+np.nan          # accumulation
MELT = np.zeros(np.shape(Temp))+np.nan         # melt

#### calculate accumulation from precipitation
ACC[Temp <= S_thresh] = Prec[Temp <= S_thresh] # precipitation below a threshold temperature equals snowfall
ACC[Temp > S_thresh] = 0.0                     # precipitation above a threshold temperature equals zero snowfall

#### prepare melt
MELT[Temp <= M_thresh] = 0.0                   # no melt is possible where we fall below the temperature threshold


########################################################################################################################
#### START SIMULATION
########################################################################################################################
print('STARTING SIMULATION')

#### time loop start
for tt in np.arange(0,len(Temp)-1):

    ## create empty arrays
    MELTSNOWpot = np.zeros(np.shape(MASK)) + np.nan          # potential snow melt
    MELTSNOW = np.zeros(np.shape(MASK)) + np.nan             # snow melt
    MELTICE = np.zeros(np.shape(MASK)) + np.nan              # ice melt
    Tas = np.zeros(np.shape(MASK)) + np.nan                  # temperature that has been used to melt the snow

    MELTSNOW[Temp[tt,:,:] <= M_thresh] = 0.0                 # no snow melt is possible where we fall below the temperature threshold
    MELTSNOWpot[Temp[tt, :, :] <= M_thresh] = 0.0            # no snow melt is possible where we fall below the temperature threshold
    MELTICE[Temp[tt, :, :] <= M_thresh] = 0.0                # no ice melt is possible where we fall below the temperature threshold

    ## calculate potential snow melt (independent of if there is enough snow to be melted)
    MELTSNOWpot[Temp[tt,:,:] > M_thresh] = (dt/24) * DDFsnow * Temp[tt,Temp[tt,:,:] > M_thresh]

    ## If we have more snow than potential snow melt: MELTSNOW = MELTSNOWpot, MELTICE = 0
    MELTSNOW[SNOWD[tt,:,:] >= MELTSNOWpot] = MELTSNOWpot[SNOWD[tt,:,:] >= MELTSNOWpot]
    MELTICE[SNOWD[tt,:,:] >= MELTSNOWpot] = 0.0

    ## If we have less snow than potential snow melt:
    # 1) Melt all the snow: MELTSNOW = SNOWD
    MELTSNOW[SNOWD[tt,:,:] < MELTSNOWpot] = SNOWD[tt,SNOWD[tt,:,:] < MELTSNOWpot]
    # 2) Melt the underlying ice with the "left-over" Temperature (T - Tas)
    Tas[SNOWD[tt,:,:] < MELTSNOWpot] = (24/dt) * SNOWD[tt,SNOWD[tt,:,:] < MELTSNOWpot]/DDFsnow  # get temperature that has been used to melt the snow: Tas
    MELTICE[SNOWD[tt,:,:] < MELTSNOWpot] = (dt/24) * DDFice * (Temp[tt,SNOWD[tt,:,:] < MELTSNOWpot] - Tas[SNOWD[tt,:,:] < MELTSNOWpot])  # melt the ice with the leftover Temp

    ## total melt
    MELT[tt,:,:] = MELTSNOW + MELTICE

    ## update snow depth
    SNOWD[tt+1,:,:] = SNOWD[tt,:,:] + (ACC[tt,:,:] - MELT[tt,:,:])
    SNOWD[tt+1, SNOWD[tt+1,:,:] < 0.0] = 0


#### calculate SMB for whole period
SMB = ACC - MELT

#### set non-glacierized points to 0
SMB[:,MASK == 0] = np.nan
ACC[:,MASK == 0] = np.nan
MELT[:,MASK == 0] = np.nan


########################################################################################################################
#### WRITE OUTPUT
########################################################################################################################
print('WRITING TO FILE')

dsout = xr.Dataset(
    data_vars=dict(
        HGT=(['lat', 'lon'], np.array(HGT), {'units': 'm'}),
        MASK=(['lat', 'lon'], np.array(MASK)),
        ACC=(["time", "lat", "lon"], np.array(ACC) / 1000, {'units': 'm w.e.'}),
        MELT=(["time", "lat", "lon"], np.array(MELT) / 1000, {'units': 'm w.e.'}),
        SMB=(["time", "lat", "lon"], np.array(SMB) / 1000, {'units': 'm w.e.'}),
    ),
    coords=dict(
        time=(["time"], np.array(time)),
        lat=(["lat"], np.array(lats)),
        lon=(["lon"], np.array(lons)),
    ),
)

dsout.to_netcdf(outfile)




