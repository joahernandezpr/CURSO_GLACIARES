#### Script with function for relief shading
#### including a function to calculate distances between 2 lat/lon points
#### output: grid ("illu") with 0=shaded, 1=sun)
#### 26.03.21 Franziska Temme
########################################################################################################################
#### INPUT: dem: elevation model (m) (row1 = north boundary); mask: mask with glacier points = 1
#### lats: array with latitudes (decreasing); lons: array with longitudes (increasing)
#### solh: solar elevation (deg); sdirfn: illumination direction from north (deg)
########################################################################################################################
##### START FUNCTION #####
########################################################################################################################
#### import
import numpy as np
import pandas as pd
import xarray as xr
import math

#### HAVERSINE FUNCTION for calculation of distance between two points
def haversine(lat1,lon1,lat2,lon2):
    lat1_rad = math.radians(lat1)
    lat2_rad = math.radians(lat2)
    lon1_rad = math.radians(lon1)
    lon2_rad = math.radians(lon2)
    delta_lat = lat2_rad-lat1_rad
    delta_lon = lon2_rad-lon1_rad
    a = ((math.sin(delta_lat/2))**2 + math.cos(lat1_rad)*math.cos(lat2_rad)*(math.sin(delta_lon/2))**2)**0.5
    d = 2*6371000*math.asin(a)
    return d

#### RELIEF SHADING FUNCTION
def relshad(dem,mask,lats,lons,solh,sdirfn):

    #### map features
    z = dem
    Hmax = dem*0.0   # additionally write out table with Hmax (not necessary - remove when everything works fine)
    illu = dem*0.0   # create grid that will be filled
    illu[:,:] = np.nan   # set to nans

    rmax = ((np.linalg.norm(np.max(lats)-np.min(lats)))**2 + (np.linalg.norm(np.max(lons)-np.min(lons)))**2)**0.5  # define max. radius (that covers DEM area) in degrees lat/lon
    nums = int(rmax * len(lats) / (lats[-1] - lats[0]))  # number of points for similar resolution as input (e.g. 200 m)

    ##### calculate direction to sun
    beta = math.radians(90-sdirfn)  # beta = 90-sdirfn
    dy = math.sin(beta)*rmax    # walk into sun direction (y) as far as rmax
    dx = math.cos(beta)*rmax    # walk into sun direction (x) as far as rmax


    ##### EXTRACT PROFILE TO SUN FROM EACH GRID POINT
    # but only for glacier grid cells
    for ilat in np.arange(1,len(lats)-1,1):
    #for ilat in np.arange(6,8,1):
        #print(ilat)
        for ilon in np.arange(1,len(lons)-1,1):
            if mask[ilat, ilon] == 1:
                #print(ilon)
                start = (lats[ilat], lons[ilon])
                targ = (start[0] + dy, start[1] + dx)  # find target position

                # nums = int(rmax * len(lats)/(lats[0]-lats[-1]))     # number of points for similar resolution as input (e.g. 200 m)

                ## make points along profile (lat/lon)
                lat_list = np.linspace(start[0], targ[0], nums)  # equally spread points along profile
                lon_list = np.linspace(start[1], targ[1], nums)  # equally spread points along profile

                ## don't walk outside DEM boundaries
                lat_list_short = lat_list[(lat_list < max(lats)) & (lat_list > min(lats))]
                lon_list_short = lon_list[(lon_list < max(lons)) & (lon_list > min(lons))]

                ## cut to same extent
                if (len(lat_list_short) > len(lon_list_short)):
                    lat_list_short = lat_list_short[0:len(lon_list_short)]
                if (len(lon_list_short) > len(lat_list_short)):
                    lon_list_short = lon_list_short[0:len(lat_list_short)]

                ## find indices (instead of lat/lon) at closets gridpoint
                idy = (ilat, (np.abs(lats - lat_list_short[-1])).argmin())
                idx = (ilon, (np.abs(lons - lon_list_short[-1])).argmin())

                ## make points along profile (indices)
                y_list = np.round(np.linspace(idy[0], idy[1], len(lat_list_short)))
                x_list = np.round(np.linspace(idx[0], idx[1], len(lon_list_short)))

                ########################################################################################################################
                #### Calculate ALTITUDE along profile
                zi = z[y_list.astype(np.int), x_list.astype(np.int)]

                ##### Calclulate DISTANCE along profile
                d_list = []
                for j in range(len(lat_list_short)):
                    lat_p = lat_list_short[j]
                    lon_p = lon_list_short[j]
                    dp = haversine(start[0], start[1], lat_p, lon_p)
                    d_list.append(dp)
                distance = np.array(d_list)

                #### get topography angle
                # print(zi, distance)
                Hang = np.degrees(np.arctan((zi[1:len(zi)] - zi[0]) / distance[1:len(distance)]))
                # print(Hang)
                # print(np.max(Hang))
                Hmax[ilat, ilon] = np.max(Hang)  ## remove

                if np.max(Hang) > solh:
                    illu[idy[0], idx[0]] = 0
                else:
                    illu[idy[0], idx[0]] = 1

            #else:
                #print('NO GLACIER CELL')

    return (illu,Hmax)

