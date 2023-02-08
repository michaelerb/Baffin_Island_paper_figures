#==============================================================================
# Some functions to support the other scripts.
#    author: Michael Erb
#    date  : 2/7/2023
#==============================================================================

import numpy as np
import math

# This function takes an array and computes the global-mean, after being told the axes of lat and lon.
def global_mean(variable,lats,index_lat,index_lon):
    variable_zonal = np.nanmean(variable,axis=index_lon)
    if index_lon < index_lat: index_lat = index_lat-1
    lat_weights = np.cos(np.radians(lats))
    variable_global = np.average(variable_zonal,axis=index_lat,weights=lat_weights)
    return variable_global

# This function takes a time-lat-lon variable and computes the mean for a given range of lon and lat.
def spatial_mean(variable,lats,lons,lat_min,lat_max,lon_min,lon_max,index_lat,index_lon,verbose=False):
    #
    j_selected = np.where((lats >= lat_min) & (lats <= lat_max))[0]
    i_selected = np.where((lons >= lon_min) & (lons <= lon_max))[0]
    if verbose: print('Computing spatial mean. lats='+str(lats[j_selected[0]])+'-'+str(lats[j_selected[-1]])+', lons='+str(lons[i_selected[0]])+'-'+str(lons[i_selected[-1]])+'.  Points are inclusive.')
    #
    if   index_lon == 1: variable_zonal = np.nanmean(variable[:,i_selected],    axis=1)
    elif index_lon == 2: variable_zonal = np.nanmean(variable[:,:,i_selected],  axis=2)
    elif index_lon == 3: variable_zonal = np.nanmean(variable[:,:,:,i_selected],axis=3)
    else: print('Invalid lon dimension chosen'); return None
    #
    lat_weights = np.cos(np.radians(lats))
    if index_lon < index_lat: index_lat = index_lat-1
    if   index_lat == 0: variable_mean = np.average(variable_zonal[j_selected],    axis=0,weights=lat_weights[j_selected])
    elif index_lat == 1: variable_mean = np.average(variable_zonal[:,j_selected],  axis=1,weights=lat_weights[j_selected])
    elif index_lat == 2: variable_mean = np.average(variable_zonal[:,:,j_selected],axis=2,weights=lat_weights[j_selected])
    else: print('Invalid lat dimension chosen'); return None
    #
    return variable_mean

# A function to calculate magnitude and direction for a time series set of U and V values
def calc_magnitude_and_direction(u_component,v_component,ages,age_range):
    #
    # Get data in the selected ages
    ind_ages = np.where((ages > age_range[1]) & (ages <= age_range[0]))[0]
    u_selected = u_component[ind_ages]
    v_selected = v_component[ind_ages]
    #
    # Calculate magnitude and direction (i.e., direction that the wind blows from, with north being zero)
    magnitude = np.sqrt(np.square(u_selected) + np.square(v_selected))
    direction = np.zeros((len(u_selected))); direction[:] = np.nan
    for i in range(len(u_selected)):
        direction[i] = math.degrees(math.atan2(-v_selected[i],-u_selected[i]))
    #
    # Make sure that 0 is north and the angle goes clockwise from there.
    direction = direction - 90
    direction[direction < 0] = 360 + direction[direction < 0]
    direction = 360 - direction
    #
    return magnitude,direction