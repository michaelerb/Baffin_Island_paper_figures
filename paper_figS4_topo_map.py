#==============================================================================
# This script visualizes some TraCE-21ka variables near Lake CF8 (on Baffin
# Island). This is Fig. S4 in the paper.
#    author: Michael Erb
#    date  : 2/7/2023
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import xarray as xr
import copy

save_instead_of_plot = True


#%% LOAD DATA

# Load the decadal-mean TraCE-21ka data
data_dir = '/projects/pd_lab/data/models/TraCE_21ka/TraCE/'
data_xarray = xr.open_dataset(data_dir+'trace.01-36.22000BP.cam2.PHIS.22000BP_decavg_400BCE.nc',decode_times=False)
lat  = data_xarray['lat'].values
lon  = data_xarray['lon'].values
year = np.around(-1000*data_xarray['time'].values)


#%% CALCULATIONS

# Convert units
geopotential_height = data_xarray['PHIS']/9.80665  # Calculate geopotential height (see https://glossary.ametsoc.org/wiki/Geopotential_height)

# Set locations
lake_lat       = 70.55
lake_lon       = -68.95 + 360
greenland_lat  = 72.36
greenland_lon  = -41.25 + 360
laurentide_lat = 61.2
laurentide_lon = -101.25 + 360

# Find data in the interval 12 - 7 ka
ind_time = np.where((year >= 7000) & (year <= 12000))[0]

# Find the time closest to 12 ka
ind_12ka = np.argmin(np.abs(12000 - year))

# Print out the elevations of the chosen locations
j_greenland  = np.argmin(np.abs(lat - greenland_lat))
i_greenland  = np.argmin(np.abs(lon - greenland_lon))
j_laurentide = np.argmin(np.abs(lat - laurentide_lat))
i_laurentide = np.argmin(np.abs(lon - laurentide_lon))
print('Greenland location, geopotential height (m):',geopotential_height[ind_12ka,j_greenland,i_greenland].values)
print('Laurentide location, geopotential height (m):',geopotential_height[ind_12ka,j_laurentide,i_laurentide].values)


#%%
# Compute the vertices of the trace grid, for plotting purposes
lon_wrap = copy.deepcopy(lon)
lon_wrap = np.insert(lon_wrap,0,lon_wrap[-1]-360)  # Add the right-most lon point to the left
lon_wrap = np.append(lon_wrap,lon_wrap[1]+360)     # Add the left-most lon point to the right
lon_edges = (lon_wrap[:-1] + lon_wrap[1:])/2

lat_edges = copy.deepcopy(lat)
lat_edges = (lat[:-1] + lat[1:])/2
lat_edges = np.insert(lat_edges,0,-90)  # Add the South Pole to the beginning
lat_edges = np.append(lat_edges,90)     # Add the North Pole to the end


#%% FIGURES
plt.style.use('ggplot')

# Make a map
plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1),(0,0),projection=ccrs.LambertConformal()); ax1.set_extent([-135,-35,45,82],ccrs.PlateCarree())
map1 = ax1.pcolormesh(lon_edges,lat_edges,geopotential_height[ind_12ka,:,:],cmap='Reds',vmin=0,vmax=2500,transform=ccrs.PlateCarree())
plt.colorbar(map1,orientation='horizontal',label='Geopotenital height (m)',ax=ax1,fraction=0.08,pad=0.02)
ax1.scatter(lake_lon,lake_lat,500,c='yellow',edgecolors='k',marker='*',alpha=1,transform=ccrs.PlateCarree())
ax1.scatter(greenland_lon, greenland_lat, 50,c='gray',marker='o',alpha=1,transform=ccrs.PlateCarree())
ax1.scatter(laurentide_lon,laurentide_lat,50,c='gray',marker='o',alpha=1,transform=ccrs.PlateCarree())
ax1.coastlines()
ax1.set_title('Annual surface height at 12 ka, TraCE-21ka',fontsize=16,loc='center')

if save_instead_of_plot:
    plt.savefig('figS4_trace_topo_12ka.png',dpi=300,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()

