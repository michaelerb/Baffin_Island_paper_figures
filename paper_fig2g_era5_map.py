#==============================================================================
# This script visualizes some ERA variables near Lake CF8 (on Baffin Island).
# This is Fig. 2g in the paper.
#    author: Michael Erb
#    date  : 2/7/2023
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.util as cutil
import xarray as xr

save_instead_of_plot = True


#%% LOAD DATA

data_dir = '/projects/pd_lab/data/modern_datasets/ERA5/'
files_to_load = [data_dir+'era5_monthly_10m_u.nc',
                 data_dir+'era5_monthly_10m_v.nc',
                 data_dir+'era5_monthly_slp.nc']

# Load the ERA5 data
data_xarray = xr.open_mfdataset(files_to_load)
lat  = data_xarray.latitude.values
lon  = data_xarray.longitude.values
time = data_xarray.time.values
# Note: Check what months are included in the data file, then decide whether to keep the last year.


#%% CALCULATIONS

# Convert units
data_xarray.msl.values = data_xarray.msl.values/100  # Convert pressure from Pa to hPa

# Compute seasonal means
month_lengths = data_xarray.time.dt.days_in_month
weights = (month_lengths.groupby("time.season") / month_lengths.groupby("time.season").sum())

# Test that the sum of weights for each season is near 1
np.testing.assert_allclose(weights.groupby("time.season").sum().values, np.ones(4))

# Compute time-weighted seasonal values
data_xarray_seasonal = (data_xarray * weights).groupby("time.season").sum(dim="time")


#%% Set lake location
lake_lat = 70.55
lake_lon = -68.95 + 360


#%% FIGURES
plt.style.use('ggplot')

# Make a map for each season
season_txt = data_xarray_seasonal.season.values

i = 1
print('Making seasonal map:',season_txt[i])
slp_season = data_xarray_seasonal.msl[i,0,:,:].values
u10_season = data_xarray_seasonal.u10[i,0,:,:].values
v10_season = data_xarray_seasonal.v10[i,0,:,:].values
range_chosen = np.arange(1010,1028,1)
x_skip = 20; y_skip = 10; vector_scale = 50

plt.figure(figsize=(10,10))
ax1 = plt.subplot2grid((1,1),(0,0),projection=ccrs.LambertConformal()); ax1.set_extent([-135,-35,45,82],ccrs.PlateCarree())

slp_cyclic,lon_cyclic = cutil.add_cyclic_point(slp_season,coord=lon)
map1 = ax1.contourf(lon_cyclic,lat,slp_cyclic,range_chosen,extend='both',transform=ccrs.PlateCarree())
wind_ann = ax1.quiver(lon[::x_skip],lat[::y_skip],u10_season[::y_skip,::x_skip],v10_season[::y_skip,::x_skip],scale=vector_scale,color='k',linewidth=3,zorder=2,transform=ccrs.PlateCarree())
ax1.quiverkey(wind_ann,X=.9,Y=1.05,labelpos='S',U=10,label='Wind speed, 10 m/s')
ax1.coastlines(linewidth=1,color='lightgray',zorder=1)
ax1.scatter(lake_lon,lake_lat,500,c='yellow',edgecolors='k',marker='*',alpha=1,transform=ccrs.PlateCarree())
plt.colorbar(map1,orientation='horizontal',label='Sea level pressure (hPa)',ax=ax1,fraction=0.08,pad=0.02)
ax1.set_title('ERA5 reanalysis, mean of 1959-2022 CE\n'+season_txt[i]+' SLP (shaded) and 10m wind (vectors)',fontsize=12,loc='center')

if save_instead_of_plot:
    plt.savefig('fig2g_era5_map_'+season_txt[i]+'.png',dpi=300,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()
