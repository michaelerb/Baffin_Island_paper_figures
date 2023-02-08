#==============================================================================
# This script visualizes some TraCE-21ka variables near Lake CF8 (on Baffin
# Island). This is Fig. S11 in the paper.
#    author: Michael Erb
#    date  : 2/7/2023
#==============================================================================

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.util as cutil
import xarray as xr
import utils
from windrose import WindroseAxes

save_instead_of_plot = True

# Select the season to plot
season = 'JJA'
#season = 'DJF'
#season = 'Annual'


#%% LOAD DATA

# Note: TraCE-21ka data can downloaded from https://www.earthsystemgrid.org/project/trace.html
data_dir = '/projects/pd_lab/data/models/TraCE_21ka/TraCE/'
if   season == 'Annual': season_txt = ''    # Annual
elif season == 'JJA':    season_txt = 'JJA' # JJA
elif season == 'DJF':    season_txt = 'DJF' # DJF
files_to_load = [data_dir+'trace.01-36.22000BP.cam2.PSL.22000BP_decavg'+season_txt+'_400BCE.nc',
                 data_dir+'trace.01-36.22000BP.cam2.PHIS.22000BP_decavg'+season_txt+'_400BCE.nc',
                 data_dir+'trace.01-36.22000BP.cam2.Q.22000BP_decavg'+season_txt+'_400BCE.nc',
                 data_dir+'trace.01-36.22000BP.cam2.U.22000BP_decavg'+season_txt+'_400BCE.nc',
                 data_dir+'trace.01-36.22000BP.cam2.V.22000BP_decavg'+season_txt+'_400BCE.nc']

# Load the decadal-mean TraCE-21ka data
data_xarray = xr.open_mfdataset(files_to_load,decode_times=False)
age = np.around(-1000*data_xarray['time'].values)
lat = data_xarray['lat'].values
lon = data_xarray['lon'].values
lev = data_xarray['lev'].values

### Calculate pressure at every location
PSL   = data_xarray['PSL'].values
P0    = data_xarray['P0'].values
A_mid = data_xarray['hyam'].values
B_mid = data_xarray['hybm'].values
A_interfaces = data_xarray['hyai'].values
B_interfaces = data_xarray['hybi'].values

pressure_mid        = (A_mid        * P0)[None,:,None,None] + (B_mid[None,:,None,None]        * PSL[:,None,:,:])
pressure_interfaces = (A_interfaces * P0)[None,:,None,None] + (B_interfaces[None,:,None,None] * PSL[:,None,:,:])

pressure_mid_globalmean = np.mean(utils.global_mean(pressure_mid,data_xarray['lat'].values,2,3),axis=0)
print(pressure_mid_globalmean)

### Calculate moisture tranport
# Checking units for the moisture flux calculation
# U,V:     m/s
# Q:       kg/kg
# ilev:    hybrid (~Pa)
# gravity: m/s2
# Moisture flux results: kg m-1 s-1

# Calculate vertically-integrated moisture transport (see https://psl.noaa.gov/data/composites/day/calculation.html)
QU = data_xarray['Q'].values * data_xarray['U'].values
QV = data_xarray['Q'].values * data_xarray['V'].values

# Vertically integrate
ind_selected = np.where(pressure_mid_globalmean >= 30000)[0]
pressure_diff = pressure_interfaces[:,1:,:,:] - pressure_interfaces[:,:-1,:,:]
QU_int = np.sum(QU[:,ind_selected,:,:] * pressure_diff[:,ind_selected,:,:],axis=1) / 9.81
QV_int = np.sum(QV[:,ind_selected,:,:] * pressure_diff[:,ind_selected,:,:],axis=1) / 9.81

# Convert units
data_xarray['PSL']  = data_xarray['PSL']/100       # Convert pressure from Pa to hPa
geopotential_height = data_xarray['PHIS']/9.80665  # Calculate geopotential height (see https://glossary.ametsoc.org/wiki/Geopotential_height)


#%% CALCULATIONS

# Find data in the interval 12 - 7 ka
ind_time = np.where((age >= 7000) & (age <= 12000))[0]

# The lake is between four model gridpoints. Average over those four gridpoints
lake_lat = 70.55
lake_lon = -68.95 + 360
lat_range = [67,74]
lon_range = [287,295]
PSL_mean                 = utils.spatial_mean(data_xarray['PSL'],   lat,lon,lat_range[0],lat_range[1],lon_range[0],lon_range[1],1,2)
geopotential_height_mean = utils.spatial_mean(geopotential_height,  lat,lon,lat_range[0],lat_range[1],lon_range[0],lon_range[1],1,2)
QU_int_mean              = utils.spatial_mean(QU_int,               lat,lon,lat_range[0],lat_range[1],lon_range[0],lon_range[1],1,2)
QV_int_mean              = utils.spatial_mean(QV_int,               lat,lon,lat_range[0],lat_range[1],lon_range[0],lon_range[1],1,2)

# Compute a SLP index between Greenland and the Laurentide
greenland_lat  = 72.36
greenland_lon  = -41.25 + 360
laurentide_lat = 61.2
laurentide_lon = -101.25 + 360
j_greenland  = np.argmin(np.abs(lat - greenland_lat))
i_greenland  = np.argmin(np.abs(lon - greenland_lon))
j_laurentide = np.argmin(np.abs(lat - laurentide_lat))
i_laurentide = np.argmin(np.abs(lon - laurentide_lon))
slp_index = data_xarray['PSL'][:,j_greenland,i_greenland] - data_xarray['PSL'][:,j_laurentide,i_laurentide]

# Compute quantities during time periods
ind_period1 = np.where((age > 11000) & (age <= 12000))[0]
ind_period2 = np.where((age > 7000)  & (age <= 8000))[0]
QU_period1 = np.mean(QU_int[ind_period1,:,:],axis=0)
QV_period1 = np.mean(QV_int[ind_period1,:,:],axis=0)
QU_period2 = np.mean(QU_int[ind_period2,:,:],axis=0)
QV_period2 = np.mean(QV_int[ind_period2,:,:],axis=0)
SLP_period1 = np.mean(data_xarray['PSL'][ind_period1,:,:],axis=0)
SLP_period2 = np.mean(data_xarray['PSL'][ind_period2,:,:],axis=0)


#%% FIGURES
plt.style.use('ggplot')

# Plot SLP index and moisture flux wind roses
plt.figure(figsize=(18,9))
ax_top = plt.subplot2grid((2,5),(0,0),colspan=5)
ax_bottom = {}

ax_top.plot(age,slp_index,linewidth=3)
ax_top.axhline(y=0,color='k',linewidth=1,linestyle='dashed',alpha=0.5)
ax_top.set_title('SLP index',fontsize=18,loc='left')
ax_top.set_ylabel('SLP index (hPa)',fontsize=16)
ax_top.set_xlim(12000,7000)
ax_top.set_xlabel('Age (yr BP)',fontsize=16)

age_change = 1000
for i,age_old in enumerate(np.arange(12000,7000,-age_change)):
    #
    ax_bottom[i] = plt.subplot2grid((2,5),(1,i),projection='windrose')
    age_new = age_old - age_change
    title_txt = 'Q-flux, '+str(int(age_old/1000))+'-'+str(int(age_new/1000))+' ka'
    bins = np.arange(0,20.1,5)
    wind_speeds,wind_directions = utils.calc_magnitude_and_direction(QU_int_mean,QV_int_mean,age,[age_old,age_new])
    ax_bottom[i].bar(wind_directions,wind_speeds,bins=bins,normed=True,opening=0.8,edgecolor='white')
    if i == 0: ax_bottom[i].set_legend()
    ax_bottom[i].set_xticklabels(['E','N-E','N','N-W','W','S-W','S','S-E'])
    ax_bottom[i].set_title(title_txt,fontsize=18,loc='left')

plt.tight_layout()
plt.suptitle(season+'-mean sea level pressure index and moisture flux (kg m$^{-1}$ s$^{-1}$) at Lake CF8',fontsize=26)
plt.subplots_adjust(top=.88)

if save_instead_of_plot:
    plt.savefig('extrafig_trace_slp_index_and_qflux_'+season+'.png',dpi=300,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()


#%% Make a two panel figure of moisture flux

plt.figure(figsize=(13,8))
ax1 = plt.subplot2grid((1,2),(0,0),projection=ccrs.LambertConformal(central_longitude=lake_lon)); ax1.set_extent([-95,-45,60,80],ccrs.PlateCarree())
ax2 = plt.subplot2grid((1,2),(0,1),projection=ccrs.LambertConformal(central_longitude=lake_lon)); ax2.set_extent([-95,-45,60,80],ccrs.PlateCarree())

psl_cyclic,lon_cyclic = cutil.add_cyclic_point(SLP_period1,coord=lon)
map1 = ax1.contourf(lon_cyclic,lat,psl_cyclic,np.arange(1010,1028,1),extend='both',transform=ccrs.PlateCarree())
plt.colorbar(map1,orientation='horizontal',label='Sea level pressure (hPa)',ax=ax1,fraction=0.08,pad=0.02)
wind_period1 = ax1.quiver(lon,lat,QU_period1,QV_period1,scale=200,color='k',width=.003,zorder=2,transform=ccrs.PlateCarree())
ax1.quiverkey(wind_period1,X=.9,Y=1.05,labelpos='S',U=20,label='20 kg m$^{-1}$ s$^{-1}$')
gl = ax1.gridlines(color='k',linewidth=1,linestyle=(0,(1,5)))
ax1.set_title('12-11 ka',fontsize=18,loc='center')

psl_cyclic,lon_cyclic = cutil.add_cyclic_point(SLP_period2,coord=lon)
map2 = ax2.contourf(lon_cyclic,lat,psl_cyclic,np.arange(1010,1028,1),extend='both',transform=ccrs.PlateCarree())
plt.colorbar(map2,orientation='horizontal',label='Sea level pressure (hPa)',ax=ax2,fraction=0.08,pad=0.02)
wind_period2 = ax2.quiver(lon,lat,QU_period2,QV_period2,scale=200,color='k',width=.003,zorder=2,transform=ccrs.PlateCarree())
ax2.quiverkey(wind_period2,X=.9,Y=1.05,labelpos='S',U=20,label='20 kg m$^{-1}$ s$^{-1}$')
gl = ax2.gridlines(color='k',linewidth=1,linestyle=(0,(1,5)))
ax2.set_title('8-7 ka',fontsize=18,loc='center')

for ax in [ax1,ax2]:
    ax.scatter(lake_lon,lake_lat,500,c='yellow',edgecolors='k',marker='*',alpha=1,transform=ccrs.PlateCarree())
    ax.coastlines(linewidth=1,color='lightgray',zorder=1)

plt.tight_layout()
plt.suptitle(season+'-mean sea level pressure and moisture flux in TraCE-21ka',fontsize=22)
plt.subplots_adjust(top=1)

if save_instead_of_plot:
    plt.savefig('figS11_trace_map_qflux_'+season+'.png',dpi=300,format='png',bbox_inches='tight')
    plt.close()
else:
    plt.show()
