####Load library
import numpy as np 
import pandas as pd 
import xarray as xr 
import cartopy 
import glob
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import os


                                        #####Loading and arranging######
#load sst anomaly data and rename and taking time mean for the conveince 
sst_anomaly=xr.open_mfdataset("/glade/campaign/univ/wyom0177/E3SM/patches_anomaly.nc", decode_times=False)
sst_anomaly = sst_anomaly.rename({'latitude': 'lat', 'longitude': 'lon'})
sst_anomaly= sst_anomaly.rename({var: var.replace('_anomaly', '') for var in sst_anomaly.data_vars if var.endswith('_anomaly')})
sst_anomaly=sst_anomaly.mean('time')

#Load control
control=xr.open_mfdataset("/glade/campaign/univ/wyom0177/E3SM/regrided/F2010.control/F2010.control.nc", decode_times=False)


# Load and aarange the patch experiments
main_dir = "/glade/campaign/uwyo/wyom0167/E3SM/regrided/"  # Update with your main directory path

# List to store matching .nc files
sst_files = []

# Dictionary to store data from each .nc file
data_dict = {}

# Traverse through all subdirectories in the main directory
for subdir, _, files in os.walk(main_dir):
    # Skip the subdirectory 'F2010.control'
    if 'F2010.control' in subdir:
        continue
    
    # Get the name of the subdirectory (without the full path)
    subdir_name = os.path.basename(subdir)
    
    # Find the .nc file that exactly matches the subdirectory's name
    matching_files = [f for f in files if f == subdir_name + '.nc']
    
    # Append the matching files to the list and read using xarray
    for file in matching_files:
        file_path = os.path.join(subdir, file)
        sst_files.append(file_path)
        
        # Open the netCDF file using xarray and store the dataset in the dictionary
        data_dict[file] = xr.open_mfdataset(file_path , decode_times=False, decode_cf=False )
updated_data_dict = {}
for key, value in data_dict.items():
    new_key = key.replace('F2010.', 'sst_patch_').replace('.nc', '').replace('n', '-')
    updated_data_dict[new_key] = value
data_dict=updated_data_dict



                                                ##calculation

# Initialize global summation arrays (zeros matching the grid shape)
sum_weighted_response = xr.zeros_like(control['FSNT'].mean('time'))
sum_sst_anomaly = xr.zeros_like(sum_weighted_response)

# Extract patch names from the sst_anomaly_list (adjust according to your naming convention)

patch_names = [var_name for var_name in sst_anomaly.data_vars]

# Loop through all patch names from the sst_anomaly_list
for patch_name in patch_names:
    # Check if this patch name exists in both sst_anomaly_list and data_dict
    if patch_name in sst_anomaly.data_vars and patch_name in data_dict:
        # Extract SST anomaly data for this patch (assuming sst_anomaly is a dictionary)
        sst_anomaly_data = sst_anomaly[patch_name]
        print(patch_name)
        # Extract TOA flux data for this patch and take weighted average
        ds = data_dict[patch_name]
        area_weighted_FSNT = (ds['FSNT'].mean('time') * control['area']).sum(dim=['lat', 'lon']) / control['area'].sum(dim=['lat', 'lon'])
        area_weighted_FLNT = (ds['FLNT'].mean('time') * control['area']).sum(dim=['lat', 'lon']) / control['area'].sum(dim=['lat', 'lon'])

        control_area_weighted_FSNT = (control['FSNT'].mean('time') * control['area']).sum(dim=['lat', 'lon']) / control['area'].sum(dim=['lat', 'lon'])
        control_area_weighted_FLNT = (control['FLNT'].mean('time') * control['area']).sum(dim=['lat', 'lon']) / control['area'].sum(dim=['lat', 'lon'])

        # Difference of area-weighted means
        delta_R_i = (area_weighted_FSNT - area_weighted_FLNT) - (control_area_weighted_FSNT - control_area_weighted_FLNT)
        
        # Calculate grid area (assuming 'control.area' exists)
        R = 6371000  # radius of Earth in meters
        a_j = control['area'] * R**2  # Convert solid angle to surface area (m^2)

        # Create binary patch mask for positive SST anomalies (mask for inside patch)
        patch_mask = (sst_anomaly_data < -0.05)   ##### change here for warming and cooling 
        ocean_mask =control.OCNFRAC.mean(dim="time")>=0.5 
        ice_free_mask = control.ICEFRAC.mean(dim="time")<=0.05 

       # Total patch area (sum of areas inside the patch) & area fraction
        grid_area_patch_with_ice = a_j * patch_mask*ice_free_mask

        # Calculate weighted SST anomaly globally 
        # Area-weighted SST anomaly (prevent division by zero)
        weighted_sst = sst_anomaly_data * grid_area_patch_with_ice 
        delta_SST_p = weighted_sst.sum() / (a_j).sum()   # scalar value for the patch

        if delta_SST_p!=0:
            partial_derivative = (delta_R_i / delta_SST_p)
        # Initialize the global summation arrays
            sum_weighted_response += partial_derivative * sst_anomaly_data*ice_free_mask   ###making ice free SST anomaly
            sum_sst_anomaly += sst_anomaly_data

# Final sensitivity map (ensure global normalization)
sensitivity_map = sum_weighted_response / sum_sst_anomaly.where(sum_sst_anomaly != 0, 1)  #### if in case there is any zero value in sum_sst_anomaly
#sensitivity_map = sensitivity_map.fillna(0) 

sensitivity_map_computed = sensitivity_map.compute()
#partial_derivative = partial_derivative.compute()

# Save sensitivity_map_computed as a NetCDF file
# Set the name of the DataArray
sensitivity_map_computed.name = "dR_dSST_global"
# Save it to a NetCDF file
sensitivity_map_computed.to_netcdf('sensitivity_map_computed_c.nc')