#Load library
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.ticker import ScalarFormatter

# Load the datasets
sensitivity_cooling = xr.open_dataset("sensitivity_map_computed_c.nc", decode_times=False)
sensitivity_warming = xr.open_dataset("sensitivity_map_computed_w.nc", decode_times=False)
warming_cooling=(sensitivity_cooling+sensitivity_warming)/2
# Set up the figure and subplots (1 row, 2 columns)
fig, axes = plt.subplots(1, 3, figsize=(16, 8), subplot_kw={'projection': ccrs.PlateCarree(central_longitude=180)})

# List of datasets and titles for each plot
datasets = [sensitivity_warming*ice_free_mask*ocean_mask, sensitivity_cooling*ice_free_mask*ocean_mask, warming_cooling*ice_free_mask*ocean_mask]  ### with necessary mask
titles = [r'Warming ', r'Cooling ', r'Warming & Cooling']

# Define the min and max for color normalization (same for both subplots)
vmin = -30 # Set the minimum value
vmax = 30   # Set the maximum value

# Loop over the datasets and plot
for i, (data, title) in enumerate(zip(datasets, titles)):
    ax = axes[i]  # Select the current axis
    ax.coastlines(linewidth=0.5, color='black')
    ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.5)  # Draw labels only on the first plot
    #ax.add_feature(cfeature.LAND, color='white', zorder=10)  # Mask land areas
    # Plot the data with the same color range for both subplots
    plot = data.dR_dSST_global.plot(
        ax=ax, 
        transform=ccrs.PlateCarree(),  # Data projection
        cmap='seismic',                # Choose a colormap
        add_colorbar=False,            # Do not add individual colorbars
        vmin=vmin, vmax=vmax           # Use the same vmin and vmax for both
    )

    # Set the title for each subplot
    ax.set_title(title, fontsize=16)

# Create a single colorbar at the bottom for both plots
cbar = fig.colorbar(
    plot, 
    ax=axes,                      
    orientation='horizontal', 
    fraction=0.05,                
    pad=0.15,                     
    aspect=40,
    extend='both'  # Adds sharp, clear ends to the colorbar

)

# Set colorbar label
cbar.set_label(r'$\frac{\partial R_i}{\partial SST_l}$ [${Wm^2K^{-1}}$]', fontsize=14)

# Format the colorbar numbers in scientific notation
cbar.formatter = ScalarFormatter(useMathText=True)
#cbar.formatter.set_powerlimits((-4, -4))  # Force e-4 scaling
cbar.update_ticks()
# Adjust layout to prevent overlap and add space for colorbar
plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.20, wspace=0.1)  # Adjust layout manually

# Show the plot
plt.show()
