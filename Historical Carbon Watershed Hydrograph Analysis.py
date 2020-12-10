# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.5.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Historical Carbon Watershed Hydrology Analysis

# +
from landlab.io import read_esri_ascii #, write_esri_ascii
from landlab.plot.imshow import imshow_grid,imshow_grid_at_node
import rioxarray
import os
import matplotlib
import numpy as np
import matplotlib.pyplot as plt

# import geopandas as gpd
# from shapely.geometry import mapping
# import matplotlib.pyplot as plt
# import landlab
# 
from landlab.components import FlowAccumulator, KinwaveImplicitOverlandFlow, SinkFiller, SinkFillerBarnes
# from landlab.components import KinwaveImplicitOverlandFlow, SinkFiller
# from landlab.utils import get_watershed_masks_with_area_threshold
# from landlab.utils import get_watershed_masks

# import hsfm

# import geopandas as gpd
# from shapely.geometry import mapping
# import seaborn as sns

# import hvplot

# from pygeotools.lib import warplib
# import gdal

# import rasterio


# -

matplotlib.rcParams['figure.figsize'] = (10.0, 8)

base_data_dir = '/data2/elilouis/'
data_dir = f'{base_data_dir}/landlab-rainier/data'

# ## Convert SfM-generated DEMs to low res ASC files
# * Use rioxarray

dem_files = [
    f'{base_data_dir}rainier_carbon/input_data/73V3/00/00/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-2.30_y+5.33_z+1.14_align.tif',
    f'{base_data_dir}rainier_carbon/input_data/79V5/10/06/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-0.73_y+1.32_z+0.23_align.tif',
    f'{base_data_dir}rainier_carbon/input_data/90V3/09/13/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-8.89_y+8.19_z+0.19_align.tif',
    f'{base_data_dir}rainier_carbon/input_data/91V3/09/09/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x+1.44_y+7.60_z+3.60_align.tif',
    f'{base_data_dir}rainier_carbon/input_data/92V3/07/28/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-0.25_y+4.56_z+0.30_align.tif',
    f'{base_data_dir}rainier_carbon/input_data/92V5/10/06/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-1.71_y+6.27_z-0.44_align.tif',
]
reference_dem_file = f'{data_dir}landlab-rainier/data/carbon_watershed_dem_30m.tif'

output_res = 30.0

for file in dem_files:
    rds = rioxarray.open_rasterio(file)
    rds = rds.rio.reproject(rds.rio.crs, resolution=output_res)
    rds = rds.where(rds != rds.attrs['_FillValue'])
    rds = rds.fillna(-9999)
    rds = rds.astype('int16')
    rds.attrs.update({
        '_FillValue': -9999
    })
    year = file.split('input_data')[1].split('/')[1]
    output_path = f'{data_dir}/{year}.asc'
    print(f'Saving {output_path}')
    rds.rio.to_raster(
        output_path,
        driver='AAIGrid',
    )

# ## Landlab It

dem_path = os.path.join(data_dir, '73V3.asc')

(rmg, z) = read_esri_ascii(dem_path, name='topographic__elevation')

rmg.set_nodata_nodes_to_closed(rmg.at_node['topographic__elevation'], -9999)

imshow_grid(rmg, 'topographic__elevation')

sfb = SinkFillerBarnes(
    rmg
)

sfb.run_one_step()

# Set watershed boundary condition - use lowest elevation if basic call fails

rmg.set_watershed_boundary_condition(z)

outlet_node = rmg.set_watershed_boundary_condition(z, return_outlet_id=True, remove_disconnected=True)
outlet_node

fa = FlowAccumulator(rmg,
                     surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
#                      runoff_rate=None,
                     depression_finder='DepressionFinderAndRouter',
#                      depression_finder=None
#                      routing='D4'
                    )

fa.run_one_step()
(da, q) = fa.accumulate_flow()

outlet_node_to_sample = np.argmax(rmg.at_node['drainage_area'])
print('Outlet Node = ' + str(outlet_node_to_sample) + '; Drainage Area= ' + 
  str(da[outlet_node_to_sample] / 1000000) + ' km^2; Elev = '+ str(round(z[outlet_node_to_sample], 1)) + ' m')

# outlet-by-drainage-area, outlet-by-lowest-elevation

outlet_node_to_sample, outlet_id

imshow_grid_at_node(rmg, 'drainage_area')

imshow_grid_at_node(rmg, 'topographic__elevation')

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 6))
ax.xaxis.set_visible(False)
ax.set_facecolor("blue")
imshow_grid(rmg, z, plot_name='Spring Creek', var_name='topographic__elevation', var_units='m', grid_units=('m', 'm'), 
          cmap='terrain', color_for_closed='white')
ax.plot(rmg.node_x[outlet_node_to_sample], rmg.node_y[outlet_node_to_sample], 'ro', label='outlet')
_ = ax.legend(loc='lower right')

# +
starting_precip_mmhr = 50    # [mm/hour] this is runoff rate uniformly distributed in space
storm_duration = 1.0         # [hour]
model_run_time = 4.0         # [hour]

n = 0.1                     # Manning's roughness coefficient, (s/m^(1/3))

dt = 300                     # time step [sec] 

#Converting units to SI [m] and [sec]
storm_duration_sec=storm_duration*3600 # [sec]
model_run_time_sec=model_run_time*3600 # [sec]
# -

discharge_at_outlet = []
hydrograph_time = []
rmg.add_zeros("surface_water__depth", at="node", clobber=True)

kw = KinwaveImplicitOverlandFlow(rmg, runoff_rate=0.0, roughness=n, depth_exp=5/3)

# ### Run/Loop Kinematic Wave Overland Flow Component

# +
elapsed_time=1;

while elapsed_time < model_run_time_sec:  
    
    if elapsed_time < storm_duration_sec:
        kw.runoff_rate=starting_precip_mmhr #This needs to be in mm/hr because the source code automatically converts to m/s
    else:
        kw.runoff_rate=1e-20                #This needs to be > 0 because of an assertion in the source code... 
        
    kw.run_one_step(dt) 
    
    # add elapsed time to the continuos time
    hydrograph_time.append(elapsed_time/3600)    
    
    # pick discharge values from nodes identified earlier and store them for plotting
    discharge_at_outlet.append(rmg.at_node['surface_water_inflow__discharge'][outlet_node_to_sample])
     
    ## output time every now and then so that you know the code
    ## is actually running
    if (elapsed_time % 100) < 2:
        print("elapsed time = ", elapsed_time)
        
    ## Updating elapsed_time  
    elapsed_time += dt

# +
plt.figure(1)
plt.plot(hydrograph_time, discharge_at_outlet, "r-", label="outlet")

plt.ylabel("Discharge (cms)")
plt.xlabel("Time (hour)")
plt.legend(loc="upper right")
title_text = "Hydrographs at three locations"
plt.title(title_text)
plt.show()
# -

# We can plot some **snapshots of water depth** on the domain. There is room to improve this plot.

# +
elapsed_time = 1.
run_time_slices = (1.2*3600, 1.25*3600)
for t in run_time_slices:
    while elapsed_time < t:
         # First, we calculate our time step.
        kw.run_one_step(dt) 

#         # Increased elapsed time
        elapsed_time += dt 
#     figure(t)
    imshow_grid(rmg, 'surface_water__depth', cmap='Blues')
# -

imshow_grid(rmg,'surface_water__depth', plot_name = 'Surface water depth', 
            var_name = 'Depth of water', var_units = 'm', grid_units = ('m','m'), 
            cmap = 'jet', limits = (0, 1.6))


