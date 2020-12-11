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
# from shapely.geometry import mapping
# import matplotlib.pyplot as plt
# import landlab
# 
from landlab.components import FlowAccumulator, OverlandFlow, \
    KinwaveImplicitOverlandFlow, SinkFiller, SinkFillerBarnes
# from landlab.components import KinwaveImplicitOverlandFlow, SinkFiller
# from landlab.utils import get_watershed_masks_with_area_threshold
# from landlab.utils import get_watershed_masks

# import hsfm

# import seaborn as sns

# import hvplot

# from pygeotools.lib import warplib
# import gdal

import rasterio


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

# ## Get common-extent polygon

from pygeotools.lib import geolib
import gdal
import shapely.wkt
import geopandas as gpd

dem_file_list = [
    '73V3.asc',
    '79V5.asc',
    '90V3.asc',
    '91V3.asc',
    '92V3.asc',
    '92V5.asc'
]


def get_polygon_from_raster_file(file_path, simplify_tolerance = 30):
    ds = gdal.Open(file_path)
    outline_ogr_polygon = geolib.get_outline(ds)
    wkt = outline_ogr_polygon.ExportToWkt()
    polygon = shapely.wkt.loads(wkt)
    polygon = polygon.simplify(simplify_tolerance) # simplify to the same tolerance as the DEMs resolution
    return polygon


polygon_list = [
    get_polygon_from_raster_file(
        os.path.join(data_dir, dem_file)
    ) 
    for dem_file in dem_file_list
]
gdf = gpd.GeoDataFrame(
    polygon_list, 
    geometry=0, 
    crs = gdal.Open(os.path.join(data_dir, dem_file_list[0])).GetProjection()
)

gdf.plot()

(rmg, z) = read_esri_ascii(dem_path, name='topographic__elevation')
imshow_grid(rmg, 'topographic__elevation')



# ## Landlab It

# 92V5 works!

dem_path = os.path.join(data_dir, '73V3.asc')

(rmg, z) = read_esri_ascii(dem_path, name='topographic__elevation')

rmg.set_nodata_nodes_to_closed(rmg.at_node['topographic__elevation'], -9999)

imshow_grid(rmg, 'topographic__elevation')

sfb = SinkFillerBarnes(
    rmg,
    fill_flat=False
)

sfb.run_one_step()

# Set watershed boundary condition - use lowest elevation if basic call fails

rmg.set_watershed_boundary_condition(z)

fa = FlowAccumulator(rmg,
#                      surface='topographic__elevation',
                     flow_director='FlowDirectorD8',
#                      runoff_rate=None,
#                      depression_finder='DepressionFinderAndRouter',
#                      routing='D4'
                    )

fa.run_one_step()
(da, q) = fa.accumulate_flow()

outlet_node_to_sample = np.argmax(rmg.at_node['drainage_area'])
print('Outlet Node = ' + str(outlet_node_to_sample) + '; Drainage Area= ' + 
  str(da[outlet_node_to_sample] / 1000000) + ' km^2; Elev = '+ str(round(z[outlet_node_to_sample], 1)) + ' m')

# For each of the links, there is a tail and head and node.
# Look at tail nodes of all links, and look at node area, find link that connects 

 #the link number that carries the largest q
outlet_link_to_sample = rmg.links_at_node[outlet_node_to_sample][3]

rmg.links_at_node[outlet_node_to_sample]

imshow_grid_at_node(rmg, 'drainage_area')

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

n = 0.02                     # Manning's roughness coefficient, (s/m^(1/3))

dt = 1000                     # time step [sec] 

#Converting units to SI [m] and [sec]
storm_duration_sec=storm_duration*3600 # [sec]
model_run_time_sec=model_run_time*3600 # [sec]
# -

discharge_at_outlet = []
hydrograph_time = []
rmg.add_zeros("surface_water__depth", at="node", clobber=True)

# ### Run/Loop Overland Flow Component

alpha = 0.15           # time-step factor (nondimensional; from Bates et al., 2010) range 0.15-0.7
starting_precip_ms = starting_precip_mmhr * (2.77778 * 10 ** -7)

of = OverlandFlow(rmg, alpha=alpha, mannings_n=n, steep_slopes=True)

# +
elapsed_time = 1.0  # s

while elapsed_time < model_run_time_sec:
    # First, we calculate our time step.
    of.dt = of.calc_time_step()
    # Now, we can generate overland flow.
     
    ## The storm starts when the model starts. While the elapsed time is less
    ## than the storm duration, we add water to the system as rainfall.
    if elapsed_time < (storm_duration_sec):
        of.rainfall_intensity = starting_precip_ms
    else:  # elapsed time exceeds the storm duration, rainfall ceases.
        of.rainfall_intensity = 0.0
   
            
    of.run_one_step()
    
    # add elapsed time to the continuos time
    hydrograph_time.append(elapsed_time/3600)
    
    # pass the surface water discharge to variable q. This includes unit q for all the links
    q = rmg.at_link["surface_water__discharge"]
    
    #q is surface water discharge - divide by x-sectional area for fluid depth
    # Search - does OverlandFlow provide depth to water
    # can also calculate shear stress
    
    #append q for each location we would like to plot hydrograph and calculate Q=q*w.
    discharge_at_outlet.append(np.abs(q[outlet_link_to_sample]) * rmg.dx)
    

     ## output time every now and then so that you know the code
    ## is actually running
    if (elapsed_time % 100) < 2:
        print("elapsed time = ", elapsed_time)
        
     ## Updating elapsed_time  
    elapsed_time += of.dt

# +
plt.figure(1)
plt.plot(hydrograph_time, discharge_at_outlet, "r-", label="outlet")

plt.ylabel("Discharge (cms)")
plt.xlabel("Time (hour)")
plt.legend(loc="upper right")
title_text = f"Hydrograph ({dem_path.split('/')[-1].split('.')[0]})"
plt.title(title_text)
plt.show()
# -

# ### Run/Loop Kinematic Wave Overland Flow Component

kw = KinwaveImplicitOverlandFlow(rmg, runoff_rate=0.0, roughness=n, depth_exp=5/3)

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


