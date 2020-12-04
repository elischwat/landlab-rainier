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

# +
import numpy as np

# import landlab
# from landlab.io import read_esri_ascii, write_esri_ascii 
# from landlab.components import FlowAccumulator
# from landlab.utils import get_watershed_masks_with_area_threshold
# from landlab.utils import get_watershed_masks
# from landlab.plot.imshow import imshow_grid

import hsfm
import rioxarray as rxr
import geopandas as gpd
from shapely.geometry import mapping
import seaborn as sns
import matplotlib.pyplot as plt
import hvplot
# -

# !find /Volumes/MyDrive/rainier_carbon/input_data/ -type f -name "*_align.tif"

files = [
    '/Volumes/MyDrive/rainier_carbon/input_data//92V3/07/28/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-0.25_y+4.56_z+0.30_align.tif',
    '/Volumes/MyDrive/rainier_carbon/input_data//94V6/09/16/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-589.73_y+267.08_z+43.33_align.tif',
    '/Volumes/MyDrive/rainier_carbon/input_data//73V3/00/00/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-2.30_y+5.33_z+1.14_align.tif',
    '/Volumes/MyDrive/rainier_carbon/input_data//92V5/10/06/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-1.71_y+6.27_z-0.44_align.tif',
    '/Volumes/MyDrive/rainier_carbon/input_data//87V1/08/21/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-9.73_y+27.82_z+4.66_align.tif',
    '/Volumes/MyDrive/rainier_carbon/input_data//90V3/09/13/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-8.89_y+8.19_z+0.19_align.tif',
    '/Volumes/MyDrive/rainier_carbon/input_data//91V3/09/09/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x+1.44_y+7.60_z+3.60_align.tif',
    '/Volumes/MyDrive/rainier_carbon/input_data//79V5/10/06/sfm/cluster_000/metashape0/pc_align/run-run-run-trans_source-DEM_dem_align/run-run-run-trans_source-DEM_reference_dem_clip_nuth_x-0.73_y+1.32_z+0.23_align.tif',
]

# +
# hsfm.plot.plot_dem_difference_from_file_name()
carbon_whole = rxr.open_rasterio(files[0], masked=True).squeeze()

carbon_whole.hvplot.image(
    cmap='RdBu'
).redim(
    value=dict(range=(-10, 10))
).opts(
    xaxis=None, 
    yaxis=None,
    width=400,
    height=500,
    show_frame=False,

)
# -

raster = rio.open_rasterio(files[0], masked=True).squeeze()

(rmg, elevations) = read_esri_ascii("data/carbon_watershed_dem_30m.asc", name="topographic__elevation")

rmg.set_nodata_nodes_to_closed(elevations, -3.4028234663852886e+38)

# +
outlet_id = rmg.core_nodes[np.argmin(rmg.at_node['topographic__elevation'][rmg.core_nodes])] # find the lowest point on our DEM?            

rmg.set_watershed_boundary_condition_outlet_id(outlet_id, elevations)    # set the lowest point as the outlet

print("Outlet ID=", outlet_id)                                        # print outlet id number
print("Outlet elevation=",rmg.at_node['topographic__elevation'][outlet_id])        # print elevation of outlet node
print("Min elevation of core nodes=", np.min(rmg.at_node['topographic__elevation'][rmg.core_nodes])) # print minimum elevation of core nodes
print("Max elevation of core nodes=", np.max(rmg.at_node['topographic__elevation'][rmg.core_nodes])) # print maximum elevation of core nodes
# -

rmg.status_at_node[outlet_id] == rmg.BC_NODE_IS_FIXED_VALUE

imshow_grid(rmg,'topographic__elevation')

imshow_grid(rmg, rmg.status_at_node, color_for_closed="blue")

fr = FlowAccumulator(
    rmg,
    surface='topographic__elevation',
    flow_director='FlowDirectorD8',
    runoff_rate=None,
    depression_finder='DepressionFinderAndRouter')

(drainage_area, discharge) = fr.accumulate_flow()
print("Watershed area above outlet [m^2]=",rmg.at_node['drainage_area'][outlet_id]) 
print("Log10 of Watershed area above outlet [m^2]=",np.log10(rmg.at_node['drainage_area'][outlet_id])) 

imshow_grid(rmg, drainage_area, plot_name = 'Upslope Catchment Area', 
            var_name = 'Upslope Catchment Area', var_units = 'm^2', grid_units = ('m','m'), 
            cmap = 'jet', limits = (0, 27000000))


imshow_grid(rmg, np.log10(drainage_area), plot_name = 'Upslope Catchment Area', 
            var_name = 'Upslope Catchment Area', var_units = 'm^2', grid_units = ('m','m'), 
            cmap = 'jet', limits = (0, 7.5))

critical_area = 10000000

rmg.status_at_node[drainage_area > critical_area] = rmg.BC_NODE_IS_FIXED_VALUE

(drainage_area, discharge) = fr.accumulate_flow()

# +
masks = get_watershed_masks(rmg)

rmg.add_field('node','mask',masks) 

imshow_grid(rmg, masks, cmap = 'jet')
# -

write_esri_ascii('./masks.asc', rmg, 'mask', clobber=True)
