import geopandas as gpd
import rasterio
import rasterio.mask
import numpy as np


shapefile_path = r"C:\Users\User\Desktop\Gaia - GIS DATA\Basemaps\CARs_PA\Shapefile\CARs_PA_SIRGAS2000_reprojetado.shp"
raster_path = r"C:\Users\User\Desktop\Gaia - GIS DATA\LULC - Rasters\MapBiomas\Raster\MapBiomas_PA_reprojetado.tif"

# Load Geospatial Data
gdf = gpd.read_file(shapefile_path)
raster = rasterio.open(raster_path)

# Double check for same projection for both layers
if gdf.crs != raster.crs:
    gdf = gdf.to_crs(raster.crs)

# Calculate area by id
def calculate_area(geometry, raster, ids):
    # Crop Raster by shp mask
    out_image, out_transform = rasterio.mask.mask(raster, [geometry], crop=True)
    out_image = out_image[0] # as we have just one band we have to get it from the array
    
    pixel_size = out_transform[0] * out_transform[4]
    area = np.sum(np.isin(out_image, ids)) * abs(pixel_size)
    return area

# Set area data on shp .dbf
gdf['mb_re_elegible'] = 0.0
gdf['mb_re_not_elegible'] = 0.0

# Iterate over all geometries
total_polygons = len(gdf)
for idx, row in gdf.iterrows():
    geometry = row.geometry
    eligible_area = calculate_area(geometry, raster, [4, 5]) / 10000
    not_eligible_area = calculate_area(geometry, raster, [8]) / 10000

    gdf.at[idx, 'mb_re_elegible'] = eligible_area
    gdf.at[idx, 'mb_re_not_elegible'] = not_eligible_area

    # Print on terminal progress
    if idx % 100 == 0: # Print every 100 geometries
        print(f"Processed {idx+1} of {total_polygons} poligons ({(idx+1)/total_polygons*100:.2f}%)")

# Salvar o resultado
gdf.to_file(shapefile_path, driver='ESRI Shapefile')

print("Finished")