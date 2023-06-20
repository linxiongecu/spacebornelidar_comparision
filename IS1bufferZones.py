### 
import pyarrow.parquet as pq

from shapely.geometry import Polygon
from shapely.geometry import Point
import  geopandas as gpd





file = '/gpfs/data1/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/is1.parquet'
# Read the Parquet file
table = pq.read_table(file)

# Convert the table to a Pandas DataFrame
df = table.to_pandas()



geometry = [Point(xy) for xy in zip(df['lon'], df['lat'])]

# Create a GeoDataFrame from the DataFrame
gdf = gpd.GeoDataFrame(df, geometry=geometry)


# Now you can work with the DataFrame as needed
print(df.head())


# Specify the center latitude and longitude
center_lat = -6.16
center_lon = -54.54

# Calculate the coordinates of the square vertices
side_length = 0.125

half_side_length = side_length / 2
min_lon = center_lon - half_side_length
max_lon = center_lon + half_side_length
min_lat = center_lat - half_side_length
max_lat = center_lat + half_side_length


# Create a GeoDataFrame with the polygon

# Print the first few rows of the GeoDataFrame
# Save the GeoDataFrame to a shapefile
#output_shapefile = 'square_polygon.shp'
#gdf.to_file(output_shapefile)


# Clip the GeoDataFrame
clipped_gdf = gpd.clip(gdf, (min_lon, min_lat, max_lon, max_lat))

# Print the clipped GeoDataFrame
print(clipped_gdf)

### buffer ????
# 0.0044915 --- 1km   diameter

grid_buffer = clipped_gdf.buffer(0.0044915/2)

### union 
#### geoseries 
# Convert GeoSeries to GeoDataFrame
gdf_buffer = gpd.GeoDataFrame(geometry=grid_buffer)

dissolved = gdf_buffer.dissolve()


dissolved.to_file('bufferfile.geojson', driver='GeoJSON', mode='w')

