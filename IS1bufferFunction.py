### 
import pyarrow.parquet as pq
from shapely.geometry import Polygon
from shapely.geometry import Point
import  geopandas as gpd
def is1_buffer(center_lon,  center_lat):
       file = '/gpfs/data1/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Data/is1.parquet'
       # Read the Parquet file
       table = pq.read_table(file)
       
       # Convert the table to a Pandas DataFrame
       df = table.to_pandas()
       geometry = [Point(xy) for xy in zip(df['lon'], df['lat'])]
       
       # Create a GeoDataFrame from the DataFrame
       gdf = gpd.GeoDataFrame(df, geometry=geometry)
       
       
       # Now you can work with the DataFrame as needed
       #print(df.head())
       
       # Specify the center latitude and longitude
       #center_lat = -6.16
       #center_lon = -54.54
       
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
       
       if clipped_gdf.empty:
                 return None 
       
       # Print the clipped GeoDataFrame
       #print(clipped_gdf)
       
       ### buffer ????
       # 0.0044915 --- 1km   diameter
       
       grid_buffer = clipped_gdf.buffer(0.0044915/2)
       
       ### union 
       #### geoseries 
       # Convert GeoSeries to GeoDataFrame
       gdf_buffer = gpd.GeoDataFrame(geometry=grid_buffer)
       
       dissolved = gdf_buffer.dissolve()
               # make tile name using bottom left corner coordinate, removing negative sign and adding cardinal direction
       if center_lon >= 0:
            lon_label = 'E_'
       else:
            lon_label = 'W_'
       if center_lat >= 0:
            lat_label = 'N.geojson'
       else:
            lat_label = 'S.geojson'

       tile = str(abs(center_lon)) + lon_label + str(abs(center_lat)) + lat_label

       # make file name (and set a path if desired) based on GEDI data product and tile
       file_out = '/gpfs/data1/vclgp/xiongl/HeightComparisonGEDI_IS2_IS1/Out/IS1buffer_'+ tile
       dissolved.to_file(file_out, driver='GeoJSON', mode='w')

import argparse

def process_arguments():
    parser = argparse.ArgumentParser(description='Process arguments')
    parser.add_argument('--lon', type=float, help='Your lon')
    parser.add_argument('--lat', type=float, help='Your lat')

    args = parser.parse_args()

    return args

def main():
    # Get arguments from command line
    args = process_arguments()

    # Access the argument values
    lon = args.lon
    lat = args.lat

    # Perform some operation using the arguments
    print("Lon:", lon)
    print("Lat:", lat)
    is1_buffer(lon,  lat)    

if __name__ == '__main__':
    main()




