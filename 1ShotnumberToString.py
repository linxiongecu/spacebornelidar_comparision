import pandas as pd
import pyarrow.parquet as pq

# Specify the Parquet file path
file_path = "calval_hansen_forest_loss_20230425.parquet"

# Read the Parquet file
parquet_data = pd.read_parquet(file_path)

# Convert the Parquet data to a string
data_string = parquet_data.to_string()

# Create a DataFrame from the data string
string_df = pd.DataFrame({'data': [data_string]})

# Specify the output Parquet file path
output_file_path = "calval_hansen_forest_loss_20230425_string.parquet"

# Write the DataFrame to a Parquet file
table = pq.from_pandas(string_df)
pq.write_table(table, output_file_path)

