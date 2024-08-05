# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 23:26:54 2024

@author: silve
"""

#%%
import pandas as pd
import glob
import os

#%%
#define base directory to my files
base_directory_flk2 = "C:/Users/silve/Desktop/FLK2_RSC"

#Define the patterns to match different sets of files
patterns_flk2 = [
    "M4/*.xlsx", "M5/*.xlsx", "M6/*.xlsx", "M7/*.xlsx", "M8/*.xlsx",
    "M9/*.xlsx"]

#create a dictionary that holds all of the DataFrames
dataframes_flk2 = {}

#Iterate trhough each pattern
for pattern_flk2 in patterns_flk2:
    #construct the search pattern
    search_pattern_flk2 = os.path.join(base_directory_flk2, pattern_flk2)
    #find excel files matching the pattern
    excel_files_flk2 = glob.glob(search_pattern_flk2)
    
    #Iterates through each file to read it into a Dataframe and store it in the dictionary
    for file_flk2 in excel_files_flk2:
        #extract folder name
        folder_name_flk2 = os.path.basename(os.path.dirname(file_flk2))
        #extract file name without extension
        file_name_flk2 = os.path.basename(file_flk2).replace('.xlsx', '')
        #determine the type of file depending on its name
        if "Results" in file_name_flk2:
            file_type = "Results"
        elif "AB" in file_name_flk2:
            file_type = "AB_coordinates"
        elif "M" in file_name_flk2:
            file_type = "M_coordinates"
        else:
            file_type = "Other"
        #Determine the hemisphere based on the file name pattern
        if file_name_flk2.endswith("1.1"):
            hemisphere = "Left"
        elif file_name_flk2.endswith("1.2"):
            hemisphere = "Right"
        else: 
            hemisphere = "Unknown"
        
        #creating a unique key for the dictionary to identify each file
        key = f"{folder_name_flk2}_{file_type}_{hemisphere}_FLK2_RSC"
        
        #read excel file into a DataFrame
        df_flk2 = pd.read_excel(file_flk2, engine = 'openpyxl')
        
        #store the DataFrame in the dictionary
        dataframes_flk2[key] = df_flk2
        
        #print confirmation
        print(f"Loaded {key} with shape {df_flk2.shape}")

print("Loaded DataFrames keys:", dataframes_flk2.keys())

for key in dataframes_flk2:
    #rename first column to 'ID'
    df_flk2 = dataframes_flk2[key]
    df_flk2.columns = ['ID'] + df_flk2.columns.tolist()[1:]
    
    #store the updated dataframe back in the dictionary
    dataframes_flk2[key] = df_flk2
#%%
#define base directory to my files
base_directory_flk2_ca1 = "C:/Users/silve/Desktop/FLK2_CA1"

#Define the patterns to match different sets of files
patterns_flk2_ca1 = [
    "M4/*.xlsx", "M5/*.xlsx", "M6/*.xlsx", "M7/*.xlsx", "M8/*.xlsx",
    "M9/*.xlsx"]

#create a dictionary that holds all of the DataFrames
dataframes_flk2_ca1 = {}

#Iterate trhough each pattern
for pattern_flk2_ca1 in patterns_flk2_ca1:
    #construct the search pattern
    search_pattern_flk2_ca1 = os.path.join(base_directory_flk2_ca1, pattern_flk2_ca1)
    #find excel files matching the pattern
    excel_files_flk2_ca1 = glob.glob(search_pattern_flk2_ca1)
    
    #Iterates through each file to read it into a Dataframe and store it in the dictionary
    for file_flk2_ca1 in excel_files_flk2_ca1:
        #extract folder name
        folder_name_flk2_ca1 = os.path.basename(os.path.dirname(file_flk2_ca1))
        #extract file name without extension
        file_name_flk2_ca1 = os.path.basename(file_flk2_ca1).replace('.xlsx', '')
        #determine the type of file depending on its name
        if "Results" in file_name_flk2_ca1:
            file_type = "Results"
        elif "AB" in file_name_flk2_ca1:
            file_type = "AB_coordinates"
        elif "M" in file_name_flk2_ca1:
            file_type = "M_coordinates"
        else:
            file_type = "Other"
        #Determine the hemisphere based on the file name pattern
        if file_name_flk2_ca1.endswith("1.1"):
            hemisphere = "Left"
        elif file_name_flk2_ca1.endswith("1.2"):
            hemisphere = "Right"
        else: 
            hemisphere = "Unknown"
        
        #creating a unique key for the dictionary to identify each file
        key = f"{folder_name_flk2_ca1}_{file_type}_{hemisphere}_FLK2_CA1"
        
        #read excel file into a DataFrame
        df_flk2_ca1 = pd.read_excel(file_flk2_ca1, engine = 'openpyxl')
        
        #store the DataFrame in the dictionary
        dataframes_flk2_ca1[key] = df_flk2_ca1
        
        #print confirmation
        print(f"Loaded {key} with shape {df_flk2_ca1.shape}")

print("Loaded DataFrames keys:", dataframes_flk2_ca1.keys())

for key in dataframes_flk2_ca1:
    #rename first column to 'ID'
    df_flk2_ca1 = dataframes_flk2_ca1[key]
    df_flk2_ca1.columns = ['ID'] + df_flk2_ca1.columns.tolist()[1:]
    
    #store the updated dataframe back in the dictionary
    dataframes_flk2_ca1[key] = df_flk2_ca1

#%% 
#denoising the results dataset by using the drop rows on conditions function. Gets rid of any endpoint lower than two and any cells whose longest shortes path is lower than the 90th quantile.
processed_dfs_flk2 = {}
for key, df_flk2 in dataframes_flk2.items():
    if "Results" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df_flk2.columns.tolist()}")
        # Ensure the column names exist in the DataFrame
        required_columns = ["ID", "# Branches", "# End-point voxels", "Average Branch Length", "Longest Shortest Path"]
        if all(col in df_flk2.columns for col in required_columns):
           new_df_flk2 = df_flk2[required_columns]
           #apply each condition to drop rows from the dataframes
           df_dropped_flk2 = drop_rows_on_conditions(new_df_flk2, conditions)
           #append the processed dataframes to the list
           processed_dfs_flk2[key] = df_dropped_flk2
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")

# Update the original dictionary with the processed DataFrames
dataframes_flk2.update(processed_dfs_flk2)
        
processed_dfs_ca1_flk2 = {}
for key, df_flk2_ca1 in dataframes_flk2_ca1.items():
    if "Results" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df_flk2_ca1.columns.tolist()}")
        # Ensure the column names exist in the DataFrame
        required_columns = ["ID", "# Branches", "# End-point voxels", "Average Branch Length", "Longest Shortest Path"]
        if all(col in df_flk2_ca1.columns for col in required_columns):
           new_df_ca1_flk2 = df_flk2_ca1[required_columns]
           #apply each condition to drop rows from the dataframes
           df_dropped_ca1_flk2 = drop_rows_on_conditions(new_df_ca1_flk2, conditions_ca1)
           #append the processed dataframes to the list
           processed_dfs_ca1_flk2[key] = df_dropped_ca1_flk2
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")
            
dataframes_flk2_ca1.update(processed_dfs_ca1_flk2)
#%%
#tidying the dataframes by getting rid of any unwanted columns
#starts a dictionary where to store the tidied dataframes
processed_m_flk2 = {}
#starts a loop which iterates through all of the dataframes in the dictionary
for key, df_flk2 in dataframes_flk2.items():
    #only selects the m and ab coordinates dataframes
    if "M_coordinates" or "AB_coordinates" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df_flk2.columns.tolist()}")
        #selects the required columns and ensures the column names exist in the DataFrame
        required_columns_m = ["ID", "Area", "X", "Y", "Feret"]
        if all(col in df_flk2.columns for col in required_columns_m):
           df_m_flk2 = df_flk2[required_columns_m]
           processed_m_flk2[key] = df_m_flk2
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")
#updates the original dictionary of dataframes with the tidied ones            
dataframes_flk2.update(processed_m_flk2)

for key, df_flk2 in dataframes_flk2.items():
    if "M_coordinates" in key:
        df_flk2['MicrogliaNumber'] = len(df_flk2)
        df_n_flk2 = df_flk2[required_columns_m + ['MicrogliaNumber']]
        dataframes_flk2[key] = df_n_flk2

processed_m_ca1_flk2 = {}
for key, df_flk2_ca1 in dataframes_flk2_ca1.items():
    if "M_coordinates" or "AB_coordinates" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df_flk2_ca1.columns.tolist()}")
        # Ensure the column names exist in the DataFrame
        required_columns_m = ["ID", "Area", "X", "Y", "Feret"]
        if all(col in df_flk2_ca1.columns for col in required_columns_m):
           df_m_ca1_flk2 = df_flk2_ca1[required_columns_m]
           processed_m_ca1_flk2[key] = df_m_ca1_flk2
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")
            
dataframes_flk2_ca1.update(processed_m_ca1_flk2)

for key, df_flk2_ca1 in dataframes_flk2_ca1.items():
    if "M_coordinates" in key:
        df_flk2_ca1['MicrogliaNumber'] = len(df_flk2_ca1)
        df_n_flk2_ca1 = df_flk2_ca1[required_columns_m + ['MicrogliaNumber']]
        dataframes_flk2_ca1[key] = df_n_flk2_ca1
#%%
dataframes_means_flk2 = {}
for key, df_flk2 in dataframes_flk2.items():
    df_mean_flk2 = df_flk2.iloc[:, 1:].mean()
    dataframes_means_flk2[key] = df_mean_flk2

dataframes_ca1_means_flk2 = {}
for key, df_flk2_ca1 in dataframes_flk2_ca1.items():
    df_mean_ca1_flk2 = df_flk2_ca1.iloc[:, 1:].mean()
    dataframes_ca1_means_flk2[key] = df_mean_ca1_flk2
#%%
dataframes_combined_flk2 = pd.concat(dataframes_means_flk2, axis = 0)

dataframes_ca1_combined_flk2 = pd.concat(dataframes_ca1_means_flk2, axis = 0)
#%%
# Create the MultiIndex Series
index_flk2 = pd.MultiIndex.from_tuples(dataframes_combined_flk2.keys(), names=['MouseID', 'Parameter'])
dfs_combined_wide_flk2 = pd.Series(dataframes_combined_flk2, index=index_flk2)

# Convert the Series to a DataFrame
df_wide_flk2 = dfs_combined_wide_flk2.unstack(level="MouseID")
print(df_wide_flk2)

# Create the MultiIndex Series
index_ca1_flk2 = pd.MultiIndex.from_tuples(dataframes_ca1_combined_flk2.keys(), names=['MouseID', 'Parameter'])
dfs_combined_wide_ca1_flk2 = pd.Series(dataframes_ca1_combined_flk2, index=index_ca1_flk2)

# Convert the Series to a DataFrame
df_wide_ca1_flk2 = dfs_combined_wide_ca1_flk2.unstack(level="MouseID")
print(df_wide_ca1_flk2)


#%%
#RSC
#identify the unique base names
base_names_flk2 = set(col.split("_")[0] for col in df_wide_flk2.columns)

#create a dictionary to store the new columns
new_columns_flk2 = {}

#loop trhough the base names to create new columns by averaging left and right values
for base_flk2 in base_names_flk2:
    left_cols_flk2 = [col for col in df_wide_flk2.columns if col.startswith(base_flk2 + "_") and "Left" in col]
    right_cols_flk2 = [col for col in df_wide_flk2.columns if col.startswith(base_flk2 + "_") and "Right" in col]
    for left_col_flk2, right_col_flk2 in zip(left_cols_flk2, right_cols_flk2):
        new_col_name_flk2 = left_col_flk2.replace("Left", "")
        new_columns_flk2[new_col_name_flk2] = (df_wide_flk2[left_col_flk2] + df_wide_flk2[right_col_flk2]) / 2

#Add the new columns to the df
for col_name_flk2, col_data_flk2 in new_columns_flk2.items():
    df_wide_flk2[col_name_flk2] = col_data_flk2

#drop the old left and right columns
df_wide_flk2 = df_wide_flk2.drop(columns = [col for col in df_wide_flk2.columns if "Left" in col or "Right" in col])

#CA1
#identify the unique base names
base_names_ca1_flk2 = set(col.split("_")[0] for col in df_wide_ca1_flk2.columns)

#create a dictionary to store the new columns
new_columns_ca1_flk2 = {}

#loop trhough the base names to create new columns by averaging left and right values
for base_ca1_flk2 in base_names_ca1_flk2:
    left_cols_ca1_flk2 = [col for col in df_wide_ca1_flk2.columns if col.startswith(base_ca1_flk2 + "_") and "Left" in col]
    right_cols_ca1_flk2 = [col for col in df_wide_ca1_flk2.columns if col.startswith(base_ca1_flk2 + "_") and "Right" in col]
    
    for left_col_ca1_flk2, right_col_ca1_flk2 in zip(left_cols_ca1_flk2, right_cols_ca1_flk2):
        new_col_name_ca1_flk2 = left_col_ca1_flk2.replace("Left", "")
        new_columns_ca1_flk2[new_col_name_ca1_flk2] = (df_wide_ca1_flk2[left_col_ca1_flk2] + df_wide_ca1_flk2[right_col_ca1_flk2]) / 2

#Add the new columns to the df
for col_name_ca1_flk2, col_data_ca1_flk2 in new_columns_ca1_flk2.items():
    df_wide_ca1_flk2[col_name_ca1_flk2] = col_data_ca1_flk2

#drop the old left and right columns
df_wide_ca1_flk2 = df_wide_ca1_flk2.drop(columns = [col for col in df_wide_ca1_flk2.columns if "Left" in col or "Right" in col])
#%%
# Transpose the DataFrame
df_transposed_flk2 = df_wide_flk2.T

# Reset the index to make the MouseID a column again
df_transposed_flk2.reset_index(inplace=True)

# Rename the columns if necessary (optional)
df_transposed_flk2.rename(columns={'index': 'MouseID'}, inplace=True)

print(df_transposed_flk2)

output_path_rsc_flk2 = r"C:\Users\silve\Desktop\FLK2_RSC\dataframes_flk2.xlsx"
df_transposed_flk2.to_excel(output_path_rsc_flk2, index=True)

# Transpose the DataFrame
df_transposed_ca1_flk2 = df_wide_ca1_flk2.T

# Reset the index to make the MouseID a column again
df_transposed_ca1_flk2.reset_index(inplace=True)

# Rename the columns if necessary (optional)
df_transposed_ca1_flk2.rename(columns={'index': 'MouseID'}, inplace=True)

print(df_transposed_ca1_flk2)

output_path_ca1_flk2 = r"C:\Users\silve\Desktop\FLK2_CA1\dataframes_ca1_flk2.xlsx"
df_transposed_ca1_flk2.to_excel(output_path_ca1_flk2, index=True)

#%%  Concatenates microglia and abeta dfs

#Combine all microglia and abeta points into two DataFrames
microglia_points_flk2 = []
abeta_points_flk2 = []

for key, df_flk2 in dataframes_flk2.items():
    if "M_coordinates" in key:
        df_flk2['Source'] = key  # Keep track of the source DataFrame
        microglia_points_flk2.append(df_flk2)
    elif "AB_coordinates" in key:
       df_flk2['Source'] = key  # Keep track of the source DataFrame
       abeta_points_flk2.append(df_flk2)

# Concatenate all microglia and abeta points into two DataFrames
microglia_df_flk2 = pd.concat(microglia_points_flk2, ignore_index=True)
abeta_df_flk2 = pd.concat(abeta_points_flk2, ignore_index=True)

#%% Calculates eucleadean distance
# Initialize lists to store results
microglia_ids_flk2 = []
abeta_ids_flk2 = []
distances_flk2 = []

# Vectorized calculation of Euclidean distances
microglia_coords_flk2 = microglia_df_flk2[['X', 'Y']].values
abeta_coords_flk2 = abeta_df_flk2[['X', 'Y']].values

# Calculate the Euclidean distance matrix
dist_matrix_flk2 = np.sqrt(((microglia_coords_flk2[:, np.newaxis, :] - abeta_coords_flk2[np.newaxis, :, :]) ** 2).sum(axis=2))

# Find the minimum distance and corresponding abeta ID for each microglia point
min_dist_indices_flk2 = np.argmin(dist_matrix_flk2, axis=1)
min_distances_flk2 = dist_matrix_flk2[np.arange(len(microglia_df_flk2)), min_dist_indices_flk2]

# Store the results
microglia_ids_flk2 = microglia_df_flk2['ID'].tolist()
abeta_ids_flk2 = abeta_df_flk2['ID'].iloc[min_dist_indices_flk2].tolist()
distances_flk2 = min_distances_flk2.tolist()

#%% Calculates number of paired abeta and microglia points 
# Create a DataFrame to store the results
df_distance_flk2 = pd.DataFrame({
    'Microglia_ID': microglia_ids_flk2,
    'Abeta_ID': abeta_ids_flk2,
    'Distance': distances_flk2
})

# Calculate the number of unique abeta points that have been paired
unique_abeta_count_flk2 = df_distance_flk2['Abeta_ID'].nunique()


# Add this count as a new column to the results_df
df_distance_flk2['Unique_Abeta_Count'] = unique_abeta_count_flk2

# Filter the results to retain only rows with distance less than 100
filtered_distances_flk2 = df_distance_flk2[df_distance_flk2['Distance'] < 100]

# Summation of distances and count of unique Microglia_ID entries in the filtered dataset
distance_sum_flk2 = filtered_distances_flk2['Distance'].sum()
unique_microglia_id_count_flk2 = filtered_distances_flk2['Microglia_ID'].nunique()

print(f"Sum of distances (excluding those over 100): {distance_sum_flk2}")
print(f"Number of unique Microglia_ID entries: {unique_microglia_id_count_flk2}")

# Save the updated dataframe to a excel workbook file
output_path_updated_flk2 = r"C:\Users\silve\Desktop\FLK2_RSC\microglia_abeta_alignment.xlsx"
df_distances_flk2.to_excel(output_path_updated_flk2, index=False)

output_path_flk2 = r"C:\Users\silve\Desktop\FLK2_RSC\microglia_abeta_filtered.xlsx"
filtered_distances_flk2.to_excel(output_path_flk2, index=False)

# Group by Abeta_ID and count the number of Microglia aligned to each Abeta
abeta_microglia_count_flk2 = filtered_distances_flk2.groupby('Abeta_ID').size().reset_index(name='Microglia_Count')

# Save the resulting dataframe to excel workbook file
output_path_abeta_count_flk2 = r"C:\Users\silve\Desktop\FLK2_RSC\abeta_microglia_count.xlsx"
abeta_microglia_count_flk2.to_excel(output_path_abeta_count_flk2, index=False)

# Calculate and print the summation of Microglia_Count and number of Abeta_ID in abeta_microglia_count
total_microglia_aligned_flk2 = abeta_microglia_count_flk2['Microglia_Count'].sum()
number_of_abeta_id_flk2 = len(abeta_microglia_count_flk2)

print(f"Total number of microglia aligned to all abeta points: {total_microglia_aligned_flk2}")
print(f"Number of unique Abeta_ID entries: {number_of_abeta_id_flk2}")
#%%  Concatenates microglia and abeta dfs

#Combine all microglia and abeta points into two DataFrames
microglia_points_flk2_ca1 = []
abeta_points_flk2_ca1 = []

for key, df_flk2_ca1 in dataframes_flk2_ca1.items():
    if "M_coordinates" in key:
        df_flk2_ca1['Source'] = key  # Keep track of the source DataFrame
        microglia_points_flk2_ca1.append(df_flk2_ca1)
    elif "AB_coordinates" in key:
       df_flk2_ca1['Source'] = key  # Keep track of the source DataFrame
       abeta_points_flk2_ca1.append(df_flk2_ca1)

# Concatenate all microglia and abeta points into two DataFrames
microglia_df_flk2_ca1 = pd.concat(microglia_points_flk2_ca1, ignore_index=True)
abeta_df_flk2_ca1 = pd.concat(abeta_points_flk2_ca1, ignore_index=True)

#%% Calculates eucleadean distance
# Initialize lists to store results
microglia_ids_flk2_ca1 = []
abeta_ids_flk2_ca1 = []
distances_flk2_ca1 = []

# Vectorized calculation of Euclidean distances
microglia_coords_flk2_ca1 = microglia_df_flk2_ca1[['X', 'Y']].values
abeta_coords_flk2_ca1 = abeta_df_flk2_ca1[['X', 'Y']].values

# Calculate the Euclidean distance matrix
dist_matrix_flk2_ca1 = np.sqrt(((microglia_coords_flk2_ca1[:, np.newaxis, :] - abeta_coords_flk2_ca1[np.newaxis, :, :]) ** 2).sum(axis=2))

# Find the minimum distance and corresponding abeta ID for each microglia point
min_dist_indices_flk2_ca1 = np.argmin(dist_matrix_flk2_ca1, axis=1)
min_distances_flk2_ca1 = dist_matrix_flk2_ca1[np.arange(len(microglia_df_flk2_ca1)), min_dist_indices_flk2_ca1]

# Store the results
microglia_ids_flk2_ca1 = microglia_df_flk2_ca1['ID'].tolist()
abeta_ids_flk2_ca1 = abeta_df_flk2_ca1['ID'].iloc[min_dist_indices_flk2_ca1].tolist()
distances_flk2_ca1 = min_distances_flk2_ca1.tolist()
#%%
# Create a DataFrame to store the results
df_distance_flk2_ca1 = pd.DataFrame({
    'Microglia_ID': microglia_ids_flk2_ca1,
    'Abeta_ID': abeta_ids_flk2_ca1,
    'Distance': distances_flk2_ca1
})

# Calculate the number of unique abeta points that have been paired
unique_abeta_count_flk2_ca1 = df_distance_flk2_ca1['Abeta_ID'].nunique()


# Add this count as a new column to the results_df
df_distance_flk2_ca1['Unique_Abeta_Count'] = unique_abeta_count_flk2_ca1

# Filter the results to retain only rows with distance less than 100
filtered_distances_flk2_ca1 = df_distance_flk2_ca1[df_distance_flk2_ca1['Distance'] < 100]

# Summation of distances and count of unique Microglia_ID entries in the filtered dataset
distance_sum_flk2_ca1 = filtered_distances_flk2_ca1['Distance'].sum()
unique_microglia_id_count_flk2_ca1 = filtered_distances_flk2_ca1['Microglia_ID'].nunique()

print(f"Sum of distances (excluding those over 100): {distance_sum_flk2_ca1}")
print(f"Number of unique Microglia_ID entries: {unique_microglia_id_count_flk2_ca1}")

# Save the updated dataframe to a excel workbook file
output_path_updated_flk2_ca1 = r"C:\Users\silve\Desktop\FLK2_CA1\microglia_abeta_alignment.xlsx"
df_distance_flk2_ca1.to_excel(output_path_updated_flk2_ca1, index=False)

output_path_flk2_ca1 = r"C:\Users\silve\Desktop\FLK2_CA1\microglia_abeta_filtered.xlsx"
filtered_distances_flk2_ca1.to_excel(output_path_flk2_ca1, index=False)

# Group by Abeta_ID and count the number of Microglia aligned to each Abeta
abeta_microglia_count_flk2_ca1 = filtered_distances_flk2_ca1.groupby('Abeta_ID').size().reset_index(name='Microglia_Count')

# Save the resulting dataframe to excel workbook file
output_path_abeta_count_flk2_ca1 = r"C:\Users\silve\Desktop\FLK2_CA1\abeta_microglia_count.xlsx"
abeta_microglia_count_flk2_ca1.to_excel(output_path_abeta_count_flk2_ca1, index=False)

# Calculate and print the summation of Microglia_Count and number of Abeta_ID in abeta_microglia_count
total_microglia_aligned_flk2_ca1 = abeta_microglia_count_flk2_ca1['Microglia_Count'].sum()
number_of_abeta_id_flk2_ca1 = len(abeta_microglia_count_flk2_ca1)

print(f"Total number of microglia aligned to all abeta points: {total_microglia_aligned_flk2_ca1}")
print(f"Number of unique Abeta_ID entries: {number_of_abeta_id_flk2_ca1}")