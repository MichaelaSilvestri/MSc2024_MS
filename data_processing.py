# -*- coding: utf-8 -*-
"""
Created on Wed Jul 10 10:24:52 2024

@author: silve
"""
#%%
import pandas as pd
import glob
import os
import numpy as np

#importing raw data into a pandas data frame
mbody = pd.read_excel('RSC\MG0\Results_MG0_1.1.xlsx', engine = 'openpyxl')
print(mbody)

    
#%%
#denoising the dataset by removing everything with less than 2 endpoints
new_mbody = mbody[mbody["# End-point voxels"] > 2]

#%%
#Calculating the 90th quantile
#first ordering the column in ascending order
longshort_sorted = new_mbody.sort_values(by = 'Longest Shortest Path')

#then calculating the position of the 90th quantile
percentile_90 = longshort_sorted['Longest Shortest Path'].quantile(0.90)

#using the 90th quantile as a cut-off point to delete noise to make sure it works
#filtered_mbody = longshort_sorted[longshort_sorted["Longest Shortest Path"] > percentile_90]   

#%%
mbodyr = pd.read_excel('RSC\MG0\Results_MG0_1.2.xlsx', engine = 'openpyxl')

new_mbodyr = mbodyr[mbodyr["# End-point voxels"] > 2]

#%%
#Calculating the 90th quantile
#first ordering the column in ascending order
longshort_sortedr = new_mbodyr.sort_values(by = 'Longest Shortest Path')

#then calculating the position of the 90th quantile
percentile_90r = longshort_sortedr['Longest Shortest Path'].quantile(0.90)

#using the 90th quantile as a cut-off point to delete noise
#filtered_mbodyr = longshort_sortedr[longshort_sortedr["Longest Shortest Path"] > percentile_90]

#%%
#create an object called cut-off to use inside the condition function
cutoff = (percentile_90 + percentile_90r) / 2

#%%
#now try write the above code in a function/for loop so that I can use it for all the datasets
#creating a function that drops rows based on a specific condition 
def drop_rows_on_conditions(df, conditions):
    """
    Drops rows in a dataframe if they meet a specified condition.

    Parameters
    ----------
    df (pd.Dataframe): The input DataFrame. 
    column_name (str): The column name which values to check the conditions against.
    conditions (list of tuple): A list of conditions where each condition is a tuple (column_name, condition_function)

    Returns
    pd.DataFrame: The DataFrame with rows droppend based on the conditions.
    """
    for column_name, condition in conditions:
        if column_name in df.columns:
            #create an object that stores the indexes of the rows to drop by applying the condition.
            rows_to_drop = df[df[column_name].apply(condition)].index   
            #uses the indexes stored in the object above to drop the rows. 
            df = df.drop(rows_to_drop)
    return df

def condition1(value):
    return value <= 2 #returns True if the value is less than 2.
def condition2(value):
    return value < cutoff #Returns True is the value is less than the cutoff.
def condition3(value):
    return value < cutoff_ca1


#create a list of conditions 
conditions = [
    ("# End-point voxels", condition1),
    ("Longest Shortest Path", condition2)]

conditions_ca1 = [
    ("# End-point voxels", condition1),
    ("Longest Shortest Path", condition3)]


#%%
#We can do the above process quicker by creating a dictionary establishing a common pattern for the directory of our files. 
#We can make the directory a bit more dynamic so that that the for loop can iterate and look for patterns in the file directory without needing to specify the fool file directory

#define base directory to my files
base_directory = "C:/Users/silve/Desktop/RSC"

#Define the patterns to match different sets of files
patterns = [
    "M03/*.xlsx", "M05/*.xlsx", "M06/*.xlsx", "M07/*.xlsx", "M09/*.xlsx",
    "M11/*.xlsx", "M12/*.xlsx", "M13/*.xlsx", "MA0/*.xlsx", "MA0M/*.xlsx",
    "MA1/*.xlsx", "MB0/*.xlsx", "MB1/*.xlsx", "MB2/*.xlsx", "MC1/*.xlsx",
    "MC2/*.xlsx", "MD0/*.xlsx", "MD0M/*.xlsx", "MD1/*.xlsx", "ME3/*.xlsx",
    "ME10/*.xlsx", "MF1/*.xlsx", "MF2/*.xlsx", "MF3/*.xlsx", "MG0/*.xlsx", 
    "MG1/*.xlsx"]

#create a dictionary that holds all of the DataFrames
dataframes = {}

#Iterate trhough each pattern
for pattern in patterns:
    #construct the search pattern
    search_pattern = os.path.join(base_directory, pattern)
    #find excel files matching the pattern
    excel_files = glob.glob(search_pattern)
    
    #Iterates through each file to read it into a Dataframe and store it in the dictionary
    for file in excel_files:
        #extract folder name
        folder_name = os.path.basename(os.path.dirname(file))
        #extract file name without extension
        file_name = os.path.basename(file).replace('.xlsx', '') #replaces first argument (old) with second (new)
        #determine the type of file depending on its name
        if "Results" in file_name:
            file_type = "Results"
        elif "Branch" in file_name:
            file_type = "Branchinfo"
        elif "AB" in file_name:
            file_type = "AB_coordinates"
        elif "M" in file_name:
            file_type = "M_coordinates"
        else:
            file_type = "Other"
        #Determine the hemisphere based on the file name pattern
        if file_name.endswith("1.1"):
            hemisphere = "Left"
        elif file_name.endswith("1.2"):
            hemisphere = "Right"
        else: 
            hemisphere = "Unknown"
        
        #creating a unique key for the dictionary to identify each file
        key = f"{folder_name}_{file_type}_{hemisphere}"
        
        #read excel file into a DataFrame
        df = pd.read_excel(file, engine = 'openpyxl')
        
        #store the DataFrame in the dictionary
        dataframes[key] = df
        
        #print confirmation
        print(f"Loaded {key} with shape {df.shape}")

print("Loaded DataFrames keys:", dataframes.keys())

for key in dataframes:
    #rename first column to 'ID'
    df = dataframes[key]
    df.columns = ['ID'] + df.columns.tolist()[1:]
    
    #store the updated dataframe back in the dictionary
    dataframes[key] = df
        
 #%%
results_dfs = {key: df for key, df in dataframes.items() if "Results" in key}
processed_dataframes = {}
for key, df in results_dfs.items():
    #could use a print statement to verify if the code is working 
    print(f"Processing Dataframe {key} with columns: {df.columns.tolist()}")
    #apply each condition to drop rows from the dataframes
    df = drop_rows_on_conditions(df, conditions)
    #append the processed dataframes to the list
    processed_dataframes[key] = df
    

#%%  Concatenates microglia and abeta dfs

#Combine all microglia and abeta points into two DataFrames
microglia_points = []
abeta_points = []

for key, df in dataframes.items():
    if "M_coordinates" in key:
        df['Source'] = key  # Keep track of the source DataFrame
        microglia_points.append(df)
    elif "AB_coordinates" in key:
       df['Source'] = key  # Keep track of the source DataFrame
       abeta_points.append(df)

# Concatenate all microglia and abeta points into two DataFrames
microglia_df = pd.concat(microglia_points, ignore_index=True)
abeta_df = pd.concat(abeta_points, ignore_index=True)

#%% Calculates eucleadean distance
# Initialize lists to store results
microglia_ids = []
abeta_ids = []
distances = []

# Vectorized calculation of Euclidean distances
microglia_coords = microglia_df[['X', 'Y']].values
abeta_coords = abeta_df[['X', 'Y']].values

# Calculate the Euclidean distance matrix
dist_matrix = np.sqrt(((microglia_coords[:, np.newaxis, :] - abeta_coords[np.newaxis, :, :]) ** 2).sum(axis=2))

# Find the minimum distance and corresponding abeta ID for each microglia point
min_dist_indices = np.argmin(dist_matrix, axis=1)
min_distances = dist_matrix[np.arange(len(microglia_df)), min_dist_indices]

# Store the results
microglia_ids = microglia_df['ID'].tolist()
abeta_ids = abeta_df['ID'].iloc[min_dist_indices].tolist()
distances = min_distances.tolist()

#%% Calculates number of paired abeta and microglia points 
# Create a DataFrame to store the results
df_distance = pd.DataFrame({
    'Microglia_ID': microglia_ids,
    'Abeta_ID': abeta_ids,
    'Distance': distances
})

# Calculate the number of unique abeta points that have been paired
unique_abeta_count = df_distance['Abeta_ID'].nunique()

#need to fix this part of the code. 
#for each abeta id we want to know how many direct connections to microglia points there are. 
#not how many direct connections in total.

# Add this count as a new column to the results_df
df_distance['Unique_Abeta_Count'] = unique_abeta_count

# Filter the results to retain only rows with distance less than 100
filtered_distances = df_distance[df_distance['Distance'] < 100]

# Summation of distances and count of unique Microglia_ID entries in the filtered dataset
distance_sum = filtered_distances['Distance'].sum()
unique_microglia_id_count = filtered_distances['Microglia_ID'].nunique()

print(f"Sum of distances (excluding those over 100): {distance_sum}")
print(f"Number of unique Microglia_ID entries: {unique_microglia_id_count}")

# Save the updated dataframe to a excel workbook file
output_path_updated = r"C:\Users\silve\Desktop\RSC\microglia_abeta_alignment.xlsx"
df_distance.to_excel(output_path_updated, index=False)

# Group by Abeta_ID and count the number of Microglia aligned to each Abeta
abeta_microglia_count = filtered_distances.groupby('Abeta_ID').size().reset_index(name='Microglia_Count')

# Save the resulting dataframe to excel workbook file
output_path_abeta_count = r"C:\Users\silve\Desktop\RSC\abeta_microglia_count.xlsx"
abeta_microglia_count.to_excel(output_path_abeta_count, index=False)

# Calculate and print the summation of Microglia_Count and number of Abeta_ID in abeta_microglia_count
total_microglia_aligned = abeta_microglia_count['Microglia_Count'].sum()
number_of_abeta_id = len(abeta_microglia_count)

print(f"Total number of microglia aligned to all abeta points: {total_microglia_aligned}")
print(f"Number of unique Abeta_ID entries: {number_of_abeta_id}")
#%%
#importing raw data into a pandas data frame
ca1mbody = pd.read_excel('CA1\MG0\Results_MG0_1.1.xlsx', engine = 'openpyxl')
print(ca1mbody)

    
#%%
#denoising the dataset by removing everything with less than 2 endpoints
ca1mbody_new = ca1mbody[ca1mbody["# End-point voxels"] > 2]

#%%
#Calculating the 90th quantile
#first ordering the column in ascending order
longshort_sorted_ca1 = ca1mbody_new.sort_values(by = 'Longest Shortest Path')

#then calculating the position of the 90th quantile
percentile_90_ca1 = longshort_sorted_ca1['Longest Shortest Path'].quantile(0.90)

#using the 90th quantile as a cut-off point to delete noise
#filtered_mbody_ca1 = longshort_sorted_ca1[longshort_sorted_ca1["Longest Shortest Path"] > percentile_90_ca1]   

#%%
ca1mbodyr = pd.read_excel('CA1\MG0\Results_MG0_1.2.xlsx', engine = 'openpyxl')

ca1mbodyr_new = ca1mbodyr[ca1mbodyr["# End-point voxels"] > 2]

#%%
#Calculating the 90th quantile
#first ordering the column in ascending order
longshort_sortedr_ca1 = ca1mbodyr_new.sort_values(by = 'Longest Shortest Path')

#then calculating the position of the 90th quantile
percentile_90r_ca1 = longshort_sortedr_ca1['Longest Shortest Path'].quantile(0.90)

#using the 90th quantile as a cut-off point to delete noise
#filtered_mbodyr_ca1 = longshort_sortedr_ca1[longshort_sortedr_ca1["Longest Shortest Path"] > percentile_90_ca1]

#%%
#create an object called cut-off to use inside the condition function
cutoff_ca1 = (percentile_90_ca1 + percentile_90r_ca1) / 2

#%%
#define base directory to my files
base_directory_ca1 = "C:/Users/silve/Desktop/CA1"

#Define the patterns to match different sets of files
patterns_ca1 = [
    "M03/*.xlsx", "M05/*.xlsx", "M06/*.xlsx", "M07/*.xlsx", "M09/*.xlsx",
    "M11/*.xlsx", "M12/*.xlsx", "M13/*.xlsx", "MA0/*.xlsx", "MA0M/*.xlsx",
    "MA1/*.xlsx", "MB0/*.xlsx", "MB1/*.xlsx", "MB2/*.xlsx", "MC1/*.xlsx",
    "MC2/*.xlsx", "MD0/*.xlsx", "MD0M/*.xlsx", "MD1/*.xlsx", "ME3/*.xlsx",
    "ME10/*.xlsx", "MF1/*.xlsx", "MF2/*.xlsx", "MF3/*.xlsx", "MG0/*.xlsx", 
    "MG1/*.xlsx"]

#create a dictionary that holds all of the DataFrames
dataframes_ca1 = {}

#Iterate trhough each pattern
for pattern_ca1 in patterns_ca1:
    #construct the search pattern
    search_pattern_ca1 = os.path.join(base_directory_ca1, pattern_ca1)
    #find excel files matching the pattern
    excel_files_ca1 = glob.glob(search_pattern_ca1)
    
    #Iterates through each file to read it into a Dataframe and store it in the dictionary
    for file_ca1 in excel_files_ca1:
        #extract folder name
        folder_name_ca1 = os.path.basename(os.path.dirname(file_ca1))
        #extract file name without extension
        file_name_ca1 = os.path.basename(file_ca1).replace('.xlsx', '') #replaces first argument (old) with second (new)
        #determine the type of file depending on its name
        if "Results" in file_name_ca1:
            file_type = "Results"
        elif "Branch" in file_name_ca1:
            file_type = "Branchinfo"
        elif "AB" in file_name_ca1:
            file_type = "AB_coordinates" 
        elif "M" in file_name_ca1:
            file_type = "M_coordinates"
        else:
            file_type = "Other"
        #Determine the hemisphere based on the file name pattern
        if file_name_ca1.endswith("1.1"):
            hemisphere = "Left"
        elif file_name_ca1.endswith("1.2"):
            hemisphere = "Right"
        else: 
            hemisphere = "Unknown"
        
        #creating a unique key for the dictionary to identify each file
        key = f"{folder_name_ca1}_{file_type}_{hemisphere}"
        
        #read excel file into a DataFrame
        df_ca1 = pd.read_excel(file_ca1, engine = 'openpyxl')
        
        #store the DataFrame in the dictionary
        dataframes_ca1[key] = df_ca1
        
        #print confirmation
        print(f"Loaded {key} with shape {df.shape}")

print("Loaded DataFrames keys:", dataframes_ca1.keys())

for key in dataframes_ca1:
    #rename first column to 'ID'
    df_ca1 = dataframes_ca1[key]
    df_ca1.columns = ['ID'] + df_ca1.columns.tolist()[1:]
    
    #store the updated dataframe back in the dictionary
    dataframes_ca1[key] = df_ca1
        
 #%%

#create an object with only the results dataframes 
results_ca1 = {key: df_ca1 for key, df_ca1 in dataframes_ca1.items() if "Results" in key}

#starts a dictionary to store all the processed results dataframes
processed_dataframes_ca1 = {}

#then create a loop which drops the rows in all the dataframes at once based on the specified conditions
for key, df_ca1 in results_ca1.items():
    #could use a print statement to verify if the code is working 
    print(f"Processing Dataframe {key} with columns: {df_ca1.columns.tolist()}")
    #apply each condition to drop rows from the dataframes
    df_ca1 = drop_rows_on_conditions(df_ca1, conditions_ca1)
    #append the processed dataframes to the list
    processed_dataframes_ca1[key] = df_ca1
    
#%%  Concatenates microglia and abeta dfs

#Combine all microglia and abeta points into two DataFrames
microglia_points_ca1 = []
abeta_points_ca1 = []

for key, df_ca1 in dataframes_ca1.items():
    if "M_coordinates" in key:
        df_ca1['Source'] = key  # Keep track of the source DataFrame
        microglia_points_ca1.append(df_ca1)
    elif "AB_coordinates" in key:
       df_ca1['Source'] = key  # Keep track of the source DataFrame
       abeta_points_ca1.append(df_ca1)

# Concatenate all microglia and abeta points into two DataFrames
microglia_df_ca1 = pd.concat(microglia_points_ca1, ignore_index=True)
abeta_df_ca1 = pd.concat(abeta_points_ca1, ignore_index=True)

#%% Calculates eucleadean distance
# Initialize lists to store results
microglia_ids_ca1 = []
abeta_ids_ca1 = []
distances_ca1 = []

# Vectorized calculation of Euclidean distances
microglia_coords_ca1 = microglia_df_ca1[['X', 'Y']].values
abeta_coords_ca1 = abeta_df_ca1[['X', 'Y']].values

# Calculate the Euclidean distance matrix
dist_matrix_ca1 = np.sqrt(((microglia_coords_ca1[:, np.newaxis, :] - abeta_coords_ca1[np.newaxis, :, :]) ** 2).sum(axis=2))

# Find the minimum distance and corresponding abeta ID for each microglia point
min_dist_indices_ca1 = np.argmin(dist_matrix_ca1, axis=1)
min_distances_ca1 = dist_matrix_ca1[np.arange(len(microglia_df_ca1)), min_dist_indices_ca1]

# Store the results
microglia_ids_ca1 = microglia_df_ca1['ID'].tolist()
abeta_ids_ca1 = abeta_df_ca1['ID'].iloc[min_dist_indices_ca1].tolist()
distances_ca1 = min_distances_ca1.tolist()
  
#%% Calculates number of paired abeta and microglia points 
# Create a DataFrame to store the results
df_distance_ca1 = pd.DataFrame({
    'Microglia_ID': microglia_ids_ca1,
    'Abeta_ID': abeta_ids_ca1,
    'Distance': distances_ca1
})

# Calculate the number of unique abeta points that have been paired
unique_abeta_count_ca1 = df_distance_ca1['Abeta_ID'].nunique()

#NEED TO FIX THIS CODE. FOR ECH ABETA POINTS WE WANT TO KNOW THE NUMBER OF MICROGLIA DIRECTLY CONNECTED TO IT. 
#NOT HOW MANY CONNECTIONS THERE ARE IN TOTAL.

# Add this count as a new column to the results_df
df_distance_ca1['Unique_Abeta_Count'] = unique_abeta_count_ca1

# Filter the results to retain only rows with distance less than 100
filtered_distances_ca1 = df_distance_ca1[df_distance_ca1['Distance'] < 100]

# Summation of distances and count of unique Microglia_ID entries in the filtered dataset
distance_sum_ca1 = filtered_distances_ca1['Distance'].sum()
unique_microglia_id_count_ca1 = filtered_distances_ca1['Microglia_ID'].nunique()

print(f"Sum of distances (excluding those over 100): {distance_sum_ca1}")
print(f"Number of unique Microglia_ID entries: {unique_microglia_id_count_ca1}")

# Save the updated dataframe to a excel workbook file
output_path_ca1 = r"C:\Users\silve\Desktop\CA1\microglia_abeta_alignment.xlsx"
df_distance_ca1.to_excel(output_path_ca1, index=False)

# Group by Abeta_ID and count the number of Microglia aligned to each Abeta
abeta_microglia_count_ca1 = filtered_distances_ca1.groupby('Abeta_ID').size().reset_index(name='Microglia_Count')

# Save the resulting dataframe to excel workbook file
output_path_abeta_count_ca1 = r"C:\Users\silve\Desktop\CA1\abeta_microglia_count.xlsx"
abeta_microglia_count_ca1.to_excel(output_path_abeta_count_ca1, index=False)

# Calculate and print the summation of Microglia_Count and number of Abeta_ID in abeta_microglia_count
total_microglia_aligned_ca1 = abeta_microglia_count_ca1['Microglia_Count'].sum()
number_of_abeta_id_ca1 = len(abeta_microglia_count_ca1)

print(f"Total number of microglia aligned to all abeta points: {total_microglia_aligned_ca1}")
print(f"Number of unique Abeta_ID entries: {number_of_abeta_id_ca1}")  


