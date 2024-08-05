# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 16:13:54 2024

@author: silve
"""
#%%
import pandas as pd
import numpy as np
import glob 
import os
#%%
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
        key = f"{folder_name}_{file_type}_{hemisphere}_RSC"
        
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
        key = f"{folder_name_ca1}_{file_type}_{hemisphere}_CA1"
        
        #read excel file into a DataFrame
        df_ca1 = pd.read_excel(file_ca1, engine = 'openpyxl')
        
        #store the DataFrame in the dictionary
        dataframes_ca1[key] = df_ca1
        
        #print confirmation
        print(f"Loaded {key} with shape {df_ca1.shape}")

print("Loaded DataFrames keys:", dataframes_ca1.keys())

for key in dataframes_ca1:
    #rename first column to 'ID'
    df_ca1 = dataframes_ca1[key]
    df_ca1.columns = ['ID'] + df_ca1.columns.tolist()[1:]
    #store the updated dataframe back in the dictionary
    dataframes_ca1[key] = df_ca1

#%% 
#denoising the results dataset by using the drop rows on conditions function. Gets rid of any endpoint lower than two and any cells whose longest shortes path is lower than the 90th quantile.
processed_dfs = {}
for key, df in dataframes.items():
    if "Results" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df.columns.tolist()}")
        # Ensure the column names exist in the DataFrame
        required_columns = ["ID", "# Branches", "# End-point voxels", "Average Branch Length", "Longest Shortest Path"]
        if all(col in df.columns for col in required_columns):
           new_df = df[required_columns]
           #apply each condition to drop rows from the dataframes
           df_dropped = drop_rows_on_conditions(new_df, conditions)
           #append the processed dataframes to the list
           processed_dfs[key] = df_dropped
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")

# Update the original dictionary with the processed DataFrames
dataframes.update(processed_dfs)
        
processed_dfs_ca1 = {}
for key, df_ca1 in dataframes_ca1.items():
    if "Results" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df_ca1.columns.tolist()}")
        # Ensure the column names exist in the DataFrame
        required_columns = ["ID", "# Branches", "# End-point voxels", "Average Branch Length", "Longest Shortest Path"]
        if all(col in df_ca1.columns for col in required_columns):
           new_df_ca1 = df_ca1[required_columns]
           #apply each condition to drop rows from the dataframes
           df_dropped_ca1 = drop_rows_on_conditions(new_df_ca1, conditions_ca1)
           #append the processed dataframes to the list
           processed_dfs_ca1[key] = df_dropped_ca1
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")
            
dataframes_ca1.update(processed_dfs_ca1)
#%%
#tidying the dataframes by getting rid of any unwanted columns
#starts a dictionary where to store the tidied dataframes
processed_m = {}
#starts a loop which iterates through all of the dataframes in the dictionary
for key, df in dataframes.items():
    #only selects the m and ab coordinates dataframes
    if "M_coordinates" or "AB_coordinates" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df.columns.tolist()}")
        #selects the required columns and ensures the column names exist in the DataFrame
        required_columns_m = ["ID", "Area", "X", "Y", "Feret"]
        if all(col in df.columns for col in required_columns_m):
           df_m = df[required_columns_m]
           processed_m[key] = df_m
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")
#updates the original dictionary of dataframes with the tidied ones            
dataframes.update(processed_m)
            
for key, df in dataframes.items():
    if "M_coordinates" in key:
        df['MicrogliaNumber'] = len(df)
        df_n = df[required_columns_m + ['MicrogliaNumber']]
        dataframes[key] = df_n


processed_m_ca1 = {}
for key, df_ca1 in dataframes_ca1.items():
    if "M_coordinates" or "AB_coordinates" in key:
        #could use a print statement to verify if the code is working 
        print(f"Processing Dataframe {key} with columns: {df_ca1.columns.tolist()}")
        # Ensure the column names exist in the DataFrame
        required_columns_m = ["ID", "Area", "X", "Y", "Feret"]
        if all(col in df_ca1.columns for col in required_columns_m):
           df_m_ca1 = df_ca1[required_columns_m]
           processed_m_ca1[key] = df_m_ca1
        else:
            print(f"Skipping DataFrame {key} because required columns are missing")
            
dataframes_ca1.update(processed_m_ca1)

for key, df_ca1 in dataframes_ca1.items():
    if "M_coordinates" in key:
        df_ca1['MicrogliaNumber'] = len(df_ca1)
        df_n_ca1 = df_ca1[required_columns_m + ['MicrogliaNumber']]
        dataframes_ca1[key] = df_n_ca1
        
#%%
filtered_dfs = {}
keys_to_delete = []

for key, df in dataframes.items():
    if "Branchinfo" in key:
        filtered_dfs[key] = df
        keys_to_delete.append(key)

# Delete the collected keys from the original dictionary
for key in keys_to_delete:
    del dataframes[key]

# Print the filtered DataFrames keys to verify
print("Filtered DataFrames keys:", filtered_dfs.keys())
       
filtered_dfs_ca1 = {}
keys_to_delete_ca1 = []

for key, df_ca1 in dataframes_ca1.items():
    if "Branchinfo" in key:
        filtered_dfs_ca1[key] = df_ca1
        keys_to_delete_ca1.append(key)

# Delete the collected keys from the original dictionary
for key in keys_to_delete_ca1:
    del dataframes_ca1[key]

# Print the filtered DataFrames keys to verify
print("Filtered DataFrames keys:", filtered_dfs_ca1.keys())
#%%
dataframes_means = {}
for key, df in dataframes.items():
    df_mean = df.iloc[:, 1:].mean()
    dataframes_means[key] = df_mean

dataframes_ca1_means = {}
for key, df_ca1 in dataframes_ca1.items():
    df_mean_ca1 = df_ca1.iloc[:, 1:].mean()
    dataframes_ca1_means[key] = df_mean_ca1   
#%%
dataframes_combined = pd.concat(dataframes_means, axis = 0)
#output_path_rsc = r"C:\Users\silve\Desktop\RSC\dataframes.xlsx"
#dataframes_combined.to_excel(output_path_rsc, index=True)

dataframes_ca1_combined = pd.concat(dataframes_ca1_means, axis = 0)
#output_path_ca1 = r"C:\Users\silve\Desktop\CA1\dataframes_ca1.xlsx"
#dataframes_ca1_combined.to_excel(output_path_ca1, index=True)
#%%
# Create the MultiIndex Series
index = pd.MultiIndex.from_tuples(dataframes_combined.keys(), names=['MouseID', 'Parameter'])
dataframes_combined_wide = pd.Series(dataframes_combined, index=index)

# Convert the Series to a DataFrame
df_wide = dataframes_combined_wide.unstack(level="MouseID")
print(df_wide)

# Create the MultiIndex Series
index_ca1 = pd.MultiIndex.from_tuples(dataframes_ca1_combined.keys(), names=['MouseID', 'Parameter'])
dataframes_combined_wide_ca1 = pd.Series(dataframes_ca1_combined, index=index_ca1)

# Convert the Series to a DataFrame
df_wide_ca1 = dataframes_combined_wide_ca1.unstack(level="MouseID")
print(df_wide_ca1)

#%%
#RSC
#identify the unique base names
base_names = set(col.split("_")[0] for col in df_wide.columns)

#create a dictionary to store the new columns
new_columns = {}

#loop trhough the base names to create new columns by averaging left and right values
for base in base_names:
    left_cols = [col for col in df_wide.columns if col.startswith(base + "_") and "Left" in col]
    right_cols = [col for col in df_wide.columns if col.startswith(base + "_") and "Right" in col]
    for left_col, right_col in zip(left_cols, right_cols):
        new_col_name = left_col.replace("Left", "")
        new_columns[new_col_name] = (df_wide[left_col] + df_wide[right_col]) / 2

#Add the new columns to the df
for col_name, col_data in new_columns.items():
    df_wide[col_name] = col_data

#drop the old left and right columns
df_wide = df_wide.drop(columns = [col for col in df_wide.columns if "Left" in col or "Right" in col])

#CA1
#identify the unique base names
base_names_ca1 = set(col.split("_")[0] for col in df_wide_ca1.columns)

#create a dictionary to store the new columns
new_columns_ca1 = {}

#loop trhough the base names to create new columns by averaging left and right values
for base_ca1 in base_names_ca1:
    left_cols_ca1 = [col for col in df_wide_ca1.columns if col.startswith(base_ca1 + "_") and "Left" in col]
    right_cols_ca1 = [col for col in df_wide_ca1.columns if col.startswith(base_ca1 + "_") and "Right" in col]
    
    for left_col_ca1, right_col_ca1 in zip(left_cols_ca1, right_cols_ca1):
        new_col_name_ca1 = left_col_ca1.replace("Left", "")
        new_columns_ca1[new_col_name_ca1] = (df_wide_ca1[left_col_ca1] + df_wide_ca1[right_col_ca1]) / 2

#Add the new columns to the df
for col_name_ca1, col_data_ca1 in new_columns_ca1.items():
    df_wide_ca1[col_name_ca1] = col_data_ca1

#drop the old left and right columns
df_wide_ca1 = df_wide_ca1.drop(columns = [col for col in df_wide_ca1.columns if "Left" in col or "Right" in col])
#%%
# Transpose the DataFrame
df_transposed = df_wide.T

# Reset the index to make the MouseID a column again
df_transposed.reset_index(inplace=True)

# Rename the columns if necessary (optional)
df_transposed.rename(columns={'index': 'MouseID'}, inplace=True)

print(df_transposed)

output_path_rsc = r"C:\Users\silve\Desktop\RSC\dataframes.xlsx"
df_transposed.to_excel(output_path_rsc, index=True)

# Transpose the DataFrame
df_transposed_ca1 = df_wide_ca1.T

# Reset the index to make the MouseID a column again
df_transposed_ca1.reset_index(inplace=True)

# Rename the columns if necessary (optional)
df_transposed_ca1.rename(columns={'index': 'MouseID'}, inplace=True)

print(df_transposed_ca1)

output_path_ca1 = r"C:\Users\silve\Desktop\CA1\dataframes_ca1.xlsx"
df_transposed_ca1.to_excel(output_path_ca1, index=True)
#%%
df_transposed.rename(columns = {'MouseID' : 'MouseID_Type'}, inplace=True)
#extract mouse id and type
df_transposed[['MouseID', 'Type']] = df_transposed['MouseID_Type'].str.extract(r'^(M\w+)(?:_(M_coordinates|AB_coordinates|Results).*)$')
#drop the original mouseid:type column
df_transposed.drop(columns = ['MouseID_Type'], inplace = True)

df_melted = df_transposed.melt(id_vars = ['MouseID', 'Type'], var_name = 'Parameter', value_vars = ['Area', 'X', 'Y', 'Feret'])
#pivot to get the format with separate columns for microglia and abeta
df_pivoted = df_melted.pivot_table(index = ['MouseID', 'Parameter'], columns = 'Type', values = 'value').reset_index()

# Flatten the columns and rename them appropriately
df_pivoted.columns = [f"{col[1]}_{col[0]}" if col[1] else col[0] for col in df_pivoted.columns]
df_pivoted.rename(columns={'M_Value': 'microglia', 'AB_Value': 'abeta'}, inplace=True)

# Pivot again to have parameters as columns with microglia and abeta differentiated
df_final = df_pivoted.pivot(index='MouseID', columns='Parameter', values=['microglia', 'abeta'])

# Flatten the columns and rename them appropriately
df_final.columns = [f'{param}_{measure}' for measure, param in df_final.columns]
df_final.reset_index(inplace=True)

# Print the final DataFrame
print(df_final)
#%%

