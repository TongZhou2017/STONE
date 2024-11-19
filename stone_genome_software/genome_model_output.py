#!/usr/bin/env python3

import os
import sys
import numpy as np
import pandas as pd
import pickle
import argparse
from sklearn.preprocessing import MinMaxScaler

# Load the model (modified to take user input for the model file path)
def load_model(model_path):
    """Load the AutoML model from the specified path."""
    return pickle.load(open(model_path, 'rb'))

# Define percentile normalization function
def normalize_by_percentile(series, lower_percentile=5, upper_percentile=95):
    """Normalize a pandas Series by clipping at specified percentiles and scaling it."""
    lower = np.percentile(series.dropna(), lower_percentile)
    upper = np.percentile(series.dropna(), upper_percentile)
    series_clipped = np.clip(series, lower, upper)
    series_normalized = (series_clipped - lower) / (upper - lower)
    return series_normalized

def percentile_normalize(series, lower_percentile=5, upper_percentile=95, outlier_value=-999):
    """
    Normalize a pandas Series by clipping at specified percentiles and scaling it, excluding outlier values.
    Handles NaN values.
    """
    # Filter out outlier values
    filtered_series = series[series != outlier_value]
    
    # Drop NaN values
    filtered_series = filtered_series.dropna()
    
    # Calculate the percentile boundaries
    lower = np.percentile(filtered_series, lower_percentile)
    upper = np.percentile(filtered_series, upper_percentile)
    
    # Avoid division by zero if the upper and lower boundaries are the same
    if upper == lower:
        return pd.Series(np.zeros(len(series)), index=series.index)
    
    # Clip and normalize the series
    series_clipped = np.clip(series, lower, upper)
    series_normalized = (series_clipped - lower) / (upper - lower)
    
    return series_normalized

# Define value remapping function
def remap_values(values):
    """Remap values according to the method described by Zarringhalam et al. (2012)."""
    remapped = np.zeros_like(values)
    remapped[values < 0.25] = values[values < 0.25] * 0.35 / 0.25
    remapped[(values >= 0.25) & (values < 0.3)] = 0.35 + (values[(values >= 0.25) & (values < 0.3)] - 0.25) * 0.2 / 0.05
    remapped[(values >= 0.3) & (values < 0.7)] = 0.55 + (values[(values >= 0.3) & (values < 0.7)] - 0.3) * 0.3 / 0.4
    remapped[values >= 0.7] = 0.85 + (values[values >= 0.7] - 0.7) * 0.15 / 0.3
    return remapped

# Define data processing function
def process_testdata(file_path, depth_threshold=50, rf_mutation_Count=0, pipe_truncation_count=0, rate_mut=0.25, rate_stop=1, acc_threshold=0, filtered_bases=None, filter_acc=None, custom_names=None):
    # Read the data
    df_original = pd.read_csv(file_path)

    # Create a copy of the original data to preserve all original position information
    df_original_positions = df_original[['ChrID', 'pipe_truncation_Strand', 'pipe_truncation_ChrPos', 'rf_mutation_Base']].copy()

    # Create a copy for processing
    df_test = df_original.copy()

    # Drop specified columns
    columns_to_drop = ['ChrID', 'rf_mutation_Base', 'pipe_truncation_Strand', 'rf_mutation_AC', 'rf_mutation_AG', 'rf_mutation_AT',
                       'rf_mutation_CA', 'rf_mutation_CT', 'rf_mutation_CG', 'rf_mutation_GA', 'rf_mutation_GT', 'rf_mutation_GC',
                       'rf_mutation_TA', 'rf_mutation_TC', 'rf_mutation_TG', 'pipe_truncation_BD', 'pipe_truncation_Base', 'pipe_truncation_ChrPos']

    # Calculate new columns
    df_test['rate_A'] = df_test['base_A'] / df_test['rf_mutation_Depth']
    df_test['rate_T'] = df_test['base_T'] / df_test['rf_mutation_Depth']
    df_test['rate_C'] = df_test['base_C'] / df_test['rf_mutation_Depth']
    df_test['rate_G'] = df_test['base_G'] / df_test['rf_mutation_Depth']
    df_test['rate_stop'] = df_test['pipe_truncation_count'] / df_test['rf_mutation_Depth']
    df_test['rate_mut'] = df_test['rf_mutation_Count'] / df_test['rf_mutation_Depth']

    # Filter based on specified bases
    if filtered_bases:
        mask_base = df_test['rf_mutation_Base'].str.upper().apply(lambda x: any(base.upper() in x for base in filtered_bases))
        df_test = df_test[~mask_base.fillna(False)]

    # Apply depth threshold filter
    df_test = df_test[df_test['rf_mutation_Depth'] >= depth_threshold]
    df_test = df_test[df_test['rate_stop'] <= rate_stop]
    df_test = df_test[df_test['rate_mut'] < rate_mut]
    df_test_unfiltered = df_test.copy()  # Keep unfiltered version for reference
    df_test = df_test[df_test['rf_mutation_Count'] >= rf_mutation_Count]
    df_test = df_test[df_test['pipe_truncation_count'] >= pipe_truncation_count]

    # Normalize columns using percentile normalization
    columns_to_normalize = ['rate_mut', 'rate_stop', 'rate_A', 'rate_T', 'rate_C', 'rate_G']
    for col in columns_to_normalize:
        df_test[col] = normalize_by_percentile(df_test[col])
        df_test_unfiltered[col] = normalize_by_percentile(df_test_unfiltered[col])

    # Remap normalized values
    for col in columns_to_normalize:
        df_test[col] = remap_values(df_test[col])
        df_test_unfiltered[col] = remap_values(df_test_unfiltered[col])

    # Filter 'acc' column if required
    if filter_acc and 'acc' in df_test.columns:
        df_test = df_test[(df_test['acc'] >= acc_threshold) | df_test['acc'].isnull()]
        df_test = df_test.drop('acc', axis=1)  # Drop the filtered 'acc' column

    # Drop unnecessary columns
    df_test = df_test.drop(columns=columns_to_drop, errors='ignore')
    df_test_cleaned = df_test

    # Get custom column names if provided
    if custom_names:
        df_test_cleaned.columns = custom_names

    # Create feature matrix X and label vector y
    if 'modified_string' in df_test_cleaned.columns:
        X = df_test_cleaned.drop(['modified_string'], axis=1)
        y = df_test_cleaned['modified_string']
    else:
        X = df_test_cleaned
        y = None
    y_mut_score = df_test_cleaned['rate_mut']
    y_stop_score = df_test_cleaned['rate_stop']

    return X, y, y_mut_score, y_stop_score, df_test_unfiltered, df_original_positions

# Define function to process and save the data
def process_and_save(file_path, output_file_path, model):
    # Call the function to process data
    X, y, y_mut_score, y_stop_score, df_cleaned, df_original_positions = process_testdata(
        file_path, depth_threshold, rf_mutation_Count, pipe_truncation_count, rate_mut, rate_stop,
        acc_threshold, filtered_bases=filtered_base, filter_acc=filter_acc)

    # Use the model to make predictions
    proba_results = model.predict_proba(X)
    proba_results_second_column = proba_results[:, 1]

    # Add the prediction, mut_score, and stop_score to the data
    df_cleaned['predict'] = proba_results_second_column
    df_cleaned['mut_score'] = y_mut_score
    df_cleaned['stop_score'] = y_stop_score
    
    # Normalize 'predict' column by percentile within each 'ChrID' group
    df_cleaned['norm_model'] = df_cleaned.groupby('ChrID')['predict'].transform(
        lambda x: percentile_normalize(x, lower_percentile=5, upper_percentile=95, outlier_value=-999)
    )

    # Merge the cleaned data with the original position data and fill missing positions with NaN
    df_output = df_original_positions.merge(
        df_cleaned[['ChrID', 'pipe_truncation_Strand', 'pipe_truncation_ChrPos', 'predict', 'mut_score', 'stop_score', 'norm_model']],
        how='left', 
        on=['ChrID', 'pipe_truncation_Strand', 'pipe_truncation_ChrPos'])
    
    # Save the result to a CSV file
    df_output.to_csv(output_file_path, index=False)
    print(f"Results saved to {output_file_path}")

# Define default parameters
depth_threshold = 10
rf_mutation_Count = 0
pipe_truncation_count = 0
rate_mut = 0.25
rate_stop = 1
acc_threshold = 0
filtered_base = False
filter_acc = True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process and predict data using SHAPE model.')
    parser.add_argument('-i', '--input_folder', required=True, help='Path to the folder containing input CSV files.')
    parser.add_argument('-o', '--output_folder', required=True, help='Path to the folder where results will be saved.')
    parser.add_argument('-m', '--model_path', required=True, help='Path to the SHAPE model file (e.g., .sav file).')

    args = parser.parse_args()

    input_folder = args.input_folder
    output_folder = args.output_folder
    model_path = args.model_path

    # Load the model
    automl = load_model(model_path)

    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Process all CSV files in the input folder
    for filename in os.listdir(input_folder):
        if filename.endswith(".csv"):
            input_file_path = os.path.join(input_folder, filename)
            output_file_path = os.path.join(output_folder, filename)  # Keep the same output file name as input
            process_and_save(input_file_path, output_file_path, automl)
