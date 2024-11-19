import pandas as pd
import numpy as np

def normalize_by_percentile(series, lower_percentile=5, upper_percentile=95):
    """Normalize a pandas Series by clipping at specified percentiles and scaling."""
    lower = np.percentile(series.dropna(), lower_percentile)
    upper = np.percentile(series.dropna(), upper_percentile)
    series_clipped = np.clip(series, lower, upper)
    series_normalized = (series_clipped - lower) / (upper - lower)
    return series_normalized

def remap_values(values):
    """Remap values according to Zarringhalam et al., 2012."""
    remapped = np.zeros_like(values)
    remapped[values < 0.25] = values[values < 0.25] * 0.35 / 0.25
    remapped[(values >= 0.25) & (values < 0.3)] = 0.35 + (values[(values >= 0.25) & (values < 0.3)] - 0.25) * 0.2 / 0.05
    remapped[(values >= 0.3) & (values < 0.7)] = 0.55 + (values[(values >= 0.3) & (values < 0.7)] - 0.3) * 0.3 / 0.4
    remapped[values >= 0.7] = 0.85 + (values[values >= 0.7] - 0.7) * 0.15 / 0.3
    return remapped

def process_testdata(df_T_path, df_D_path, depth_threshold=50, rf_mutation_Count=0, pipe_truncation_count=0, rate_mut=0.25, rate_stop=1, acc_threshold=0, filtered_bases=None, filter_acc=None, custom_names=None):
    # Load data
    df_original_T = pd.read_csv(df_T_path)
    df_original_D = pd.read_csv(df_D_path)
    # Create a copy for processing
    df_test_T = df_original_T.copy()
    df_test_D = df_original_D.copy()

    # Delete specified column
    columns_to_drop = ['ChrID', 'rf_mutation_Base', 'pipe_truncation_Strand', 'rf_mutation_AC', 'rf_mutation_AG', 'rf_mutation_AT', 'rf_mutation_CA', 'rf_mutation_CT', 'rf_mutation_CG', 'rf_mutation_GA', 'rf_mutation_GT', 'rf_mutation_GC', 'rf_mutation_TA', 'rf_mutation_TC', 'rf_mutation_TG', 'pipe_truncation_BD', 'pipe_truncation_Base', 'pipe_truncation_ChrPos']
    
    # Calculate a new column
    df_test_T['rate_A'] = df_test_T['base_A'] / df_test_T['rf_mutation_Depth']
    df_test_T['rate_T'] = df_test_T['base_T'] / df_test_T['rf_mutation_Depth']
    df_test_T['rate_C'] = df_test_T['base_C'] / df_test_T['rf_mutation_Depth']
    df_test_T['rate_G'] = df_test_T['base_G'] / df_test_T['rf_mutation_Depth']
    df_test_T['rate_stop'] = df_test_T['pipe_truncation_count'] / df_test_T['rf_mutation_Depth']
    df_test_T['rate_mut'] = df_test_T['rf_mutation_Count'] / df_test_T['rf_mutation_Depth']

    df_test_D['rate_A'] = df_test_D['base_A'] / df_test_D['rf_mutation_Depth']
    df_test_D['rate_T'] = df_test_D['base_T'] / df_test_D['rf_mutation_Depth']
    df_test_D['rate_C'] = df_test_D['base_C'] / df_test_D['rf_mutation_Depth']
    df_test_D['rate_G'] = df_test_D['base_G'] / df_test_D['rf_mutation_Depth']
    df_test_D['rate_stop'] = df_test_D['pipe_truncation_count'] / df_test_D['rf_mutation_Depth']
    df_test_D['rate_mut'] = df_test_D['rf_mutation_Count'] / df_test_D['rf_mutation_Depth']

    # Specify columns to process
    columns_to_adjust = ['rate_A', 'rate_T', 'rate_C', 'rate_G', 'rate_stop', 'rate_mut']

    # Perform a subtraction operation on each column and set values less than zero to zero
    for col in columns_to_adjust:
        df_diff = df_test_T[col] - df_test_D[col]
        df_test_T[col] = df_diff.apply(lambda x: max(x, 0))
    
    # Select the input bases
    if filtered_bases:
        mask_base = df_test_T['rf_mutation_Base'].str.upper().apply(lambda x: any(base.upper() in x for base in filtered_bases))
        df_test_T = df_test_T[~mask_base.fillna(False)]

    # Threshold filtering
    df_test_T = df_test_T[df_test_T['rf_mutation_Depth'] >= depth_threshold]
    df_test_T = df_test_T[df_test_T['rate_stop'] <= rate_stop]
    df_test_T = df_test_T[df_test_T['rate_mut'] < rate_mut]

    df_test_T1 = df_test_T.copy()

    df_test_T = df_test_T[df_test_T['rf_mutation_Count'] >= rf_mutation_Count]
    df_test_T = df_test_T[df_test_T['pipe_truncation_count'] >= pipe_truncation_count]

    # Normalize by percentile and remap
    for col in columns_to_adjust:
        df_test_T[col] = normalize_by_percentile(df_test_T[col])
        df_test_T[col] = remap_values(df_test_T[col])
        df_test_T1[col] = normalize_by_percentile(df_test_T1[col])
        df_test_T1[col] = remap_values(df_test_T1[col])

    # Filter the 'acc' column
    if filter_acc and 'acc' in df_test_T.columns:
        df_test_T = df_test_T[(df_test_T['acc'] >= acc_threshold) | df_test_T['acc'].isnull()]
        df_test_T = df_test_T.drop('acc', axis=1) #Delete the filtered 'acc' column
    
    df_test_T_final = df_test_T.drop(columns=columns_to_drop)
    # Clean data
    df_test_cleaned = df_test_T_final.replace([np.inf, -np.inf], np.nan).dropna()

    # Retrieve custom column names
    if custom_names:
        df_test_cleaned.columns = custom_names

    # Create feature matrix X and label vector y
    X = df_test_cleaned.drop(['modified_string'], axis=1)
    y = df_test_cleaned['modified_string']
    y_mut_score = df_test_cleaned['rate_mut']
    y_stop_score = df_test_cleaned['rate_stop']

    return X, y, y_mut_score, y_stop_score, df_test_T1
