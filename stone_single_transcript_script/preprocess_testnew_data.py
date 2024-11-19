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

def process_testdata(file_path, depth_threshold=50, rf_mutation_Count=0, pipe_truncation_count=0, rate_mut=0.25, rate_stop=1, acc_threshold=0, filtered_bases=None, filter_acc=None, custom_names=None):

    df_original = pd.read_csv(file_path)
    df_test = df_original.copy()
    columns_to_drop = ['ChrID', 'rf_mutation_Base', 'pipe_truncation_Strand','rf_mutation_AC','rf_mutation_AG','rf_mutation_AT','rf_mutation_CA','rf_mutation_CT','rf_mutation_CG','rf_mutation_GA','rf_mutation_GT','rf_mutation_GC','rf_mutation_TA','rf_mutation_TC','rf_mutation_TG','pipe_truncation_BD','pipe_truncation_Base','pipe_truncation_ChrPos']
    # Calculate a new column
    df_test['rate_A'] = df_test['base_A'] / df_test['rf_mutation_Depth']
    df_test['rate_T'] = df_test['base_T'] / df_test['rf_mutation_Depth']
    df_test['rate_C'] = df_test['base_C'] / df_test['rf_mutation_Depth']
    df_test['rate_G'] = df_test['base_G'] / df_test['rf_mutation_Depth']
    df_test['rate_stop'] = df_test['pipe_truncation_count'] / df_test['rf_mutation_Depth']
    df_test['rate_mut'] = df_test['rf_mutation_Count'] / df_test['rf_mutation_Depth']
  
    # df_test['rate_A'] = df_test['base_A'] / df_test['depth']
    # df_test['rate_T'] = df_test['base_T'] / df_test['depth']
    # df_test['rate_C'] = df_test['base_C'] / df_test['depth']
    # df_test['rate_G'] = df_test['base_G'] / df_test['depth']

    # df_test['rate_stop'] = df_test['pipe_truncation_count'] / df_test['depth']
    # df_test['rate_mut'] = df_test['rf_mutation_Count'] / df_test['depth']

    if filtered_bases:
        mask_base = df_test['rf_mutation_Base'].str.upper().apply(lambda x: any(base.upper() in x for base in filtered_bases))
        df_test = df_test[~mask_base.fillna(False)]

    df_test = df_test[df_test['rf_mutation_Depth'] >= depth_threshold]
    df_test = df_test[df_test['rate_stop'] <= rate_stop]
    df_test = df_test[df_test['rate_mut'] < rate_mut]
    df_test1 = df_test.copy()
    df_test = df_test[df_test['rf_mutation_Count'] >= rf_mutation_Count]
    df_test = df_test[df_test['pipe_truncation_count'] >= pipe_truncation_count]

    columns_to_normalize = ['rate_mut', 'rate_stop', 'rate_A', 'rate_T', 'rate_C', 'rate_G']
    for col in columns_to_normalize:
        df_test[col] = normalize_by_percentile(df_test[col])
        df_test1[col] = normalize_by_percentile(df_test1[col])

    for col in columns_to_normalize:
        df_test[col] = remap_values(df_test[col])
        df_test1[col] = remap_values(df_test1[col])


    if filter_acc and 'acc' in df_test.columns:
        df_test = df_test[(df_test['acc'] >= acc_threshold) | df_test['acc'].isnull()]
        df_test = df_test.drop('acc', axis=1)  

    df_test = df_test.drop(columns=columns_to_drop)
    df_test_cleaned = df_test
    if custom_names:
        df_test_cleaned.columns = custom_names

    # Create feature matrix X and label vector y
    X = df_test_cleaned.drop(['modified_string'], axis=1)
    y = df_test_cleaned['modified_string']
    y_mut_score = df_test_cleaned['rate_mut']
    y_stop_score = df_test_cleaned['rate_stop']

    return X, y, y_mut_score, y_stop_score, df_test1
