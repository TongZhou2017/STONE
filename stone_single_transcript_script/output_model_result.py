import pandas as pd
def merge_and_save_data(chr_position, proba_results_second_column, y_mut_score, y_stop_score, sequence_length, csv_path):
    # Merge Data
    y_list = proba_results_second_column.tolist()
    y_mut = y_mut_score.tolist()
    y_stop = y_stop_score.tolist()
    
    df_merged = pd.DataFrame({'ChrPos': chr_position, 'predict': y_list, 'mut_score': y_mut, 'stop_score': y_stop})

    full_sequence = pd.DataFrame({'ChrPos': range(1, sequence_length + 1)})

    merged_data = pd.merge(full_sequence, df_merged, on='ChrPos', how='left')

    # Convert column type to string
    merged_data['predict'] = merged_data['predict'].astype(str)
    merged_data['mut_score'] = merged_data['mut_score'].astype(str)
    merged_data['stop_score'] = merged_data['stop_score'].astype(str)

    # Fill NaN values
    merged_data['predict'] = merged_data['predict'].fillna("NULL")
    merged_data['mut_score'] = merged_data['mut_score'].fillna("NULL")
    merged_data['stop_score'] = merged_data['stop_score'].fillna("NULL")

    # Write to CSV file
    merged_data.to_csv(csv_path, index=False)

    # print(f"CSV file '{csv_path}' has been created.")
