import numpy as np
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import sklearn
def calculate_auc(X, y, automl, df_original, rf_mutation_Count1, pipe_truncation_count1,title='ROC Curve',save_path=None):
    # Make predictions using the model.
    proba_results_XG = automl.predict_proba(X)
    proba_results_second_column = proba_results_XG[:, 1]
    
    # Add the prediction results 'merge_result' to df_original.
    df_original['merge_result'] = np.nan  # Initialize the merge_result column.
    df_original.loc[X.index, 'merge_result'] = proba_results_second_column
    # df_original['merge_result'] = proba_results_second_column

    # Fill the merge_result column with 0 for rows that meet the condition.
    # condition = (df_original['rf_mutation_Count'] == rf_mutation_Count1) | (df_original['pipe_truncation_count'] == pipe_truncation_count1)
    # df_original.loc[condition, 'merge_result'] = 0

    # Clean the NaN values in the true labels and extract the corresponding prediction results.
    df_cleaned = df_original.dropna(subset=['modified_string', 'merge_result'])
    y_label = df_cleaned['modified_string']
    y_model = df_cleaned['merge_result']

    y_mut_score_row =df_cleaned['rate_mut']
    y_stop_score_row=df_cleaned['rate_stop']
    position=df_cleaned['pipe_truncation_ChrPos']
    
    # calculate the Area Under the Curve (AUC)
    auc_model = roc_auc_score(y_label, y_model)
    auc_stop =roc_auc_score(y_label, y_stop_score_row)
    auc_mut=roc_auc_score(y_label, y_mut_score_row)
    N=len(y_label)

    fpr_model, tpr_model, _ = sklearn.metrics.roc_curve(y_label, y_model)
    roc_auc_model = auc(fpr_model, tpr_model)
    N_model=len(y_label)
    fpr_mut, tpr_mut, _ = sklearn.metrics.roc_curve(y_label, y_mut_score_row)
    roc_auc_mut = auc(fpr_mut, tpr_mut)
    N_mut=len(y_label)
    fpr_stop, tpr_stop, _ = sklearn.metrics.roc_curve(y_label, y_stop_score_row)
    roc_auc_stop = auc(fpr_stop, tpr_stop)
    N_stop=len(y_label)
    # Plot ROC curves
    plt.figure(figsize=(7, 7))
    plt.plot(fpr_model, tpr_model, '-', color='seagreen', label=f"AUC_model={roc_auc_model:.3f},N={N_model}")
    plt.plot(fpr_mut, tpr_mut, '-', color='#E4A689', label=f"AUC_mut={roc_auc_mut:.3f},N={N_mut}")
    plt.plot(fpr_stop, tpr_stop, '-', color='purple', label=f"AUC_stop={roc_auc_stop:.3f},N={N_stop}")
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='grey', alpha=0.8)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc='lower right')
    plt.grid(False)
    if save_path:
        plt.savefig(save_path, format='pdf')
    plt.show()
    return auc_model,auc_stop,auc_mut,position,y_model,y_mut_score_row,y_stop_score_row,plt


