# -----------------------------------------------------------------------------
# 02_train_growth_stage_models.py
# -----------------------------------------------------------------------------
#
# Author:      Thomas Martini Jørgensen <tmjq@dtu.dk>
#
# Copyright (c) 2026, Thomas Martini Jørgensen / DTU Compute, Technical University of Denmark
#
# This script is part of the code accompanying the paper:
# "Identification of robust bacterial biomarkers in the wheat root microbiome for predicting soil legacy effects and crop development "
# by Frederik Bak, Thomas Martini Jørgensen, Inês Nunes, Veronika Hansen, and  Mette Haubjerg Nicolaisen
#
# GitHub Repository: https://github.com/BakDK/Soil_legacy_2026 
#
# This script performs the training of the random forest models for the growth stage analysis and performs Shapley Analysis 
# to interpret the model predictions. 

# This code is licensed under the MIT License. See the LICENSE file in the
# root directory of this repository for the full license text.
#
# -----------------------------------------------------------------------------

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import shap
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import confusion_matrix, accuracy_score, cohen_kappa_score, f1_score, roc_auc_score
from sklearn.preprocessing import label_binarize # Needed for multi-class AUC
import seaborn as sns
import os
import sys # Added for command-line arguments

# ===================================================================
# --- CORE ANALYSIS FUNCTION  --
# ===================================================================
def run_full_analysis_rf(train_csv_path, test_csv_path, output_dir, n_runs=20, seed_start=123):
    print("--- Reading data ---")
    try:
        train_df = pd.read_csv(train_csv_path)
        test_df = pd.read_csv(test_csv_path)
    except FileNotFoundError as e:
        print(f"FATAL ERROR: Data file not found: {e}. Halting execution.")
        sys.exit(1)
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # --- Data Preparation  ---
    target_col = 'response_var'
    feature_cols = [col for col in train_df.columns if col != target_col]
    train_df['target_code'] = train_df[target_col].astype('category').cat.codes
    # Ensure test set categories match training set categories
    train_categories = train_df[target_col].astype('category').cat.categories
    test_df[target_col] = pd.Categorical(test_df[target_col], categories=train_categories, ordered=True)
    test_df['target_code'] = test_df[target_col].cat.codes
    
    class_names = list(train_df[target_col].astype('category').cat.categories)
    X_train, y_train = train_df[feature_cols], train_df['target_code']
    X_test, y_test = test_df[feature_cols], test_df['target_code']
    print(f"Starting analysis for {n_runs} runs with RandomForest...")

    # --- Main Loop ) ---
    all_targeted_results_list = []
    last_run_shap_object = None
    all_general_summary_rows = []
    all_confusion_matrices = []
    
    # Initialize a list to store performance metrics from each run
    all_performance_metrics = []

    for i in range(n_runs):
        # Main loop for training and evaluation 
        
        current_seed = seed_start + i
        print(f"\n--- Running model {i+1}/{n_runs} with seed {current_seed} ---")
        model = RandomForestClassifier(n_estimators=200, max_depth=10, max_features="sqrt", n_jobs=-1, random_state=current_seed)
        model.fit(X_train, y_train)
        explainer = shap.Explainer(model, X_train)
        shap_values_run = explainer(X_test)
        
        predictions = model.predict(X_test)
        
        # Calculate performance metrics for this run
        
        # For AUC, we need prediction probabilities
        y_pred_proba = model.predict_proba(X_test)
        
        # Binarize the true labels for multi-class AUC calculation
        y_test_binarized = label_binarize(y_test, classes=model.classes_)
        
        # One-vs-Rest AUC for each class
        per_class_auc = roc_auc_score(y_test_binarized, y_pred_proba, average=None, multi_class='ovr')
        
        run_metrics = {
            'run': i,
            'accuracy': accuracy_score(y_test, predictions),
            'kappa': cohen_kappa_score(y_test, predictions),
            'f1_macro': f1_score(y_test, predictions, average='macro'),
            'f1_weighted': f1_score(y_test, predictions, average='weighted'),
            'auc_macro': roc_auc_score(y_test_binarized, y_pred_proba, average='macro', multi_class='ovr'),
            'auc_weighted': roc_auc_score(y_test_binarized, y_pred_proba, average='weighted', multi_class='ovr')
        }
        
        # Add the per-class AUC scores to the dictionary
        for class_idx, class_name in enumerate(class_names):
            run_metrics[f'auc_{class_name}'] = per_class_auc[class_idx]
            
        all_performance_metrics.append(run_metrics)
        
        cm_run = confusion_matrix(y_test, predictions, labels=model.classes_)
        all_confusion_matrices.append(cm_run)

        mean_abs_shaps_this_run = np.mean(np.abs(shap_values_run.values), axis=0)
        
        for class_idx, class_name in enumerate(class_names):
            for feature_idx, feature_name in enumerate(feature_cols):
                all_general_summary_rows.append({
                    'Run': i, 'Feature': feature_name, 'Class': class_name,
                    'Mean_Abs_SHAP_per_run': mean_abs_shaps_this_run[feature_idx, class_idx]
                })

        for j, class_name in enumerate(class_names):
            tp_mask_series = (predictions == j) & (y_test == j)
            tp_mask_numpy = tp_mask_series.values
            num_tps_in_run = tp_mask_numpy.sum()

            if num_tps_in_run > 0:
                X_test_tps = X_test[tp_mask_numpy]
                X_test_non_tps = X_test[~tp_mask_numpy]
                shap_values_tps = shap_values_run[tp_mask_numpy, :, j]
                
                for k, feature_name in enumerate(feature_cols):
                    feature_abund_in_tps = X_test_tps[feature_name]
                    high_abund_mask = (feature_abund_in_tps > 0).values
                    num_abundant_tps = high_abund_mask.sum()
                    
                    if num_abundant_tps > 0:
                        feature_shap_in_tps = shap_values_tps.values[:, k]
                        sum_shaps_in_high_abund = feature_shap_in_tps[high_abund_mask].sum()
                        
                        if sum_shaps_in_high_abund > 0:
                            mean_abund_in_tps = feature_abund_in_tps.mean()
                            mean_abund_in_non_tps = X_test_non_tps[feature_name].mean()
                            fold_change = (mean_abund_in_tps + 1e-9) / (mean_abund_in_non_tps + 1e-9)
                            
                            potency_shap = sum_shaps_in_high_abund / num_abundant_tps
                            overall_impact_shap = sum_shaps_in_high_abund / num_tps_in_run
                            mean_abund_in_high_abund = feature_abund_in_tps[high_abund_mask].mean()
                            prop_tps = num_abundant_tps / num_tps_in_run
                            
                            all_targeted_results_list.append({
                                'Run': i, 'Feature': feature_name, 'Class': class_name,
                                'Potency_SHAP': potency_shap, 'Overall_Impact_SHAP': overall_impact_shap,
                                'Mean_Abund_in_High_Abund_TPs': mean_abund_in_high_abund,
                                'Proportion_of_TPs_with_Feature': prop_tps,
                                'Fold_Change_TP_vs_nonTP': fold_change
                            })
        
        if i == n_runs - 1:
            last_run_shap_object = shap_values_run

    # --- Aggregation and Saving 
    print("\n--- All runs complete. Aggregating and saving results. ---")
    
    # Aggregate and save the performance metrics
    if all_performance_metrics:
        # Convert list of dicts to a DataFrame
        metrics_df = pd.DataFrame(all_performance_metrics)
        
        # Calculate mean and std dev for each metric
        # .T transposes the result to get metrics as rows
        summary_df = metrics_df.drop(columns=['run']).agg(['mean', 'std']).T
        summary_df = summary_df.rename(columns={'mean': 'Mean', 'std': 'Std_Dev'})
        
        # Save to CSV
        output_path = os.path.join(output_dir, "rf_performance_metrics_summary.csv")
        summary_df.to_csv(output_path, index=True) # index=True to keep metric names as a column
        print(f"Performance metrics summary saved to {output_path}")
    
    if all_general_summary_rows:
        general_df_raw = pd.DataFrame(all_general_summary_rows)
        general_df_agg = general_df_raw.groupby(['Feature', 'Class']).agg(Mean_Abs_SHAP_Importance=('Mean_Abs_SHAP_per_run', 'mean'), SD_Abs_SHAP_Importance=('Mean_Abs_SHAP_per_run', 'std'), n_runs_feat_present=('Mean_Abs_SHAP_per_run', lambda x: (x > 1e-9).sum())).reset_index().sort_values(by=['Class', 'Mean_Abs_SHAP_Importance'], ascending=[True, False])
        
        # Save to output_dir
        output_path = os.path.join(output_dir, "rf_general_shap_summary.csv")
        general_df_agg.to_csv(output_path, index=False)
        print(f"Full RF general SHAP summary saved to {output_path}")
        
        top_15_df = general_df_agg.groupby('Class').head(15).reset_index(drop=True)
        output_path_top15 = os.path.join(output_dir, "rf_top_15_features.csv")
        top_15_df.to_csv(output_path_top15, index=False)
        print(f"Top 15 features summary for RF saved to {output_path_top15}")

    if all_targeted_results_list:
        targeted_df_raw = pd.DataFrame(all_targeted_results_list)
        targeted_df_agg = targeted_df_raw.groupby(['Feature', 'Class']).agg(Avg_Potency_SHAP=('Potency_SHAP', 'mean'), Avg_Overall_Impact_SHAP=('Overall_Impact_SHAP', 'mean'), Avg_Mean_Abund_in_High_Abund_TPs=('Mean_Abund_in_High_Abund_TPs', 'mean'), Avg_Proportion_of_TPs_with_Feature=('Proportion_of_TPs_with_Feature', 'mean'), Avg_Fold_Change_TP_vs_nonTP=('Fold_Change_TP_vs_nonTP', 'mean'), Num_Runs_Found=('Run', 'nunique')).reset_index().sort_values(by=['Class', 'Avg_Overall_Impact_SHAP'], ascending=[True, False])
        
        # Save to output_dir
        output_path = os.path.join(output_dir, "rf_targeted_shap_summary.csv")
        targeted_df_agg.to_csv(output_path, index=False)
        print(f"Aggregated RF targeted analysis summary saved to {output_path}")

    if all_confusion_matrices:
        all_cm_np = np.array(all_confusion_matrices)
        avg_cm = np.mean(all_cm_np, axis=0)
        std_cm = np.std(all_cm_np, axis=0)
        desired_class_order = ['early', 'middle', 'late']
        original_class_order = class_names 
        reorder_idx = [original_class_order.index(c) for c in desired_class_order]
        reordered_avg_cm = avg_cm[np.ix_(reorder_idx, reorder_idx)]
        reordered_std_cm = std_cm[np.ix_(reorder_idx, reorder_idx)]
        annot_labels = np.array([f"{avg:.1f}\n(±{std:.1f})" for avg, std in zip(reordered_avg_cm.ravel(), reordered_std_cm.ravel())]).reshape(reordered_avg_cm.shape)
        plt.figure(figsize=(10, 8))
        sns.heatmap(reordered_avg_cm, annot=annot_labels, fmt='', cmap='Greens', xticklabels=desired_class_order, yticklabels=desired_class_order, cbar_kws={'label': 'Average Count'})
        plt.title(f'Average Confusion Matrix over {n_runs} Runs - RandomForest\n(Mean ± SD)', fontsize=16)
        plt.ylabel('True Label')
        plt.xlabel('Predicted Label')
        
        # Save to output_dir
        output_path = os.path.join(output_dir, "rf_avg_confusion_matrix.png")
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close() # Use plt.close() instead of plt.show() for non-interactive scripts
        print(f"Confusion matrix plot saved to {output_path}")
        
    if last_run_shap_object is not None:
        print("\n--- Generating SHAP Beeswarm Plot (from last run, top features only) ---")
        current_seed = seed_start + n_runs - 1
        num_features_to_show = 15
        for i, class_name in enumerate(class_names):
            print(f"  - Plotting for class: {class_name}")
            class_shap_values = last_run_shap_object[:, :, i]
            mean_abs_shaps = np.mean(np.abs(class_shap_values.values), axis=0)
            top_indices = np.argsort(mean_abs_shaps)[-num_features_to_show:]
            top_features_shap_values = class_shap_values[:, top_indices]
            shap.plots.beeswarm(top_features_shap_values, max_display=num_features_to_show + 1, show=False)
            
            fig = plt.gcf()
            fig.suptitle(f"RandomForest SHAP Beeswarm for Class: {class_name}\n(Top {num_features_to_show} Features from Last Run, Seed {current_seed})", fontsize=16)
            fig.tight_layout(rect=[0, 0.03, 1, 0.95])
            
            # Save to output_dir
            filename = f"rf_shap_beeswarm_{class_name}_last_run.png"
            output_path = os.path.join(output_dir, filename)
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close() # Use plt.close() instead of plt.show()
        print("RF Beeswarm plots for the last run have been saved.")

# ===================================================================
# --- SCRIPT EXECUTION ORCHESTRATOR ---
# This block makes the script runnable from the command line.
# ===================================================================
if __name__ == '__main__':
    # USAGE:
    # python <script_name>.py <input_dir> <output_dir>
    #
    # EXAMPLE:
    # python src/python/02...py outputs/growth_stage_cv13-24/01_datasets outputs/growth_stage_cv13-24/02_models_and_shap
    #
    
    # --- 1. Argument Parsing and Validation ---
    if len(sys.argv) != 3:
        print("Error: Invalid number of arguments.")
        print("Usage: python <script_name>.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory not found at '{input_dir}'")
        sys.exit(1)
    
    # --- 2. Define Input File Paths ---
    # The R script produces consistent filenames
    train_file = os.path.join(input_dir, "growth_stage_train_data.csv")
    test_file = os.path.join(input_dir, "growth_stage_test_data.csv")
    
    # --- 3. Call the Main Analysis Function ---
    run_full_analysis_rf(
        train_csv_path=train_file,
        test_csv_path=test_file,
        output_dir=output_dir,
        n_runs=20 # This can remain hardcoded or be made an argument later if needed
    )
    
    print("\n\nAll Python analysis for this fold is complete!")