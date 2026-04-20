# -----------------------------------------------------------------------------
# 02_RandomForestTrainingAndShapleyAnalysisForLegacyData.py
# -----------------------------------------------------------------------------
#
# Author:      Thomas Martini Jørgensen <tmjq@dtu.dk>
#
# Copyright (c) 2026, Thomas Martini Jørgensen / DTU Compute, Technical University of Denmark
#
# This script is part of the code accompanying the paper:
# "Identification of robust bacterial biomarkers in the wheat root microbiome for predicting soil legacy and crop development "
# by Frederik Bak, Thomas Martini Jørgensen, Inês Nunes, Veronika Hansen, and  Mette Haubjerg Nicolaisen
#
# GitHub Repository: https://https://github.com/BakDK/Future_cropping
#
# This script performs the training of the random forest models for the legacy analysis and performs Shapley Analysis 
# to interpret the model predictions. It includes functionality to evaluate the models on various test sets from several time points.
#
# This code is licensed under the MIT License. See the LICENSE file in the
# root directory of this repository for the full license text.
#
# -----------------------------------------------------------------------------



import pandas as pd
from sklearn.ensemble import RandomForestClassifier
import shap
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report, cohen_kappa_score, roc_auc_score
from sklearn.preprocessing import label_binarize
import seaborn as sns
import os
import pickle
import sys  # Added for command-line arguments
import re   # Added for regular expressions


def train_and_save_models(train_csv_path, model_output_dir, n_runs=20, seed_start=123):
    print(f"\n{'='*20}\nPHASE 1: TRAINING AND SAVING MODELS\n{'='*20}")
    # The master script now creates this directory, but os.makedirs with exist_ok=True is safe
    os.makedirs(model_output_dir, exist_ok=True) 
    print("--- Reading and preparing training data ---")
    try:
        train_df = pd.read_csv(train_csv_path, index_col=0)
    except FileNotFoundError as e:
        print(f"FATAL ERROR: Training file not found: {e}. Halting execution.")
        sys.exit(1) # Use sys.exit() for a clean stop
        
    target_col = 'response_var'
    feature_cols = [col for col in train_df.columns if col != target_col]
    train_df[target_col] = train_df[target_col].astype('category')
    train_df['target_code'] = train_df[target_col].cat.codes
    X_train, y_train = train_df[feature_cols], train_df['target_code']
    for i in range(n_runs):
        current_seed = seed_start + i
        print(f"--- Training model {i+1}/{n_runs} with seed {current_seed} ---")
        model = RandomForestClassifier(n_estimators=200, max_depth=10, max_features="sqrt", n_jobs=-1, random_state=current_seed)
        model.fit(X_train, y_train)
     
        model_path = os.path.join(model_output_dir, f"model_run_{i}_seed_{current_seed}.pkl")
        
        explainer = shap.Explainer(model, X_train)
        model_payload = {'model': model, 'explainer': explainer}

        with open(model_path, 'wb') as f:
            pickle.dump(model_payload, f)
        print(f"   ... Saved to {model_path}")
    print("\n--- Model training phase complete. ---")



def evaluate_models_on_test_set(test_csv_path, model_dir, results_output_dir, output_prefix, class_names, n_runs=20, seed_start=123):
    print(f"\n{'='*20}\nPHASE 2: EVALUATING ON: {output_prefix}\n{'='*20}")
    # Ensure the output directory for this evaluation exists
    os.makedirs(results_output_dir, exist_ok=True)
    
    try:
        test_df = pd.read_csv(test_csv_path, index_col=0)
    except FileNotFoundError as e:
        print(f"Error: Test file not found: {e}. Skipping this analysis.")
        return
    target_col = 'response_var'
    test_df[target_col] = test_df[target_col].astype('category').cat.set_categories(class_names)
    test_df['target_code'] = test_df[target_col].cat.codes
    feature_cols = [col for col in test_df.columns if col not in [target_col, 'target_code']]
    X_test, y_test = test_df[feature_cols], test_df['target_code']

   
    all_general_summary_rows = []
    all_confusion_matrices = []
    all_performance_metrics = []
    all_targeted_results_list = []

    for i in range(n_runs):
        current_seed = seed_start + i
       
        model_path = os.path.join(model_dir, f"model_run_{i}_seed_{current_seed}.pkl")
    
        try:
            with open(model_path, 'rb') as f:
                model_payload = pickle.load(f)
            model = model_payload['model']
            explainer = model_payload['explainer']
        except FileNotFoundError:
            print(f"   ERROR: Model file not found at {model_path}. Skipping run {i+1}.")
            continue

        predictions = model.predict(X_test)
        prediction_probabilities = model.predict_proba(X_test)
        
        report = classification_report(y_test, predictions, output_dict=True, zero_division=0)
        
        run_metrics = {'accuracy': accuracy_score(y_test, predictions),'kappa': cohen_kappa_score(y_test, predictions),'f1_macro': report['macro avg']['f1-score'],'f1_weighted': report['weighted avg']['f1-score']}
        
        if len(np.unique(y_test)) > 1:
            run_metrics['auc_macro'] = roc_auc_score(y_test, prediction_probabilities, multi_class='ovr', average='macro')
            run_metrics['auc_weighted'] = roc_auc_score(y_test, prediction_probabilities, multi_class='ovr', average='weighted')
            y_test_binarized = label_binarize(y_test, classes=np.arange(len(class_names)))
            for idx, class_name in enumerate(class_names):
                if len(np.unique(y_test_binarized[:, idx])) > 1:
                    run_metrics[f'auc_{class_name}'] = roc_auc_score(y_test_binarized[:, idx], prediction_probabilities[:, idx])
                else:
                    run_metrics[f'auc_{class_name}'] = np.nan
        else:
            run_metrics['auc_macro'] = np.nan
            run_metrics['auc_weighted'] = np.nan
            for class_name in class_names:
                run_metrics[f'auc_{class_name}'] = np.nan
        
        all_performance_metrics.append(run_metrics)
        
        shap_values_run = explainer(X_test)
        all_confusion_matrices.append(confusion_matrix(y_test, predictions, labels=np.arange(len(class_names))))
        mean_abs_shaps_this_run = np.mean(np.abs(shap_values_run.values), axis=0)
        for class_idx, class_name in enumerate(class_names):
            for feature_idx, feature_name in enumerate(feature_cols):
                all_general_summary_rows.append({'Run': i, 'Feature': feature_name, 'Class': class_name, 'Mean_Abs_SHAP_per_run': mean_abs_shaps_this_run[feature_idx, class_idx]})
        for j, class_name in enumerate(class_names):
            tp_mask = (predictions == j) & (y_test == j)
            num_tps_in_run = tp_mask.sum()
            if num_tps_in_run > 0:
                X_test_tps = X_test[tp_mask]
                X_test_non_tps = X_test[~tp_mask]
                shap_values_tps = shap_values_run[tp_mask.values, :, j]
                for k, feature_name in enumerate(feature_cols):
                    feature_abund_in_tps = X_test_tps[feature_name]
                    high_abund_mask = (feature_abund_in_tps > 0)
                    num_abundant_tps = high_abund_mask.sum()
                    if num_abundant_tps > 0:
                        feature_shap_in_tps = shap_values_tps.values[:, k]
                        sum_shaps_in_high_abund = feature_shap_in_tps[high_abund_mask].sum()
                        if sum_shaps_in_high_abund > 0:
                            mean_abund_in_tps = feature_abund_in_tps.mean()
                            mean_abund_in_non_tps = X_test_non_tps[feature_name].mean()
                            fold_change = (mean_abund_in_tps + 1e-9) / (mean_abund_in_non_tps + 1e-9)
                            all_targeted_results_list.append({'Run': i, 'Feature': feature_name, 'Class': class_name,'Potency_SHAP': sum_shaps_in_high_abund / num_abundant_tps,'Overall_Impact_SHAP': sum_shaps_in_high_abund / num_tps_in_run,'Mean_Abund_in_High_Abund_TPs': feature_abund_in_tps[high_abund_mask].mean(),'Proportion_of_TPs_with_Feature': num_abundant_tps / num_tps_in_run,'Fold_Change_TP_vs_nonTP': fold_change,'Avg_Mean_Abund_in_High_Abund_TPs': feature_abund_in_tps[high_abund_mask].mean()})

    print(f"\n--- Aggregating and saving results for {output_prefix} ---")
    if all_performance_metrics:
        metrics_df = pd.DataFrame(all_performance_metrics)
        metrics_summary_df = pd.DataFrame({'Mean': metrics_df.mean(), 'Std_Dev': metrics_df.std()})
  
        output_path = os.path.join(results_output_dir, f"rf_performance_metrics_{output_prefix}.csv")
        metrics_summary_df.to_csv(output_path)
        print(f"Saved performance metrics to: {output_path}")

    if all_general_summary_rows:
        general_df_raw = pd.DataFrame(all_general_summary_rows)
        general_df_agg = general_df_raw.groupby(['Feature', 'Class']).agg(Mean_Abs_SHAP_Importance=('Mean_Abs_SHAP_per_run', 'mean'),SD_Abs_SHAP_Importance=('Mean_Abs_SHAP_per_run', 'std')).reset_index().sort_values(by=['Class', 'Mean_Abs_SHAP_Importance'], ascending=[True, False])
    
        output_path = os.path.join(results_output_dir, f"rf_general_shap_summary_{output_prefix}.csv")
        general_df_agg.to_csv(output_path, index=False)
        print(f"Saved: {output_path}")

    if all_targeted_results_list:
        targeted_df_raw = pd.DataFrame(all_targeted_results_list)
        targeted_df_agg = targeted_df_raw.groupby(['Feature', 'Class']).agg(Avg_Potency_SHAP=('Potency_SHAP', 'mean'),Avg_Overall_Impact_SHAP=('Overall_Impact_SHAP', 'mean'),Avg_Mean_Abund_in_High_Abund_TPs=('Mean_Abund_in_High_Abund_TPs', 'mean'),Avg_Proportion_of_TPs_with_Feature=('Proportion_of_TPs_with_Feature', 'mean'),Avg_Fold_Change_TP_vs_nonTP=('Fold_Change_TP_vs_nonTP', 'mean'),Num_Runs_Found=('Run', 'nunique')).reset_index().sort_values(by=['Class', 'Avg_Overall_Impact_SHAP'], ascending=[True, False])
  
        output_path = os.path.join(results_output_dir, f"rf_targeted_shap_summary_{output_prefix}.csv")
        targeted_df_agg.to_csv(output_path, index=False)
        print(f"Saved: {output_path}")

    if all_confusion_matrices:
     
        all_cm_np = np.array(all_confusion_matrices)
        avg_cm = np.mean(all_cm_np, axis=0)
        std_cm = np.std(all_cm_np, axis=0)
        avg_cm_df = pd.DataFrame(avg_cm, index=class_names, columns=class_names)
        avg_cm_df.index.name = 'True Label'

        csv_path = os.path.join(results_output_dir, f"rf_avg_confusion_matrix_{output_prefix}.csv")
        avg_cm_df.to_csv(csv_path)
        print(f"Saved average confusion matrix to: {csv_path}")


        annot_labels = np.array([f"{avg:.1f}\n(±{std:.1f})" for avg, std in zip(avg_cm.ravel(), std_cm.ravel())]).reshape(avg_cm.shape)
        plt.figure(figsize=(10, 8))
        sns.heatmap(avg_cm, annot=annot_labels, fmt='', cmap='Greens', xticklabels=class_names, yticklabels=class_names, cbar_kws={'label': 'Average Count'})
        plt.title(f'Average Confusion Matrix for {output_prefix} over {len(all_confusion_matrices)} Runs\n(Mean ± SD)', fontsize=16)
        plt.ylabel('True Label')
        plt.xlabel('Predicted Label')
   
        plot_path = os.path.join(results_output_dir, f"rf_avg_confusion_matrix_{output_prefix}.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved confusion matrix heatmap to: {plot_path}")

# ===================================================================
# --- SCRIPT EXECUTION ORCHESTRATOR 
# ===================================================================
if __name__ == '__main__':
    # This script is designed to be called with two command-line arguments.
    # It automatically finds and evaluates test sets, saving all outputs
    # to the specified directories.
    #
    # Usage:
    #   python  02_RandomForestTrainingAndShapleyAnalysisForLegacyData.py <input_dir> <output_dir>
    #
    
    # --- 1. ARGUMENT PARSING AND VALIDATION ---
    if len(sys.argv) != 3:
        print("Error: Invalid number of arguments.")
        print("Usage: RandomForestTrainingAndShapleyAnalysisForLegacyData.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_dir = sys.argv[2] # This directory is for BOTH models and results

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory not found at '{input_dir}'")
        sys.exit(1)
    
    # --- 2. SETUP ---
    N_RUNS = 20
    print(f"*** RUNNING PYTHON ANALYSIS ***")
    print(f"    Reading data from: {input_dir}")
    print(f"    Saving models and results to: {output_dir}")
    os.makedirs(output_dir, exist_ok=True)
    
    # --- 3. PHASE 1: Train and save all models ---
    train_file = os.path.join(input_dir, "train_data_fertilizer.csv")
    train_and_save_models(
        train_csv_path=train_file,
        model_output_dir=output_dir, # Models are saved in the main output dir
        n_runs=N_RUNS
    )
    
    # --- 4. Get class names from training data ---
    train_df = pd.read_csv(train_file, index_col=0)
    class_names = list(train_df['response_var'].astype('category').cat.categories)
    
    # --- 5. PHASE 2: Evaluate on the main test set (with a fixed name) ---
    print("\nEvaluating on main test set...")
    main_test_file = os.path.join(input_dir, "test_data_main_fertilizer.csv")
    if os.path.exists(main_test_file):
        evaluate_models_on_test_set(
            test_csv_path=main_test_file,
            model_dir=output_dir,
            results_output_dir=output_dir, # Results go to the same dir
            output_prefix="main",
            class_names=class_names,
            n_runs=N_RUNS
        )
    else:
        print(f"  - Warning: Main test file not found at '{main_test_file}', skipping.")

    # --- 6. PHASE 2 (cont.): Evaluate on future time-point test sets ---
    print("\nEvaluating on future time-point test sets...")
    
    for filename in os.listdir(input_dir):
        match = re.search(r'(\d+DAS)', filename)
        if match:
            test_file_path = os.path.join(input_dir, filename)
            time_point = match.group(1)
            output_prefix = f"test_{time_point}"
            
            print(f"  - Evaluating {filename} -> output prefix: '{output_prefix}'")
            evaluate_models_on_test_set(
                test_csv_path=test_file_path,
                model_dir=output_dir,
                results_output_dir=output_dir, # Results go to the same dir
                output_prefix=output_prefix,
                class_names=class_names,
                n_runs=N_RUNS
            )

    print("\n\nAll Python training and evaluation complete!")