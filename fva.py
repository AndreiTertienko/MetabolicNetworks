import os
import shutil
import subprocess
import time
import csv
import re
from pathlib import Path
import itertools
import json

#thresholds = ['0.85', '0.90', '0.95']
waiting_time = 0
thresholds = ['0.95']
output_basename = 'fva_results'
output_suffices = []
for threshold in thresholds:
    output_suffices.append(f'{output_basename}{str(threshold)}.csv')

def check_for_outputs(base_dir, sample_ids):
    """
    1) Goes into base_dir
    2) For each sampleID, enters the matching folder
    3) Iterates over all subfolders in that folder (bins)
    4) Checks for fva output files
    5) Counts them
    """
    fva_missing = 0
    bin_counter = 0 
    print('Counting the existing fva outputs')
    for cohort, sample_id_list in sample_ids.items():
        for sample_id in sample_id_list:
            sample_path = os.path.join(base_dir, "ecGEMs", cohort, sample_id)
            if not os.path.isdir(sample_path):
                print(f"Sample folder not found: {sample_path}")
                continue
            
            models = [ name for name in os.listdir(sample_path) if  name.endswith('.yml')]

            output_names = []
            for m in models:
                for s in output_suffices:
                    output_names.append(f'{m}_{s}') 
            fva_paths = [os.path.join(base_dir, 'fva_results', cohort, sample_id, n) for n in output_names]
            bin_counter += len(models)

            for path in fva_paths:
                if not os.path.isfile(path):
                    fva_missing += 1

    fva_to_do = bin_counter * len(thresholds)
    return fva_to_do, fva_missing

def run_fva_for_samples(base_dir, sample_ids, fva_source_path, waiting_time):
    """
    1) Goes into base_dir (typically scratch/ecGEMs).
    2) For each sampleID, enters the matching folder.
    3) Iterates over all subfolders in that folder (bins).
    4) Edits the corresponding Adapter file line (obj.params.path).
    5) Copies fva.m into the bin folder and updates the first two lines (bin_folder and adapter_file).
    6) Runs fva.m in MATLAB (no display / GUI).
    """
    taxon_dict = {}
    abund_dict = {}
    taxonomies = {}
    for cohort, sample_id_list in sample_ids.items():
        for sample_id in sample_id_list:
            sample_path = os.path.join(base_dir, "ecGEMs", cohort, sample_id)
            if not os.path.isdir(sample_path):
                print(f"Sample folder not found: {sample_path}")
                continue
            
            taxonomy_filename = f'{sample_id}_gtdbtk_summary.tsv'
            taxonomy_path = os.path.join(base_dir, "GTDBtk", sample_id, taxonomy_filename)

            taxonomies[sample_id] = taxonomy_path

            models = []  # dynamic list for any number of models
            for root, dirs, files in os.walk(sample_path):
                for file in files:
                    if file.endswith(".yml"):
                        model_name = Path(file).stem
                        models.append(model_name)

            if not os.path.isfile(fva_source_path):
                print(f"No fva.m file found: {fva_source_path}")
                continue
        
            target_fva = os.path.join(sample_path, "fva.m")
            shutil.copy(fva_source_path, target_fva)

            with open(target_fva, 'r') as f:
                fvalines = f.readlines()

            models_str = ", ".join([f'"{k}"' for k in models])
            models_line = f"models = [{models_str}];\n"
            fvalines[0] = models_line

            with open(target_fva, 'w') as f:
                f.writelines(fvalines)

            # STEP 6: Run fva.m in MATLAB if the output files do not exist yet
            output_names = []
            for m in models:
                for s in output_suffices:
                    output_names.append(f'{m}.yml_{s}') 
            fva_paths = [os.path.join(base_dir, 'fva_results', cohort, sample_id, n) for n in output_names]
            all_files_exist = True
            for path in fva_paths:
                if not os.path.isfile(path):
                    all_files_exist = False
                    break
            
            if not all_files_exist:
                fva_total, fva_missing = check_for_outputs(base_dir, sample_ids)
                fva_done = fva_total - fva_missing
                print(f'{fva_done} fva results already created out of {fva_total}')
                print(f'left to do: {fva_missing}')

                matlab_cmd = [
                    "matlab", "-nodisplay", "-nosplash", "-nodesktop",
                    "-r", f"cd('{sample_path}'); fva; exit"
                ]
                print("Running MATLAB command:", " ".join(matlab_cmd))
                subprocess.run(matlab_cmd)
                time.sleep(waiting_time)
            elif all_files_exist:
                print(f"all files exist for {sample_path}")

            taxon_dict, abund_dict = merge_fva_results(base_dir, cohort, sample_id, models, taxonomies, taxon_dict, abund_dict)
    return taxon_dict, abund_dict

def merge_fva_results(base_dir, cohort, sample_id, list_of_models, taxonomy_paths, taxon_dict, abund_dict):
    for i in range(len(list_of_models)):
        model = list_of_models[i]
        binname = model.split('Ready_ecGEM')[0]
        parts = re.split('(\d+)', binname)
        bin_name = f'{parts[0]}.{parts[1]}.{parts[2]}'
        taxonomy_path = taxonomy_paths[sample_id]
        if not os.path.exists(taxonomy_path):
            print(f'taxonomy file for {sample_id} not found at {taxonomy_path}')
            return None
        species = ''
        with open(taxonomy_path, 'r', encoding='utf-8') as f:
            for line in f:
                parts = line.split('\t')
                if bin_name == parts[0]:
                    taxon = parts[1]
                    species = taxon.split(';')[-1]
                    break
            if species == 's__':
                continue
            if not species:
                print(f'taxonomy for {sample_id}, {bin_name} not found')

        path_to_fvas = os.path.join(base_dir, "fva_results", cohort, sample_id)
        if os.path.isdir(path_to_fvas): 
            for file in os.listdir(path_to_fvas):
                if file.endswith(".csv") and model in file:
                    full_path = os.path.join(path_to_fvas, file)
                    if taxon not in taxon_dict:
                        taxon_dict[taxon] = {}
                        taxon_dict[taxon][cohort] = [full_path]
                    elif taxon in taxon_dict and not cohort in taxon_dict[taxon]:
                        taxon_dict[taxon][cohort] = [full_path]
                    elif taxon in taxon_dict and cohort in taxon_dict[taxon]:
                        taxon_dict[taxon][cohort].append(full_path)
                    else:
                        print('Something is wrong with the taxon_dict logic')
        else:
            print(f'directory {path_to_fvas} does not exist')

        if taxon not in abund_dict:
            abund_dict[taxon] = {}
            abund_dict[taxon][cohort] = [list((sample_id, bin_name))]
        elif taxon in abund_dict and not cohort in abund_dict[taxon]:
            abund_dict[taxon][cohort] = [list((sample_id, bin_name))]
        elif taxon in abund_dict and cohort in abund_dict[taxon]:
            abund_dict[taxon][cohort].append(list((sample_id, bin_name)))
        else:
            print('Something is wrong with the abund_dict logic')

    return taxon_dict, abund_dict

def save_merged_results(base_dir, taxon_dict, thresholds):
    output_dir = os.path.join(base_dir, "merged_fva")
    os.makedirs(output_dir, exist_ok=True)
    out_paths = []
    zeros_options = ("_zeros", "_nozeros", "_unfiltered")
    # Now merge files for each (taxon, filename) pair
    for threshold in thresholds:
        for zeros in zeros_options:
            for taxon, cohort_dict in taxon_dict.items():
                data_rows_all = []
                for cohort, csv_paths in cohort_dict.items():
                    data_rows_cohort = []
                    for csv_path in csv_paths:
                        fname = os.path.basename(csv_path)

                        if not csv_path:
                            print("csv does not exist")
                            continue
                        
                        if zeros != "_unfiltered":
                            if not zeros in fname:
                                #print(f"no {zeros} in {fname}")
                                continue
                        else:
                            if zeros_options[0] in fname or zeros_options[1] in fname:
                                #it is a filtered result that we want to skip because we not search for unfiltered
                                continue

                        if not threshold in fname:
                            #print(f"no {threshold} in {fname}")
                            continue

                        header = None
                        with open(csv_path, 'r', newline='', encoding='utf-8') as f:
                            reader = csv.reader(f)
                            file_header = next(reader, None)  # First row
                            if file_header and header is None:
                                header = file_header
                            elif file_header and file_header != header:
                                print(f"Warning: CSV headers differ in '{csv_path}'. Merging may be misaligned.")
                            for row in reader:
                                data_rows_cohort.append(row)
                                data_rows_all.append(row)

                    # Build output name as "<taxon>_<original filename>"
                    if not data_rows_cohort:
                        #print(f"nothing found for {cohort}, {taxon}, keep zeros = {zeros}, {threshold}")
                        continue
                    else:
                        data_rows_cohort.sort()
                        list(data_rows_cohort for data_rows_cohort,_ in itertools.groupby(data_rows_cohort))
                    out_name = f"{taxon}_{threshold}{zeros}.csv"
                    out_path = os.path.join(output_dir, cohort, out_name)
                    with open(out_path, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f)
                        if header:
                            writer.writerow(header)
                        writer.writerows(data_rows_cohort)

                        #print(f"Merged file for TaxID {taxon} -> {out_path}")
                        out_paths.append(out_path)

                if not data_rows_all:
                    #print(f"nothing found for {cohort}, {taxon}, keep zeros = {zeros}, {threshold}")
                    continue
                else:
                    data_rows_all.sort()
                    list(data_rows_all for data_rows_all,_ in itertools.groupby(data_rows_all))
                    out_name = f"{taxon}_{threshold}{zeros}.csv"
                    out_path = os.path.join(output_dir, 'merged_cohorts', out_name)
                    with open(out_path, 'w', newline='', encoding='utf-8') as f:
                        writer = csv.writer(f)
                        if header:
                            writer.writerow(header)
                        writer.writerows(data_rows_all)
                        out_paths.append(out_path)
   
    return output_dir, out_paths

def create_networks(fva_dir, cohorts, input_csvs, keep_zeros, thresholds):
    """
    Merges all *_filtered_zeros.csv or *_filtered_nozeros.csv files (depending on keep_zeros) 
    by threshold (the two digits in *.0.??_filtered*.csv). Only columns 0 (Reaction) and 3 (Direction)
    are included in the final merged file. The original files are left intact.
    """
    # Find the files to include in the network
    for cohort in cohorts:
        if keep_zeros == "zeros":
            merged_fvas = [f for f in input_csvs if (f.endswith("_zeros.csv") and cohort in f)]
        elif keep_zeros == "nozeros":
            merged_fvas = [f for f in input_csvs if (f.endswith("_nozeros.csv") and cohort in f)]
        elif keep_zeros == "unfiltered":
            merged_fvas = [f for f in input_csvs if (f.endswith("_unfiltered.csv") and cohort in f)]
        elif keep_zeros == "complete":
            merged_fvas = [f for f in input_csvs if (f.endswith("_unfiltered.csv") and cohort in f)]

        
        # For each threshold, merge all files containing that threshold
        for th in thresholds:
            all_rows = []
            for file_path in merged_fvas:
                file_name = os.path.basename(file_path)
                if th in file_name:
                    taxon = file_name.split('_0.')[0]
                    with open(file_path, 'r', newline='', encoding='utf-8') as f:
                        reader = csv.reader(f)
                        header = next(reader, None)
                        if not header:
                            continue
                        # Gather only columns 0 and 3 (Reaction, Direction)
                        for row in reader:
                            if len(row) < 4:
                                continue
                            if float(row[1]) >= 0 and float(row[2]) >= 0:
                                positive = True
                            elif float(row[1]) <= 0 and float(row[2]) <= 0:
                                positive = False
                            if keep_zeros == "unfiltered":
                                if not (float(row[1]) > 0 and float(row[2]) > 0) and not (float(row[1]) < 0 and float(row[2])) < 0:
                                    all_rows.append([taxon, row[0], "unspecific"])
                                elif row[3] == "absorbed" and positive:
                                    all_rows.append([taxon, row[0], "secreted"])
                                elif row[3] == "absorbed" and not positive:
                                    all_rows.append([taxon, row[0], "uptaken"])
                                elif row[3] == "excreted" and positive:
                                    all_rows.append([taxon, row[0], "uptaken"])
                                elif row[3] == "excreted" and not positive:
                                    all_rows.append([taxon, row[0], "secreted"])
                                else:
                                    print(f"confusion with secrete/uptake, skipping {file_path}") 
                                positive = None
                            elif keep_zeros == "complete":
                                if not (float(row[1]) > 0 and float(row[2]) > 0) and not (float(row[1]) < 0 and float(row[2])) < 0:
                                    all_rows.append([taxon, row[0], row[1], row[2], "unspecific"])
                                elif row[3] == "absorbed" and positive:
                                    all_rows.append([taxon, row[0], row[1], row[2], "secreted"])
                                elif row[3] == "absorbed" and not positive:
                                    all_rows.append([taxon, row[0], row[1], row[2], "uptaken"])
                                elif row[3] == "excreted" and positive:
                                    all_rows.append([taxon, row[0], row[1], row[2], "uptaken"])
                                elif row[3] == "excreted" and not positive:
                                    all_rows.append([taxon, row[0], row[1], row[2], "secreted"])
                                else:
                                    print(f"confusion with secrete/uptake, skipping {file_path}") 
                                positive = None
                            else:
                                if row[3] == "excreted" and positive:
                                    all_rows.append([row[0], taxon])
                                elif row[3] == "excreted" and not positive:
                                    all_rows.append([taxon, row[0]])
                                elif row[3] == "absorbed" and positive:
                                    all_rows.append([taxon, row[0]])
                                elif row[3] == "absorbed" and not positive:
                                    all_rows.append([row[0], taxon])
                                else:
                                    print(f"confusion with secrete/uptake, skipping {file_path}")
                                positive = None

            all_rows = [list(x) for x in dict.fromkeys(map(tuple, all_rows))]

            # Write out the merged result if we have any data
            if all_rows:
                out_name = f"network_{th}_{keep_zeros}.csv"
                out_path = os.path.join(fva_dir, cohort, out_name)
                with open(out_path, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    # The columns are already named Reaction and Direction in the CSV data,
                    # but we need a header row in the merged file:
                    writer.writerows(all_rows)
                print(f"Merged {keep_zeros} files for threshold {th} -> {out_path}")

def diff_abundance(base_dir, abund_dict):
    from scipy.stats import ranksums, ttest_ind, mannwhitneyu
    from statsmodels.stats.multitest import multipletests
    p_values = {}
    for taxon, cohorts_dict in abund_dict.items():
        abund_data = {}
        for cohort, list_of_tuples in cohorts_dict.items():
            for i in range(len(list_of_tuples)):
                sample_id = list_of_tuples[i][0]
                bin_name = list_of_tuples[i][1]
                abund_path = os.path.join(base_dir, 'abundance', f'{sample_id}.abund')
                with open(abund_path, 'r', encoding='utf-8') as f:
                    for line in f:
                        parts = line.split('\t')
                        if f'{bin_name}' in parts[0]:
                            if cohort in abund_data:
                                abund_data[cohort].append(parts[2].strip('\n'))
                            else:
                                abund_data[cohort] = [parts[2].strip('\n')]
                            break
                list_of_tuples[i].append(parts[2].strip('\n'))                          

        # If there are exactly two cohorts, perform the unpaired Wilcoxon (rank-sum) test
        if len(abund_data) == 2:
            cohs = list(abund_data.keys())
            group1 = [float(x) for x in abund_data[cohs[0]]]
            group2 = [float(x) for x in abund_data[cohs[1]]]
            # Mann-Whitney U
            stat, mwu_p = mannwhitneyu(group1, group2, alternative='two-sided')
            # Wilcoxon rank-sum test
            stat, p = ranksums(group1, group2)
            if p <= 0.05:
                print(f'{taxon} (Wilcoxon): {p}')
            # Student's t-test with unequal variances
            t_stat, t_p = ttest_ind(group1, group2, equal_var=False)
            if taxon == 's__CAG-217':
                print(mwu_p)
                print(p)
                print(t_p)
            if t_p <= 0.05:
                print(f'{taxon} (t-test): {t_p}')
            # Store both p-values in a dictionary for the current taxon
            p_values[taxon] = {'wilcox': p, 'ttest': t_p}
    
    # Perform BH correction on the collected p-values (separately for Wilcoxon and t-test)
    if p_values:
        taxon_list = list(p_values.keys())
        wilcox_vals = [p_values[t]['wilcox'] for t in taxon_list]
        ttest_vals  = [p_values[t]['ttest'] for t in taxon_list]
        _, wilcox_corr, _, _ = multipletests(wilcox_vals, alpha=0.05, method='fdr_bh')
        _, ttest_corr, _, _  = multipletests(ttest_vals, alpha=0.05, method='fdr_bh')
        for idx, taxon in enumerate(taxon_list):
            p_values[taxon]['wilcox_corrected'] = wilcox_corr[idx]
            p_values[taxon]['ttest_corrected'] = ttest_corr[idx]
            # Optionally print corrected p-values when significant
            if wilcox_corr[idx] <= 0.05:
                print(f'{taxon} (Wilcoxon BH-corrected): {wilcox_corr[idx]}')
            if ttest_corr[idx] <= 0.05:
                print(f'{taxon} (t-test BH-corrected): {ttest_corr[idx]}')
    
    return abund_dict, p_values

def plot_abundances(abund_dict, p_values, use_corrected=True):
    """
    Generates several plots based on abundance data.
      1. Heatmap with dendrogram: For each taxon, shows significance as -log10(p) 
         (derived from Wilcoxon test) and uses hierarchical clustering so a dendrogram 
         appears on the y-axis.
      2. Boxplot for significant taxa only.
      3. Density (KDE) plot overlaid for both cohorts.
    
    Inputs:
      abund_dict: Dictionary with structure {taxon: {cohort: [[sample_id, bin_name, abundance], ...], ...}, ...}
      p_values: Dictionary with keys = taxon and values including 'wilcox' and 'ttest' p-values
                plus (optionally) 'wilcox_corrected' and 'ttest_corrected'.
      use_corrected: Boolean flag to select corrected p-values if True, otherwise uncorrected.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    import numpy as np

    sns.set_style("whitegrid")
    sns.set_palette("colorblind")

    # Prepare DataFrame for boxplot and density plots: one row per measurement.
    rows = []
    for taxon, cohort_dict in abund_dict.items():
        for cohort, records in cohort_dict.items():
            for rec in records:
                try:
                    abundance = float(rec[2])
                except (ValueError, IndexError):
                    continue
                rows.append({'Taxon': taxon, 'Cohort': cohort, 'Abundance': abundance})
    df = pd.DataFrame(rows)
    if df.empty:
        print("No abundance data available for plotting.")
        return

    # Determine which p-values to use.
    get_val = (lambda taxon, key: p_values[taxon].get(key + '_corrected')
               if use_corrected and key + '_corrected' in p_values[taxon]
               else p_values[taxon].get(key, 1))

    # 1. Heatmap with dendrogram: Hierarchical clustering the taxa.
    p_rows = []
    for taxon, vals in p_values.items():
        p_val = get_val(taxon, 'wilcox')
        p_rows.append({'Taxon': taxon, 'p_value': p_val})
    df_p = pd.DataFrame(p_rows).set_index('Taxon')
    # Replace zeros with a small value to avoid -inf in transformation.
    df_p['-log10(p)'] = -np.log10(df_p['p_value'].replace(0, 1e-10))
    df_p = df_p.drop('p_value', axis = 1)
    # Use clustermap to display hierarchical clustering (dendrogram on y-axis).
    '''cg = sns.clustermap(df_p[['-log10(p)']],
                        cmap='viridis',
                        annot=True,
                        fmt='.2f',
                        cbar_kws={'label': '-log10(p)'},
                        dendrogram_ratio=(.2, .2),
                        metric='euclidean',
                        method='average',
                        cluster_cols=False)'''
    cg = sns.clustermap(df_p,
                        annot=True,
                        metric='euclidean',
                        method='average')
    plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    plt.title("Hierarchical Clustering Heatmap (Wilcoxon)" + (" [Corrected]" if use_corrected else ""))
    plt.show()

    # 2. Boxplot for significant taxa only (Wilcoxon p<=0.05)
    sig_taxa = [taxon for taxon, vals in p_values.items() if get_val(taxon, 'wilcox') <= 0.05]
    df_sig = df[df['Taxon'].isin(sig_taxa)]
    if not df_sig.empty:
        plt.figure(figsize=(12, 6))
        sns.boxplot(x='Taxon', y='Abundance', hue='Cohort', data=df_sig)
        plt.title("Microbe Abundance Boxplot (Significant Taxa: Wilcoxon p<=0.05)" +
                  (" [Corrected]" if use_corrected else ""))
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.show()
    else:
        print("No significant taxa found for the boxplot.")

    # 3. Density (KDE) plot: Overlaid for both cohorts.
    plt.figure(figsize=(8, 6))
    for cohort in df['Cohort'].unique():
        subset = df[df['Cohort'] == cohort]
        sns.kdeplot(subset['Abundance'], label=cohort, shade=True)
    plt.title("Density Plot of Microbe Abundances by Cohort")
    plt.xlabel("Abundance")
    plt.ylabel("Density")
    plt.legend()
    plt.tight_layout()
    plt.show()

    '''# 4. Scatter plot: x-axis Student's t-test p, y-axis Wilcoxon p for each taxon.
    scatter_rows = []
    for taxon in p_values.keys():
        t_p = get_val(taxon, 'ttest')
        w_p = get_val(taxon, 'wilcox')
        scatter_rows.append({'Taxon': taxon, 'ttest': t_p, 'wilcox': w_p})
    print(scatter_rows)
    df_scatter = pd.DataFrame(scatter_rows)
    plt.figure(figsize=(6, 6))
    plt.scatter(df_scatter['ttest'], df_scatter['wilcox'], color='blue', alpha=0.7)
    # Add cutoff lines at p = 0.05.
    plt.axvline(x=0.05, color='red', linestyle='--')
    plt.axhline(y=0.05, color='red', linestyle='--')
    plt.xlabel("Student's t-test p-value" + (" [Corrected]" if use_corrected else ""))
    plt.ylabel("Wilcoxon p-value" + (" [Corrected]" if use_corrected else ""))
    plt.title("Scatter Plot of p-values Across Taxa")
    plt.tight_layout()
    plt.show()'''
    
    # 5. Histogram: Distribution of abundance value counts per taxon.
    # First, prepare data: for each taxon and for each cohort, count the number of abundance values,
    # plus prepare a combined count "All" for each taxon.
    count_data = []
    for taxon, cohort_dict in abund_dict.items():
        total = 0
        for cohort, records in cohort_dict.items():
            count = len(records)
            count_data.append({'Taxon': taxon, 'Cohort': cohort, 'Count': count})
            total += count
        count_data.append({'Taxon': taxon, 'Cohort': 'All', 'Count': total})
    
    df_count = pd.DataFrame(count_data)
    # For each Cohort group, compute the frequency: number of taxa having each abundance count.
    freq_df = df_count.groupby("Cohort")["Count"].value_counts().rename("NumTaxa").reset_index()
    
    # Create a continuous range for x-axis, from 0 to the global maximum count.
    global_max = int(df_count["Count"].max())
    cohorts_list = df_count["Cohort"].unique()
    continuous_dfs = []
    for cohort in cohorts_list:
        # Get frequency data for the cohort and reindex to a continuous range.
        sub = freq_df[freq_df["Cohort"] == cohort].set_index("Count")[["NumTaxa"]]
        full_range = pd.DataFrame({"Count": np.arange(0, global_max + 1)})
        full_range = full_range.merge(sub, how="left", left_on="Count", right_index=True)
        full_range["NumTaxa"] = full_range["NumTaxa"].fillna(0)
        full_range["Cohort"] = cohort
        continuous_dfs.append(full_range)
    df_cont = pd.concat(continuous_dfs)
    
    plt.figure(figsize=(8, 6))
    sns.barplot(data=df_cont, x="Count", y="NumTaxa", hue="Cohort", palette="colorblind")
    plt.xscale("log")
    plt.title("Distribution of Taxa Abundance Counts by Cohort Group (Histogram, Log Scale)")
    plt.xlabel("Number of Abundance Values (log scale)")
    plt.ylabel("Number of Taxa")
    plt.tight_layout()
    plt.show()

def create_venn_metabolites(network_files, output_file):
    """
    Creates a Venn diagram showing how many metabolites are in each network file.
    
    Parameters:
       network_files: list of paths to network CSV files. The CSV files have no header.
                      The first two columns are checked for metabolite names starting with "EX_".
       circle_labels: list of labels (strings) to label each circle manually.
       output_file: full file path to which the diagram should be saved.
       
    The function extracts metabolite names by:
       - Reading the first two columns of each file.
       - Keeping names that start with "EX_"
       - Removing the "EX_" prefix.
    Only 2 or 3 network files are supported.
    """
    import csv
    import os
    import matplotlib.pyplot as plt

    # Import the appropriate venn diagram function based on number of sets
    if len(network_files) == 2:
        from matplotlib_venn import venn2
    elif len(network_files) == 3:
        from matplotlib_venn import venn3
    else:
        raise ValueError("Only 2 or 3 network files are supported for this Venn diagram.")
    
    # For each network file, create a set of metabolites
    metabolite_sets = []
    labels = []
    for fpath in network_files:
        mets = set()
        if not os.path.exists(fpath):
            print(f"File not found: {fpath}")
            metabolite_sets.append(mets)
            continue
        with open(fpath, "r", newline='', encoding="utf-8") as f:
            reader = csv.reader(f)
            for row in reader:
                # Check first two columns for names starting with 'EX_'
                for item in row[:2]:
                    if item.startswith("EX_"):
                        mets.add(item.replace("EX_", "", 1))
        metabolite_sets.append(mets)
        label = f'{os.path.basename(fpath).split("_")[1]}: {len(mets)}'
        labels.append(label)
        print(label)
    
    # Plot Venn diagram based on the number of sets
    plt.figure(figsize=(8,6))
    if len(metabolite_sets) == 2:
        v = venn2(subsets = metabolite_sets)
    elif len(metabolite_sets) == 3:
        v = venn3(subsets = metabolite_sets)
        v.get_label_by_id('001').set_text(labels[2])
        v.get_label_by_id('001').set_fontsize(8)
        v.get_label_by_id('011').set_text(labels[1])
        v.get_label_by_id('011').set_fontsize(8)
        v.get_label_by_id('111').set_text(labels[0])
        v.get_label_by_id('111').set_fontsize(8)
        for text in v.set_labels:
            text.set_fontsize(1)

    plt.title("Metabolite Distribution Across Networks")
    # Save the diagram
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.show()
    plt.savefig(output_file, dpi=300)
    plt.close()
    print(f"Venn diagram saved to: {output_file}")

def save_dictionaries(abund_dict, taxonomy_dict, p_values, output_dir=r"C:\Users\A\Documents\MSc_Thesis\output"):
    """
    Saves the abundance dictionary, taxonomy dictionary, and p-values dictionary to JSON files.
    
    Parameters:
      abund_dict: The abundance dictionary.
      taxonomy_dict: The taxonomy dictionary.
      p_values: Dictionary containing p-values for each taxon.
      output_dir: Directory where the JSON files will be saved.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    abund_path = os.path.join(output_dir, "abund_dict.json")
    taxonomy_path = os.path.join(output_dir, "taxonomy_dict.json")
    pvalues_path = os.path.join(output_dir, "p_values.json")
    
    with open(abund_path, "w") as f:
        json.dump(abund_dict, f, indent=4)
    with open(taxonomy_path, "w") as f:
        json.dump(taxonomy_dict, f, indent=4)
    with open(pvalues_path, "w") as f:
        json.dump(p_values, f, indent=4)
    
    print(f"Abundance dictionary saved to: {abund_path}")
    print(f"Taxonomy dictionary saved to: {taxonomy_path}")
    print(f"P-values dictionary saved to: {pvalues_path}")

def create_absolute_abundance_matrix(base_dir, taxon_dict, output_dir=None):
    """
    Creates a matrix of absolute abundance values from reads.tsv.
    Uses the taxon_dict to map bin names to taxa.
    
    Parameters:
        base_dir: Base directory where reads.tsv is found
        taxon_dict: Dictionary mapping taxa to cohorts and file paths
        output_dir: Directory to save output file (defaults to base_dir/output)
        
    Returns:
        DataFrame with taxa as rows and samples as columns
    """
    import pandas as pd
    import re
    
    if output_dir is None:
        output_dir = os.path.join(base_dir, "output")
    os.makedirs(output_dir, exist_ok=True)
    
    # Read the reads.tsv file
    reads_file = os.path.join(base_dir, "abundance", "reads.tsv")
    if not os.path.exists(reads_file):
        print(f"Error: File not found: {reads_file}")
        return None
        
    # Read the TSV file (assuming tab-separated with no header)
    reads_df = pd.read_csv(reads_file, sep='\t', header=None, 
                          names=['sample_id', 'bin_name', 'read_count'])
    
    print(f"Read {len(reads_df)} records from reads.tsv")
    
    # Create a mapping from (sample_id, bin_name) to taxon
    # The taxonomy dict contains file paths with sample_id and bin name embedded
    bin_to_taxon = {}
    
    for taxon, cohort_dict in taxon_dict.items():
        for cohort, file_paths in cohort_dict.items():
            for path in file_paths:
                # Extract sample_id from file path
                match = re.search(r'Feces_[^\\]+', path)
                if match:
                    sample_id = match.group(0)
                    
                    # Extract bin name from file path
                    bin_match = path.split('\\')[-1].split('Ready')[0]
                    bin_name = '.'.join(re.split('(\d+)', bin_match))
                    # Create key for the bin_to_taxon dictionary
                    bin_to_taxon[(sample_id, bin_name)] = taxon

    # Debug - how many bins were mapped
    print(f"Mapped {len(bin_to_taxon)} sample-bin combinations to taxa")
    
    # Process each reads.tsv record and map to taxa
    abundance_matrix = {}
    matched_count = 0
    
    for _, row in reads_df.iterrows():
        sample_id = row['sample_id']
        bin_name = row['bin_name']
        read_count = int(row['read_count'])

        bin_prefix = bin_name.split('.')[-1]
        
        if bin_prefix:
            # Look up the taxon for this sample_id and bin
            taxon = bin_to_taxon.get((sample_id, bin_name))
            
            if taxon:
                matched_count += 1
                if taxon not in abundance_matrix:
                    abundance_matrix[taxon] = {}
                abundance_matrix[taxon][sample_id] = abundance_matrix[taxon].get(sample_id, 0) + read_count
    
    print(f"Matched {matched_count} records from reads.tsv to taxa")
    
    # Convert to DataFrame - this is the matrix
    all_taxa = sorted(abundance_matrix.keys())
    all_samples = sorted(set(sum([list(abundance_matrix[t].keys()) for t in all_taxa], [])))
    
    # Initialize dataframe with zeros
    matrix_df = pd.DataFrame(0, index=all_taxa, columns=all_samples)
    
    # Fill in the counts
    for taxon, samples in abundance_matrix.items():
        for sample, count in samples.items():
            matrix_df.at[taxon, sample] = count
    
    # Save to file for R
    output_path = os.path.join(output_dir, "absolute_abundance_matrix.tsv")
    matrix_df.to_csv(output_path, sep='\t')
    print(f"Saved absolute abundance matrix to: {output_path}")
    
    return matrix_df


if __name__ == "__main__":
    base_directory = r"C:\Users\A\Documents\MSc_Thesis"
    sampleIDs = {"day0":["Feces_10010_Visit1_S106", "Feces_20004_Visit1_S5", "Feces_20030_Visit1_S8", \
                        "Feces_20054_Visit1_S86", "Feces_20073_Visit1_S51", "Feces_20100_Visit1_S67", \
                        "Feces_10022_Visit1_S62", "Feces_20006_Visit1_S71", "Feces_20035_Visit1_S76", \
                        "Feces_20059_Visit1_S50", "Feces_20074_Visit1_S50", "Feces_20110_Visit1_S55", \
                        "Feces_10028_Visit1_S98", "Feces_20008_Visit1_S29", "Feces_20046_Visit1_S12", \
                        "Feces_20070_Visit1_S48", "Feces_20085_Visit1_S142", "Feces_20116_Visit1_S27", \
                        "Feces_20002_Visit1_S27", "Feces_20010_Visit1_S24", "Feces_20049_Visit1_S29", \
                        "Feces_20071_Visit1_S41", "Feces_20097_Visit1_S15"], \
                "day84":["Feces_10010_Visit3_S103", "Feces_20004_Visit3_S9", "Feces_20030_Visit3_S85", \
                        "Feces_20054_Visit3_S14", "Feces_20073_Visit3_S2", "Feces_20100_Visit3_S9", \
                        "Feces_10022_Visit3_S111", "Feces_20006_Visit3_S7", "Feces_20035_Visit3_S45", \
                        "Feces_20059_Visit3_S36", "Feces_20074_Visit3_S25", "Feces_20110_Visit3_S15", \
                        "Feces_10028_Visit3_S100", "Feces_20008_Visit3_S14", "Feces_20046_Visit3_S82", \
                        "Feces_20070_Visit3_S49", "Feces_20085_Visit3_S29", "Feces_20116_Visit3_S25", \
                        "Feces_20002_Visit3_S87", "Feces_20010_Visit3_S6","Feces_20049_Visit3_S37", \
                        "Feces_20071_Visit3_S121", "Feces_20097_Visit3_S10"]}
    #sampleIDs = {"day84": ["Feces_10010_Visit3_S103"]}
    fva_script_path = r"C:\Users\A\Documents\MSc_Thesis\codes\fva.m"
    cohorts = list(sampleIDs.keys())
    cohorts.append("merged_cohorts")
    keep_zeros = ("complete", "zeros", "nozeros", "unfiltered")

    taxon_dict, abund_dict = run_fva_for_samples(base_directory, sampleIDs, fva_script_path, waiting_time)
    merged_dir, merged_files = save_merged_results(base_directory, taxon_dict, thresholds)
    
    for zeros_option in keep_zeros:
        create_networks(merged_dir, cohorts, merged_files, zeros_option, thresholds)
    
    #abund_dict, p_values = diff_abundance(base_directory, abund_dict)
    #plot_abundances(abund_dict, p_values)

    network_files = [
        r"C:\Users\A\Documents\MSc_Thesis\merged_fva\merged_cohorts\network_0.85_nozeros.csv",
        r"C:\Users\A\Documents\MSc_Thesis\merged_fva\merged_cohorts\network_0.90_nozeros.csv", 
        r"C:\Users\A\Documents\MSc_Thesis\merged_fva\merged_cohorts\network_0.95_nozeros.csv"
    ]
    output_file = r"C:\Users\A\Documents\MSc_Thesis\output\nozeros_venn.png"

    #create_venn_metabolites(network_files, output_file)

    #save_dictionaries(abund_dict, taxon_dict, p_values)
    #create_absolute_abundance_matrix(base_directory, taxon_dict, output_dir=None)
