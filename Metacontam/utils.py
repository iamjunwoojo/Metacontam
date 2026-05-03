import collections
import pandas as pd
import numpy as np
import os
import glob
import subprocess


def parse_metadata(metadata_path):
    metadata_parsed=[]
    with open(metadata_path, 'r') as file:
        for line in file:
            sample_name=line.split("\t")[0].strip()
            sample_type=line.split("\t")[1].strip()
            sample1_path=line.split("\t")[2].strip()
            sample2_path=line.split("\t")[3].strip()
            metadata_parsed.append((sample_name,sample_type,sample1_path,sample2_path))
    return metadata_parsed

def parse_metadata_for_preval(META_DATA):
    # parsing metadata into multiple dict()
    # store metadata for each sample
    meta_data = collections.defaultdict(list)
    # store metadata for each sample type
    meta_type = collections.defaultdict(list)
    with open(META_DATA, "r") as meta:
        for line in meta.readlines():
            sample_id = line.strip().split("\t")[0]
            sample_type = line.strip().split("\t")[1]
            sample_read_1_abs_dir = line.strip().split("\t")[2]
            sample_read_2_abs_dir = line.strip().split("\t")[3]
            meta_data[sample_id] = [sample_type, sample_read_1_abs_dir, sample_read_2_abs_dir]
            meta_type[sample_type].append(sample_id)

    return meta_data, meta_type




def process_sample_report(input_report_path):
    """
    Reads the given Kraken2 report, extracts information for species class 'S',
    and collects data for samples with taxid '0' or '1', returning the results.

    Parameters:
    input_report_path (str): Path to the Kraken2 report file

    Returns:
    pd.DataFrame: DataFrame with counts for each taxid related to the sample
    dict: Dictionary with the total read count for each sample
    """
    # Open the file
    with open(input_report_path) as file:
        total_read_dict = collections.defaultdict(int)  # Store total reads per sample
        dict_count = {}  # Store counts for taxid when species class is 'S'

        # Extract sample name from the report file path
        sample_name = input_report_path.split('/')[-1].split('.')[0]

        # Process each line in the report
        for line in file:
            fields = line.split("\t")
            count = int(fields[1])  # Make sure count is an integer
            species_class = fields[3]
            taxid = fields[4]
            sciname = fields[5]

            # Add the count to total reads if taxid is '0' or '1'
            if taxid in ('0', '1'):
                total_read_dict[sample_name] += count

            # Store the count in dict_count if species_class is 'S'
            if species_class == 'S':
                dict_count[taxid] = [count]

    # Return the results as a DataFrame and total read dictionary
    return pd.DataFrame(dict_count, index=[sample_name]), total_read_dict


def combine_reports(input_report_paths):
    """
    Processes multiple Kraken2 reports and combines the information into a single DataFrame.

    Parameters:
    input_report_paths (list): List of Kraken2 report file paths

    Returns:
    pd.DataFrame: Combined DataFrame with all the taxid-related information
    """
    # Initialize a list to hold DataFrames for each report
    df_list = []
    for report_path in input_report_paths:
        df, _ = process_sample_report(report_path)
        df_list.append(df)

    # Concatenate the DataFrames and fill NaN values with 0
    combined_df = pd.concat(df_list).fillna(0)

    # Convert all columns to int to avoid float (.0) values
    combined_df = combined_df.astype(int)

    return combined_df











def get_prevalence_type(meta_data, meta_type, taxid2prevalence, sample_types,
                        output_dir ,min_reads, min_abundance=0, bracken_dir=None):
    prevalence_type = dict()
    for taxid in taxid2prevalence:
        prevalence_type[taxid] = np.zeros(len(sample_types), dtype=float)

    taxon_prevalence_dict = collections.defaultdict(float)
    taxon_name_dict = dict()
    if  bracken_dir:
        Bracken_dir=bracken_dir
    elif not  bracken_dir:
        Bracken_dir=os.path.join(output_dir, "Bracken_dir")

    for sample_id in meta_data:
        sample_dir = os.path.join(Bracken_dir, f"{sample_id}_bracken_species.report")
        sample_type = meta_data[sample_id][0]
        #print(sample_types.index(sample_type))
        with open(sample_dir, "r") as report:
            lines = report.readlines()
            total_reads = 0
            for line in lines:
                taxon_id = line.split("\t")[4]
                read_count = int(line.split("\t")[1])
                taxon_name = line.split("\t")[-1].strip()
                taxon_rank = line.split("\t")[3]

                # get total number of reads by adding classified and unclassified reads
                if taxon_id == "0":
                    total_reads += read_count
                elif taxon_id == "1":
                    total_reads += read_count
                elif (read_count >= min_reads or read_count/total_reads >= min_abundance) \
                    and taxon_id in taxid2prevalence:
                    taxon_name_dict[taxon_id] = taxon_name
                    #print(taxon_id, prevalence_type[taxon_id][sample_types.index(sample_type)])
                    prevalence_type[taxon_id][sample_types.index(sample_type)] += 1

    for taxon_id in prevalence_type:
        for i in range(len(sample_types)):
            prevalence_type[taxon_id][i] /= len(meta_type[sample_types[i]])

    return prevalence_type



def filtering_matrix(input_dataframe,  min_reads, preval_threshold):
    input_dataframe = input_dataframe.transpose()
    list_bolean = []
    for i in range(input_dataframe.shape[0]):
        list_input = input_dataframe.iloc[i,:].values
        if len([i for i in list_input if float(i) >= min_reads])/len(list_input) > float(preval_threshold):
            list_bolean.append(True)
        else:
            list_bolean.append(False)
    return input_dataframe[list_bolean]






def filtering_count_matrix(input_dataframe, output, min_reads, preval_threshold):
    # 1. Remove taxid 9606 (human) column
    col_str = input_dataframe.columns.astype(str)
    if "9606" in col_str:
        input_dataframe = input_dataframe.loc[:, col_str != "9606"]
        print("9606 column removed.")
    else:
        print("9606 column not found.")

    # 2. Filter by prevalence threshold
    input_dataframe = input_dataframe.transpose()
    list_bolean = []
    for i in range(input_dataframe.shape[0]):
        list_input = input_dataframe.iloc[i, :].values
        if len([i for i in list_input if float(i) >= min_reads]) / len(list_input) > float(preval_threshold):
            list_bolean.append(True)
        else:
            list_bolean.append(False)
    input_dataframe[list_bolean].to_csv(os.path.join(output, "kraken_filtered_matrix.txt"), sep="\t")
















def Make_adjacent_matrix(Rscript_path, tmp_matrix_path, edge_path, coefficient_threshold=0):
    threshold = str(coefficient_threshold)
    output_dir = os.path.dirname(edge_path)
    os.makedirs(output_dir, exist_ok=True)
    subprocess.run(
        ["Rscript", Rscript_path, tmp_matrix_path, edge_path, threshold],
        check=True
    )





def find_best_threshold(input_df, min_reads, blacklist, Rscript_path,
                        top_n=100, shift_constant=0.7, output=".", 
                        tmp_matrix_name="tmp_netcomi_input.tsv", 
                        tmp_edge_name="tmp_netcomi_edges.tsv",
                        debug=False):

    # 1. Remove taxid 9606 (human) column
    col_str = input_df.columns.astype(str)
    print("columns[:10]:", col_str[:10])
    if "9606" in col_str:
        filtered_input_df = input_df.loc[:, col_str != "9606"]
        print("9606 column removed.")
    else:
        filtered_input_df = input_df.copy()
        print("9606 column not found.")

    # 2. Preprocess: transpose so rows=species, cols=samples
    df = filtered_input_df.transpose()
    n_samples = df.shape[1]

    # 3. Compute prevalence per taxon
    prevalence_all = (df >= min_reads).sum(axis=1) / n_samples
    print("=== Prevalence Check (sample) ===")
    sample_taxa = list(prevalence_all.index[:5])
    for taxid in sample_taxa:
        pres = prevalence_all.loc[taxid]
        print(f"Species/Taxid {taxid}: prevalence={pres:.3f}")
    print(f"Total species (taxa): {len(prevalence_all)}")

    # 4. Blacklist: used for threshold estimation only, not for network matrix
    blacklist = [str(t) for t in blacklist]
    blacklist_in_df = [t for t in blacklist if t in df.index]
    if not blacklist_in_df:
        print("❌ No blacklist taxa found in data.")
        return None

    prevalence_black = prevalence_all.loc[blacklist_in_df]
    top_preval_black = prevalence_black.sort_values(ascending=False).head(top_n)
    median_black = top_preval_black.median()

    # 5. Select top N most prevalent species (irrespective of blacklist) for network
    top_taxa = prevalence_all.sort_values(ascending=False).head(top_n).index.tolist()
    filtered_df = df.loc[top_taxa]
    tmp_matrix_path = os.path.join(output, tmp_matrix_name)
    filtered_df.to_csv(tmp_matrix_path, sep="\t", header=True, index=True)

    tmp_edge_path = os.path.join(output, tmp_edge_name)

    try:
        # 6. Run NetCoMi Rscript to generate edge file
        Make_adjacent_matrix(Rscript_path, tmp_matrix_path, tmp_edge_path, coefficient_threshold=0)

        # 7. Compute median correlation from edge file (positive edges only)
        edge_df = pd.read_csv(tmp_edge_path, sep="\t")
        if 'asso' in edge_df.columns:
            pos_corr = edge_df[edge_df['asso'] >= 0]
            if len(pos_corr) == 0:
                median_corr = 0.0
                print("No positive correlations found in edge file, median_corr set to 0.")
            else:
                median_corr = pos_corr['asso'].median()
        else:
            print(f"edge file columns: {list(edge_df.columns)}")
            raise ValueError("No 'asso' column found in edge file. Please check Rscript output or community_detection() parsing.")

        # 8. Compute shift and prevalence threshold
        shift = round(shift_constant * median_corr * 10 / np.sqrt(n_samples), 2)
        threshold = min(max(median_black + shift, 0.0), 0.99)

        print("---------- Find Best prevalence threshold ----------")
        print(f"✅ Samples: {n_samples}")
        print(f"✅ Blacklist taxa found: {len(blacklist_in_df)}")
        print(f"✅ Median of top {top_n} blacklist prevalence = {round(median_black, 3)}")
        print(f"✅ Median Pearson corr (top {top_n} species, >=0 only) = {round(median_corr, 3)}")
        print(f"✅ Shift = shift_constant × median_corr × 10 / sqrt(n) = {round(shift, 3)}")
        print(f"✅ Preliminary prevalence threshold: {threshold}")

        # 9. Count species remaining after threshold
        species_remaining = (prevalence_all >= threshold).sum()
        print(f"✅ Species remaining after threshold {threshold}: {species_remaining}")

        # If too many species remain or threshold==1, return None to trigger min_reads increment
        if species_remaining >= 500 or threshold == 1:
            print("⚠ Threshold leads to too many species (>=500) or threshold==1 → Returning None")
            threshold = None

        if debug:
            print(f"tmp_matrix_path = {tmp_matrix_path}")
            print(f"tmp_edge_path = {tmp_edge_path}")
            print(edge_df.head())

    finally:
        # 10. Remove temporary files
        if os.path.exists(tmp_matrix_path):
            os.remove(tmp_matrix_path)
        if os.path.exists(tmp_edge_path):
            os.remove(tmp_edge_path)

    return threshold







def merging_instrain_compare_output(output):
    file_pattern = os.path.join(output, "IScompare/*/output/*comparisonsTable.tsv")
    file_list = glob.glob(file_pattern)

    dataframes = []

    for file in file_list:
        df = pd.read_csv(file, sep="\t")
        if df.shape[0] > 0:
            dataframes.append(df)

    merged_df = pd.concat(dataframes, ignore_index=True)
    output_file = os.path.join(output, "merged_IS_compare_Table.tsv")
    merged_df.to_csv(output_file, sep="\t", index=False)

    print(f" Successfully merged all comparisonsTable.tsv files into '{output_file}'!")




def get_abundance(output_dir,meta_data,bracken_dir=None):
    if bracken_dir:
        Bracken_dir=bracken_dir
    else:
        Bracken_dir=os.path.join(output_dir, "Bracken_dir")
    report_files=[]
    for sample_id in meta_data:
        sample_dir = os.path.join(Bracken_dir, f"{sample_id}_bracken_species.report")
        report_files.append(sample_dir)

    df_list = []
    for file in report_files:
        df = pd.read_csv(file, sep="\t", header=None, usecols=[0, 3, 4],
                         names=["abundance", "rank", "taxid"], dtype={"taxid": str})
        df = df[df["rank"] == "S"]
        df = df[["taxid", "abundance"]].set_index("taxid")
        df.rename(columns={"abundance": file}, inplace=True)
        df_list.append(df)

    merged_df = pd.concat(df_list, axis=1).fillna(0)
    taxid_mean_abundance = merged_df.mean(axis=1)


    final_df = pd.DataFrame(taxid_mean_abundance, columns=["mean_abundance"])
    final_df.index.name = None

    print(final_df.head(10))

    final_df.to_csv(os.path.join(output_dir,"Abundance.txt"), sep="\t", header=["mean_abundance"], index=True, float_format="%.6f")    
        

