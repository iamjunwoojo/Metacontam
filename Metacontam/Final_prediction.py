import pandas as pd
import numpy as np
from collections import defaultdict
import sys
import os

def process_df_with_blacklist(df, base_percentile, blacklist):
    taxid_to_conANI = defaultdict(list)

    for _, row in df.iterrows():
        if row['compared_bases_count'] > 50:
            taxid = str(row['scaffold'].split('|')[1])  # must be string
            con_ani = row['conANI']
            # Prevent NA values
            try:
                if not pd.isnull(con_ani):
                    taxid_to_conANI[taxid].append(float(con_ani))
            except Exception:
                continue

    # Compute mean ANI per taxid
    mean_con = {taxid: np.mean(values) for taxid, values in taxid_to_conANI.items() if values}
    mean_con_values = list(mean_con.values())

    # Sets of all species and blacklist species
    all_species = set(taxid_to_conANI.keys())
    blacklist = set(str(x) for x in blacklist)  # convert all to string
    bl_species = all_species & blacklist
    bl_ratio = len(bl_species) / len(all_species) if all_species else 0

    # Compute effective percentile
    effective_percentile = base_percentile - 0.6 * bl_ratio
#    if effective_percentile < 0.01:
#        effective_percentile = 0.01

    threshold = np.percentile(mean_con_values, effective_percentile * 100)

    # Sort by mean ANI descending
    sorted_con = sorted(mean_con.items(), key=lambda x: x[1], reverse=True)
    return sorted_con, threshold, bl_ratio, effective_percentile




def final_prediction(output_dir, kraken_db, blacklist_list, base_percentile=0.6):
    print("\n\n---------- Final prediction Start ----------\n")
    IS_compare_file = os.path.join(output_dir, "merged_IS_compare_Table.tsv")
    if not os.path.exists(IS_compare_file):
        print(f"❌ File not found: {IS_compare_file}")
        sys.exit(1)

    IS_compare_df = pd.read_csv(IS_compare_file, sep="\t")
    sorted_con, threshold, bl_ratio, effective_percentile = process_df_with_blacklist(
        IS_compare_df, base_percentile, blacklist_list
    )

    Final_predicted_taxa = [i for i, l in sorted_con if l >= threshold]
    print(f'Total species: {len(sorted_con)}')
    print(f'Blacklist species: {len(set(blacklist_list) & set([i for i, _ in sorted_con]))}')
    print(f'Blacklist ratio: {bl_ratio:.4f}')
    print(f'Effective percentile: {effective_percentile:.4f}')
    print(f'Threshold ANI value: {threshold:.4f}')
    print(f'{len(Final_predicted_taxa)} species were classified as contaminants.')

    # Build result dataframe
    taxa = [i for i, l in sorted_con]
    scores = [l for i, l in sorted_con]
    df = pd.DataFrame({
        "Taxid": taxa,
        "Mean-Pairwise-ANI": scores
    })

    df["contamination_status"] = df["Mean-Pairwise-ANI"].apply(
        lambda x: "Contaminant" if x >= threshold else "Non-contaminant"
    )

    # Attach scientific names from names.dmp
    names_dmp_path = os.path.join(kraken_db, "taxonomy", "names.dmp")
    if not os.path.exists(names_dmp_path):
        print(f"❌ File not found: {names_dmp_path}")
        sys.exit(1)

    names_df = pd.read_csv(
        names_dmp_path,
        sep="\t\|\t",
        header=None,
        names=["Taxid", "name_txt", "unique_name", "name_class"],
        engine="python"
    )
    names_df['name_class'] = names_df['name_class'].str.strip("\t|")
    names_df['Taxid'] = names_df['Taxid'].astype(str)
    sci_names = names_df[names_df["name_class"] == "scientific name"].copy()
    sci_names["Taxid"] = sci_names["Taxid"].astype(str)
    df["Taxid"] = df["Taxid"].astype(str)
    df = df.merge(sci_names[["Taxid", "name_txt"]], on="Taxid", how="left")
    df.rename(columns={"name_txt": "scientific_name"}, inplace=True)

    # Save results
    out_path = os.path.join(output_dir, "Final_prediction.txt")
    df.to_csv(out_path, sep="\t", index=False)
    print(f"\nResult saved: {out_path}")
    print("\n---------- Final prediction End ----------\n")
    return df

# -------------------------------
# Usage example (when used as a module)
# -------------------------------
# if __name__ == "__main__":
#     output_dir = "/your/output/path"
#     kraken_db = "/your/kraken/db/path"
#     blacklist_list = ["12345", "67890", ...]  # list of taxid strings
#     final_prediction(
#         output_dir=output_dir,
#         kraken_db=kraken_db,
#         blacklist_list=blacklist_list,
#         base_percentile=0.6
#     )

