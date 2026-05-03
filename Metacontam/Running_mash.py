import numpy as np
import pandas as pd
import random


import subprocess
import os

def run_mash(meta_data, output, THREADS, kmer_size=21, min_kmer_copy_num=2, sketch_size=100000):
    sketch_dir = os.path.join(output,"mash_sketches")
     
    if not os.path.exists(sketch_dir):
        os.mkdir(sketch_dir)

    for sample_id in meta_data:
        sketch_file = os.path.join(sketch_dir, sample_id)
        temp_fastq = os.path.join(output, "temp.fastq")

        read_file_1 = meta_data[sample_id][1]
        read_file_2 = meta_data[sample_id][2]
        
        os.system(f"cat {read_file_1} {read_file_2} > {temp_fastq}")

        subprocess.run(["mash", "sketch",
                        "-s", str(sketch_size),
                        "-k", str(kmer_size),
                        "-m", str(min_kmer_copy_num),
                        "-r", temp_fastq,
              #          "-p", str(THREADS),
                        "-o", sketch_file])

        if os.path.exists(temp_fastq):
            os.remove(temp_fastq)




def calc_mash_dist(meta_data, output):
    sketch_dir = os.path.join(output, "mash_sketches")
    sample_list = list(meta_data.keys())

    # Build distance matrix
    with open(os.path.join(output, "dist_matrix.txt"), "w") as dist_matrix:
        # Write header
        dist_matrix.write("\t" + "\t".join(sample_list) + "\n")

        for sample_id1 in sample_list:
            dist_matrix.write(sample_id1)  # Write row name
            for idx, sample_id2 in enumerate(sample_list):
                dry_run = subprocess.check_output(["mash", "dist",
                                                   os.path.join(sketch_dir, f"{sample_id1}.msh"),
                                                   os.path.join(sketch_dir, f"{sample_id2}.msh")])

                lines = dry_run.splitlines()
                dist = lines[0].decode("utf-8").split("\t")[2]

                # No tab after last column
                if idx == len(sample_list) - 1:
                    dist_matrix.write(f"\t{dist}")
                else:
                    dist_matrix.write(f"\t{dist}")
            dist_matrix.write("\n")  # Move to next row





def bin_based_stratified_sampling(input_mt, num_bins=10, total_samples_pair=1000, seed=None):
    # Set the seed: if not provided, default to 4
    if seed is None:
        seed = 4
    random.seed(seed)
    np.random.seed(seed)

    # Extract the upper triangle of the matrix, excluding the diagonal
    upper_triangle = input_mt.where(np.triu(np.ones(input_mt.shape), k=1).astype(bool))
    upper_triangle_values = upper_triangle.stack()

    # Convert the upper triangle values to a sorted list of tuples
    sorted_values = sorted(upper_triangle_values.items(), key=lambda x: x[1], reverse=True)
    
    # If the number of columns in input_mt is 50 or less, return all pairs without stratified sampling
    if input_mt.shape[1] < 46:
        print("sample pairs < 1000 ...  Skipping stratified sampling")
        return [(row, col, value) for (row, col), value in sorted_values]

    df = pd.DataFrame([(row, col, value) for (row, col), value in sorted_values],
                      columns=["row", "col", "value"])

    # Divide the values into bins (histogram-based)
    num_bins = num_bins  # Number of bins
    df["bin"] = pd.cut(df["value"], bins=np.linspace(0, 1, num_bins + 1), include_lowest=True)

    # Determine the number of samples per bin, proportional to bin counts
    bin_counts = df["bin"].value_counts()
    bin_sample_counts = (bin_counts / bin_counts.sum() * total_samples_pair).astype(int)

    # Distribute remaining samples to bins with the most data
    remaining_samples = total_samples_pair - bin_sample_counts.sum()
    if remaining_samples > 0:
        fractional_bins = bin_counts / bin_counts.sum() * remaining_samples
        fractional_bins = fractional_bins.round().astype(int)
        bin_sample_counts += fractional_bins

    # Ensure exact total by adjusting the largest bin if necessary
    adjustment = total_samples_pair - bin_sample_counts.sum()
    if adjustment != 0:
        largest_bin = bin_sample_counts.idxmax()
        bin_sample_counts[largest_bin] += adjustment

    # Sample data from each bin
    sampled_data = []
    for bin_interval, count in bin_sample_counts.items():
        bin_data = df[df["bin"] == bin_interval]
        sampled_data.append(bin_data.sample(n=min(count, len(bin_data)), replace=False, random_state=seed))

    # Combine sampled data and return the result as a list of tuples
    sampled_df = pd.concat(sampled_data)
    return [(row, col, value) for row, col, value in sampled_df[["row", "col", "value"]].to_records(index=False)]


