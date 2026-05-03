import os
import subprocess
import pandas as pd
import shutil


#def Instrain_profile(metadata, output_dir, output_bam_dir, candidate_fasta, threads):
#    if output_bam_dir:
#        output_BAM_dir =  output_bam_dir
#    else:
#        output_BAM_dir = os.path.join(output_dir,"Bamfiles")
#    output_IS_dir = os.path.join(output_dir,"ISfiles")
#    if not candidate_fasta:
#        candidate_fasta =  os.path.join( output_dir ,"Genome_dir","Candidate.fasta")
#    for sample_name, sample_type, sample1_path, sample2_path in metadata:
#        subprocess.run([
#                "inStrain",
#                "profile",
#                "--skip_mm_profiling",
#                os.path.join(output_BAM_dir,f"{sample_name}.sorted.bam"),
#                "--skip_plot_generation",
#                "--skip_genome_wide",
#                "-c", str(0),
#                candidate_fasta,
#                "-o", os.path.join(output_IS_dir,f"{sample_name}.IS"),
#                "-p",str(threads)
#                ]
#                ,check=True)


def Instrain_profile(metadata, output_dir, output_bam_dir, candidate_fasta, threads):
    if output_bam_dir:
        output_BAM_dir = output_bam_dir
    else:
        output_BAM_dir = os.path.join(output_dir, "Bamfiles")
    output_IS_dir = os.path.join(output_dir, "ISfiles")
    os.makedirs(output_IS_dir, exist_ok=True)

    if not candidate_fasta:
        candidate_fasta = os.path.join(output_dir, "Genome_dir", "Candidate.fasta")

    for sample_name, sample_type, sample1_path, sample2_path in metadata:
        is_output_dir = os.path.join(output_IS_dir, f"{sample_name}.IS")

        try:
            subprocess.run([
                "inStrain", "profile",
                "--skip_mm_profiling",
                os.path.join(output_BAM_dir, f"{sample_name}.sorted.bam"),
                "--skip_plot_generation",
                "--skip_genome_wide",
                "-c", "0",
                candidate_fasta,
                "-o", is_output_dir,
                "-p", str(threads)
            ], check=True)
        except subprocess.CalledProcessError:
            print(f"[Error] inStrain profile failed for {sample_name}")
            if os.path.exists(is_output_dir):
                shutil.rmtree(is_output_dir)


#def Instrain_compare(output_dir, candidate_fasta, threads):
#    stratified_sampling_file = pd.read_csv(os.path.join(output_dir,"pair_output.tsv"),sep="\t")
#    output_IS_dir = os.path.join(output_dir,"ISfiles")
#    if not candidate_fasta:
#        candidate_fasta =  os.path.join( output_dir ,"Genome_dir","Candidate.fasta")
#    for row in stratified_sampling_file.iterrows():
#        try:
#            sample1=row[1]['sample1']
#            sample2=row[1]['sample2']
#            subprocess.run([
#                "inStrain",
#                "compare",
#                "-p",str(threads),
#                "--breadth","0",
#                "-c","0",
#                "-f","0",
#                "-sc",candidate_fasta,
#                "--skip_plot_generation",
#                "-i", os.path.join(output_IS_dir,f"{sample1}.IS"), os.path.join(output_IS_dir,f"{sample2}.IS"),
#                "-o", os.path.join(output_dir,"IScompare",f"{sample1}_{sample2}_output")
#                ]
#                ,check=True)
#        except subprocess.CalledProcessError:
#            print(f"[Error] inStrain compare failed for {sample1} vs {sample2}")
#            # Remove failed output directory
#            if os.path.exists(compare_output_dir):
#                shutil.rmtree(compare_output_dir)




def Instrain_compare(output_dir, candidate_fasta, threads):
    stratified_sampling_file = pd.read_csv(os.path.join(output_dir, "pair_output.tsv"), sep="\t")
    output_IS_dir = os.path.join(output_dir, "ISfiles")
    if not candidate_fasta:
        candidate_fasta = os.path.join(output_dir, "Genome_dir", "Candidate.fasta")

    os.makedirs(os.path.join(output_dir, "IScompare"), exist_ok=True)

    for _, row in stratified_sampling_file.iterrows():
        sample1 = row['sample1']
        sample2 = row['sample2']
        compare_output_dir = os.path.join(output_dir, "IScompare", f"{sample1}_{sample2}_output")

        try:
            subprocess.run([
                "inStrain", "compare",
                "-p", str(threads),
                "--breadth", "0",
                "-c", "0",
                "-f", "0",
                "-sc", candidate_fasta,
                "--skip_plot_generation",
                "-i", os.path.join(output_IS_dir, f"{sample1}.IS"),
                      os.path.join(output_IS_dir, f"{sample2}.IS"),
                "-o", compare_output_dir
            ], check=True)

        except subprocess.CalledProcessError:
            print(f"[Error] inStrain compare failed for {sample1} vs {sample2}")
            if os.path.exists(compare_output_dir):
                shutil.rmtree(compare_output_dir)
