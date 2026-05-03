import os
import subprocess
import glob



def bowtie2_index( output_dir):
    '''create bowtie2 index'''
    subprocess.run([
            "bowtie2-build",
            "--quiet",
            os.path.join(output_dir, "Genome_dir", "Candidate.fasta"),
            os.path.join(output_dir, "Genome_dir", "Candidate")],
        check=True)


def mapping(metadata, output_dir, threads):
    output_bam_dir = os.path.join(output_dir,"Bamfiles")
    for sample_name, sample_type, sample1_path, sample2_path in metadata:
        r1 = sample1_path
        r2 = sample2_path
        sam = os.path.join(output_bam_dir,f"{sample_name}.sam")
        bam = os.path.join(output_bam_dir,f"{sample_name}.bam")
        sorted_bam = os.path.join(output_bam_dir,f"{sample_name}.sorted.bam")
        index = os.path.join(output_dir,f"Genome_dir","Candidate")
        print(f"{sample_name} mapping to Candidate.fasta")
        # Run Bowtie2
        subprocess.run([
        "bowtie2", "-p", str(threads), "-x", index, "-1", r1, "-2", r2, "-S", sam
        ])
        print("convert samfile to bamfile")
        # Convert SAM to BAM
        subprocess.run(["samtools", "view", "-bh", "-o", bam, sam])
        print("sorting samfile")
        # Sort BAM
        subprocess.run(["samtools", "sort", "-o", sorted_bam, bam])
        print("indexing sorted bambile")
        # Index BAM
        subprocess.run(["samtools", "index", sorted_bam])

        # Remove intermediate files
        os.remove(sam)
        os.remove(bam)
