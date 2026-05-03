import subprocess
import os
def running_kraken(metadata, db, outputdir, core):
    for sample_name in metadata:
        sample1 = metadata[sample_name][1]
        sample2 = metadata[sample_name][2]
        subprocess.run(
            [
                "kraken2",
                "--threads", str(core),
                "--db", db,
                "--memory-mapping",
                "--paired",
                sample1,sample2,
                "--report", os.path.join(outputdir, "Kraken_dir", f'{sample_name}.report'),
                "--output", os.path.join(outputdir, "Kraken_dir",f'{sample_name}.kraken_out')
            ],
            check=True
        )



# do bracken -d /media/junwoojo/18T/standard -i ${i}.report -o ${i}.bracken -r 150 -l 'S' -t 2 ;done



def running_Braken(metadata, db, kraken_dir , bracken_dir , core , read_length):
    for sample_name in metadata:
        subprocess.run(
            [
                "bracken",
                "-d", db,
                "-i", os.path.join(kraken_dir,f'{sample_name}.report'),
                "-w", os.path.join(bracken_dir,f'{sample_name}_bracken_species.report'),
                "-o", os.path.join(bracken_dir, f'{sample_name}.bracken'),
                "-r", str(read_length),
                "-l","S",
                "-t", "2"
            ],
            check=True
        )








#############################################################################################################3



def has_reads_in_kraken_report(report_path):
    try:
        with open(report_path) as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) > 1 and fields[1].isdigit():
                    if int(fields[1]) > 0:
                        return True
    except Exception as e:
        print(f"Error reading {report_path}: {e}")
    return False

def write_dummy_bracken_output(output_file):
    with open(output_file, 'w') as f:
        f.write("name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\tnew_est_reads\tfraction_total_reads\n")
        f.write("Unclassified\t0\tU\t0\t0\t0\t0.0\n")

def write_dummy_species_report(output_file):
    with open(output_file, 'w') as f:
        f.write("100.00\t0\t0\tR\t1\troot\n")

def running_Braken(metadata, db, kraken_dir, bracken_dir, core, read_length):
    for sample_name in metadata:
        report_file = os.path.join(kraken_dir, f'{sample_name}.report')
        output_file = os.path.join(bracken_dir, f'{sample_name}.bracken')
        output_species_file = os.path.join(bracken_dir, f'{sample_name}_bracken_species.report')

        if not os.path.exists(report_file):
            print(f"[SKIP] {sample_name}: Kraken report not found.")
            write_dummy_bracken_output(output_file)
            write_dummy_species_report(output_species_file)
            continue

        if not has_reads_in_kraken_report(report_file):
            print(f"[SKIP] {sample_name}: No reads in Kraken report. Writing dummy outputs.")
            write_dummy_bracken_output(output_file)
            write_dummy_species_report(output_species_file)
            continue

        print(f"[BRACKEN] Running on {sample_name}")
        try:
            subprocess.run(
                [
                    "bracken",
                    "-d", db,
                    "-i", report_file,
                    "-w", output_species_file,
                    "-o", output_file,
                    "-r", str(read_length),
                    "-l", "S",
                    "-t", "2"
                ],
                check=True
            )
        except subprocess.CalledProcessError:
            print(f"[ERROR] Bracken failed on {sample_name}. Writing dummy outputs.")
            write_dummy_bracken_output(output_file)
            write_dummy_species_report(output_species_file)
