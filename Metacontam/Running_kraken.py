import subprocess
import os

# KRAKEN_MMAP=1 adds --memory-mapping. (2026-09-15)
#
# By default kraken2 reads hash.k2d (68.4 GB for the standard database) entirely
# into process memory. On a machine with less RAM than that it is SIGKILLed by the
# OOM killer -- this happened on a 62 GB machine. --memory-mapping switches to mmap
# so a database larger than RAM still runs, and any pages left in the page cache are
# reused by the next sample, which speeds up repeated calls.
#
# The hash access pattern is random, however, so on a spinning disk this can be far
# slower than a sequential read. Enable it only when the database is on NVMe.
#
# Off by default, to keep byte-identical behaviour with the published runs. The
# classification is the same either way; there is no reason to change the setting.
_MMAP = os.environ.get("KRAKEN_MMAP", "") == "1"


def running_kraken(metadata, db, outputdir, core):
    for sample_name in metadata:
        sample1 = metadata[sample_name][1]
        sample2 = metadata[sample_name][2]
        subprocess.run(
            [
                "kraken2",
                "--threads", str(core),
                "--db", db,
            ]
            + (["--memory-mapping"] if _MMAP else [])
            + [
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
    """Run Bracken on each sample.

    Three situations all write a 0-read placeholder, but they mean different
    things. Unless they are recorded separately, "nothing was detected in this
    sample" and "this sample failed" look identical downstream.

        ok                  ran normally
        empty               no reads in the Kraken report -- genuine non-detection
        failed_no_report    the Kraken report is missing -- an earlier step failed
        failed_bracken      Bracken exited with an error -- the run failed

    The status is written to bracken_dir/bracken_status.tsv and summarised at the
    end. The pipeline still continues, since empty samples are normal in
    low-biomass data, but a failure leaves a record.
    """
    status = {}
    for sample_name in metadata:
        report_file = os.path.join(kraken_dir, f'{sample_name}.report')
        output_file = os.path.join(bracken_dir, f'{sample_name}.bracken')
        output_species_file = os.path.join(bracken_dir, f'{sample_name}_bracken_species.report')

        if not os.path.exists(report_file):
            print(f"[FAIL] {sample_name}: Kraken report not found. Writing 0-read placeholder.")
            write_dummy_bracken_output(output_file)
            write_dummy_species_report(output_species_file)
            status[sample_name] = "failed_no_report"
            continue

        if not has_reads_in_kraken_report(report_file):
            print(f"[EMPTY] {sample_name}: No reads in Kraken report (genuine non-detection).")
            write_dummy_bracken_output(output_file)
            write_dummy_species_report(output_species_file)
            status[sample_name] = "empty"
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
            status[sample_name] = "ok"
        except subprocess.CalledProcessError:
            print(f"[FAIL] {sample_name}: Bracken exited with an error. "
                  f"Writing 0-read placeholder.")
            write_dummy_bracken_output(output_file)
            write_dummy_species_report(output_species_file)
            status[sample_name] = "failed_bracken"

    # ---- write the status table ----
    status_path = os.path.join(bracken_dir, "bracken_status.tsv")
    with open(status_path, "w") as fh:
        fh.write("sample\tstatus\n")
        for k in metadata:
            fh.write(f"{k}\t{status.get(k, 'unknown')}\n")

    n_ok     = sum(1 for v in status.values() if v == "ok")
    n_empty  = sum(1 for v in status.values() if v == "empty")
    failed   = {k: v for k, v in status.items() if v.startswith("failed")}
    print(f"\n[BRACKEN SUMMARY] ok {n_ok}  |  empty (no reads detected) {n_empty}  "
          f"|  failed {len(failed)}")
    if failed:
        print("  These samples produced 0-read placeholders because a step FAILED, "
              "not because nothing was detected:")
        for k, v in sorted(failed.items()):
            print(f"    {k}\t{v}")
        print("  They are indistinguishable from empty samples downstream. "
              "Re-run them or drop them from the metadata before interpreting results.")
    print(f"  status written to {status_path}\n")
