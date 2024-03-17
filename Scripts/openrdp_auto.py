import os
import parsl
from parsl.app.app import python_app, bash_app
from parsl.configs.local_threads import config
import re
import random
import csv
import time

parsl.load(config)

@bash_app
def perform_alignment(seq_file, threads, output_file):
    return f'mafft --thread {threads} --auto --inputorder {seq_file} > {output_file}' 

@bash_app
def infer_recombination(input_file, output_file):
    # Docker version:
    # return f'cat {input_file} | docker run --rm -i openrdp > {output_file}'
    # Command-line version:
    return f'openrdp {input_file} -o {output_file}'

def extract_id(header):
    match = re.search(r'^>(\S+)', header)
    if match:
        return match.group(1)
    else:
        return None

def stratified_sampling(N, fasta_alignment, csv_stratum):
    id_strata_map = {}
    with open(csv_stratum, 'r') as csv_file:
        next(csv_file)
        for line in csv_file:
            id_, strata = line.strip().split(',')
            id_strata_map[id_] = strata
    
    strata_counts = {}
    fixed_samples_per_strata = {}
    for id_, strata in id_strata_map.items():
        strata_counts[strata] = strata_counts.get(strata, 0) + 1
    for strata, count in strata_counts.items():
        fixed_samples_per_strata[strata] = min(count, int(N / len(strata_counts)))
    
    remaining_samples = N - sum(fixed_samples_per_strata.values())
    total_count = sum(strata_counts.values())
    samples_per_strata = {strata: fixed_samples_per_strata[strata] + int(remaining_samples * count / total_count) for strata, count in strata_counts.items()}
    
    sampled_ids = {strata: random.sample([id_ for id_, s in id_strata_map.items() if s == strata], count) for strata, count in samples_per_strata.items()}
    
    output_file_name = 'sampled_DENV.fasta'
    with open(fasta_alignment, 'r') as fasta_file, open(output_file_name, 'w') as output_file:
        for line in fasta_file:
            if line.startswith('>'):
                current_id = extract_id(line)
                writing = any(current_id in sampled_ids[strata] for strata in sampled_ids)
            if writing:
                output_file.write(line)
    
    print("\nDENV taxa randomly chosen:")
    print(sampled_ids)
    print("\n\n")
    
    return output_file_name

@python_app
def perform_inference(sampled_fasta, sample_size):
    output_file = f"alignment.fasta_{sample_size}"
    alignment_file = perform_alignment(sampled_fasta, threads=2, output_file=output_file).result()
    print(f'Alignment output status----->{alignment_file}\n')
    alignment_file=output_file
    outfile= f"openrdp_results_{sample_size}"
    recombination_output = infer_recombination(alignment_file, output_file=outfile).result()
    
def main():
    N = 4
    fasta_alignment = "/home/hugodpo/OpenRDP/DENV_genomes.fasta"
    csv_stratum = "/home/hugodpo/OpenRDP/DENV_serotypes.csv"
    openrdp_results = f"/home/hugodpo/OpenRDP/openrdp_results_{N}"

    sampled_fasta = stratified_sampling(N, fasta_alignment, csv_stratum)

    if not os.path.exists(sampled_fasta):
        print(f"Error: Sampled FASTA file '{sampled_fasta}' does not exist.")
   
    start_time = time.time()

    inference_result = perform_inference(sampled_fasta, N).result()
    
    end_time = time.time()
    execution_time = int(end_time - start_time)
    
    with open(openrdp_results, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(f"• Total time of OpenRDP execution with N = {N}: {execution_time} s\n" + content)

    print("-------------------------------------------------------------")
    print(f"\n• Alignment of selected sequences was saved on 'aligned.fasta_{N}'")
    print(f"\n• Results of openrdp inference were saved on 'openrdp_results_{N}'")
    print(f"\n• Total time of OpenRDP execution with N = {N}:", execution_time, "s\n")

main()
