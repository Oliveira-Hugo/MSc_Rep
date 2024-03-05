using CSV
using DataFrames
include("DENV_curation.jl")

#Change user and directories accordingly
fasta_file = "/home/user/input.fasta"
csv_files = [
    "/home/user/GenomeDetective/GenomeDetective/DENV_1/world/BVBRC_genome.csv",
    "/home/user/GenomeDetective/GenomeDetective/DENV_2/world/BVBRC_genome.csv",
    "/home/user/GenomeDetective/GenomeDetective/DENV_3/world/BVBRC_genome.csv",
    "/home/user/GenomeDetective/GenomeDetective/DENV_4/world/BVBRC_genome.csv"
]

#Set the name and location of the output files
output_file = "/home/user/DENV_genomes.fasta"
output_file2 = "/home/user/DENV_serotypes.csv"

BVBRC = Dict{Int, Set{String}}()
for (index, csv) in enumerate(csv_files)
    BVBRC[index] = read_csv(csv)
end

unmatched_IDs = find_unmatched(fasta_file, csv_files)
sequences = read_fasta(fasta_file)
trimmed_seqs = trimm_seqs(sequences, unmatched_IDs)

create_fasta_output(fasta_file, unmatched_IDs, output_file)
create_csv_output(trimmed_seqs, BVBRC, output_file2)
