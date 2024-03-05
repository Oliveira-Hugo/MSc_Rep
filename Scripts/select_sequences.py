import csv
from Bio import SeqIO

N = 10  #Replace with the desired number of sequences for selection

fasta_file_path = "/home/hugodpo/Documents/LABINFO/Codigos/DENV_genomes.fasta"
csv_file_path = "/home/hugodpo/Documents/LABINFO/Codigos/DENV_serotypes.csv"

selected_fasta_file_path = "tests.fasta"
selected_csv_file_path = "tests.csv"

serotypes = {}
with open(csv_file_path, 'r') as csv_file:
    reader = csv.reader(csv_file)
    next(reader)
    for row in reader:
        serotypes[row[0]] = row[1]

selected_genome_sequences = []
for idx, record in enumerate(SeqIO.parse(fasta_file_path, "fasta")):
    if idx >= N:
        break
    selected_genome_sequences.append(record)

SeqIO.write(selected_genome_sequences, selected_fasta_file_path, "fasta")

with open(selected_csv_file_path, 'w', newline='') as selected_csv_file:
    writer = csv.writer(selected_csv_file)
    writer.writerow(["Sequence ID", "Serotype"])
    for record in selected_genome_sequences:
        sequence_id = record.id
        serotype = serotypes.get(sequence_id, "")
        writer.writerow([sequence_id, serotype])

print("Selected sequences saved on FASTA file:", selected_fasta_file_path)
print("Selected sequences saved on CSV file:", selected_csv_file_path)
