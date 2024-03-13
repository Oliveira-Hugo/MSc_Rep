# MSc_Rep
Repository for master's script and data storage

# General objective
Submit an alignment fasta dataset of concatenated DENV genes from all serotypes and regions to OpenRDP for recombination inference. 

# Specific objectives
1- Check the co-occurrence of output.nex DENV taxa on BVBRC_genome.csv files;  
2- Test the [OpenRDP implementation](https://github.com/PoonLab/OpenRDP);  
3- Try running OpenRDP with the whole DENV dataset. If the program is unable to do so, try with increasing random stratified sampling;  
4- Implement pipeline for i) handling data, ii) perform stratified sampling, iii) align sampled seqs., iv) perform recombination inference;  
5- Create script for curate the OpenRDP result and save it on a separated file.  

# Description
1) The "main.jl" function calls functions defined on "DENV_curation.jl" and process DENV data as briefly listed on (1);  
2) "BVBRC_x.csv" are inputs of the script, but are not present on this repository for storage reasons;   
3) "DENV_genomes.fasta" and "DENV_serotypes.csv" are the outputs of running main.jl with our dataset;  
4) The parsl script described on Specific Objectives item 4 is named "openrdp_auto.py";  
5) The file for processing OpenRDP results is named "trimming.py".  

# Current objective(s)
* **Adapt the "openrdp_auto.py" pipeline for running via Docker on SDumont**
