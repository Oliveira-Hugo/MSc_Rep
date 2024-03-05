function read_csv(csv_file)
    data = CSV.read(csv_file, DataFrame)
    return Set(filter(!ismissing, data[!, "GenBank Accessions"]))
end
    
function read_fasta(file::String)
    ids = Set{String}()
    open(file) do f
        for line in eachline(f)
            if startswith(line, '>')
                id = replace(strip(line), ">" => "")
                push!(ids, id)
            end
        end
    end
    return ids
end

function find_unmatched(fasta_file, csv_files)
    all_IDs = Set{String}()
    for csv_file in csv_files
        union!(all_IDs, read_csv(csv_file))
    end

    fasta_IDs = read_fasta(fasta_file)
    unmatched_IDs = setdiff(fasta_IDs, all_IDs)
    return(unmatched_IDs)
end

function trimm_seqs(sequences, unmatched_IDs)
    trimmed_sequences = Set{String}()
    for seq in sequences
        id = first(split(seq, '|'))
        if id in unmatched_IDs
            continue
        end
        push!(trimmed_sequences, seq)
    end
    return trimmed_sequences
end


function create_fasta_output(fasta_file, unmatched_seqs, output_file)
    input_stream = open(fasta_file, "r")
    output_stream = open(output_file, "w")

    current_id = ""
    current_seq = ""

    for line in eachline(input_stream)
        if startswith(line, '>')
            if !isempty(current_id) && !(current_id in unmatched_seqs)
                println(output_stream, ">$current_id")
                println(output_stream, current_seq)
            end
            current_id = replace(strip(line), ">" => "")
            current_seq = ""
        else
            current_seq *= line
        end
    end

    if !isempty(current_id) && !(current_id in unmatched_seqs)
        println(output_stream, ">$current_id")
        println(output_stream, current_seq)
    end

    close(input_stream)
    close(output_stream)
end
 
function create_csv_output(trimmed_seqs, BVBRC::Dict{Int, Set{String}}, output_file::String)
    result_df = DataFrame(ID = String[], Serotype = String[])
    for id in trimmed_seqs
        serotype_found = false
        for (serotype, genbank_accessions) in BVBRC
            if id in genbank_accessions
                push!(result_df, (id, "DENV_$serotype"))
                serotype_found = true
                break
            end
        end
        if !serotype_found
            println("Warning: ID $id not found in any BVBRC serotype.")
        end
    end
    CSV.write(output_file, result_df)
end
