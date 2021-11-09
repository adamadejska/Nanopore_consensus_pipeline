#########################################################################
# Given a fasta file with n reads, create a n x m matrix of kmer frequencies
# where m is the number of kmers. Specify the length of the kmer in the main().
# (We're using 5mers for now) Save the matrix to a csv file for further processing.
# Written by Nick Noll (minor changes by Ada Madejska)
# This script takes two arguments: 1. The path to the folder with file 2. the name of the file
# Each argument needs to be separated by space. No multiple file functionality.
######################################################################### 

module Kmers
using Pkg
#Pkg.add("DataFrames")
using DataFrames
using CSV

struct Sequence
    seq  :: Array{UInt8}
    name :: String
    meta :: String
end

function readfasta(io::IO)
    chan = Channel{Sequence}(0)
    @async begin
        buf=IOBuffer()
        line=readline(io)
        while !isempty(line) && line[1] == '>'
            words      = split(line[2:end])
            name, meta = words[1], join(words[2:end], " ")

            line=readline(io)

            while !isempty(line) && line[1] != '>'
                write(buf,rstrip(line))
                line=readline(io)
            end
            put!(chan, Sequence(take!(buf), name, meta))
        end

        close(buf)
        close(chan)
    end

    return chan
end

const hashtab = Dict{UInt8, UInt64}(
    UInt8('a') => 0, UInt8('A') => 0,
    UInt8('c') => 1, UInt8('C') => 1,
    UInt8('g') => 2, UInt8('G') => 2,
    UInt8('t') => 3, UInt8('T') => 3,
)

function sketch(seq::Array{UInt8}, k::Int)
    # Create a map of hashes corresponding to kmers and counts of each kmer for the particular read.
    fwd :: UInt64 = 0
    rev :: UInt64 = 0

    mask  :: UInt64 = (1 << (2*k)) - 1;
    shift :: UInt64 = 2*(k-1);

    table = Dict{UInt64, Int}()
    l = 0
    for nuc in seq
        c = get(hashtab,nuc,nothing)
        if c === nothing
            l = 0
            continue
        end

        fwd = ((fwd << 2) | c) & mask
        rev = ((rev >> 2) | ((3 ⊻ c) << shift))
        l  += 1

        if l ≥ k
            val = min(fwd, rev)
            if val ∈ keys(table)
                table[val] += 1
            else
                table[val] = 1
            end
        end
    end

    return table
end

function asmatrix(counts, file_name, root_path)
    names  = collect(keys(counts))
    hashes = sort([x for x in reduce(∪, Set(keys(c)) for c in values(counts))])

    # Make a DataFrame to store the hash counts for each read
    df = DataFrame(Reads=names)
    matrix = zeros(Int,length(names),length(hashes))
    col = zeros(Int,1, length(names))
    
    for (i,hash) in enumerate(hashes)
        for (j,name) in enumerate(names)
            col[1,j] = get(counts[name],hash,0)
        end
        col_name = string(hash)
        df[!,col_name] = vec(col)
        col = zeros(Int,1, length(names))
    end

    # Save the DataFrame as a CSV file.
    csv_name = split.(file_name,".")[1] * "_kmer_matrix.csv"
    
    CSV.write(ARGS[3]*'/'*csv_name, df)
    return matrix, names, hashes
end

const root = ARGS[1]
function main(;k=5)

    k = parse(Int64, ARGS[4])
    count = Dict{String, Dict{UInt64, Int}}()

    Threads.@threads for file in readdir(root)
        if endswith(file, ARGS[2])
            # Read each fasta file and for each read in the file calculate the kmer counts
            open("$root/$file") do io
                for record in readfasta(io)
                    count["$file/$(record.name)"] = sketch(record.seq, k)
                end
                m, n, h = asmatrix(count,"$file", "$root")
            end
        end
    end
end

# run main with default settings
main()
end
