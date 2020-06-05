#!/usr/bin/env Rscript
## Usage is canu_read_parser <cutoff> <input_file> <output_file>


library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

cutoff = as.numeric(args[1])
input_file_name = args[2]
output_file_name = args[3]

fasta_file <- readDNAStringSet(input_file_name)

reads <- sapply(strsplit(as.character(names(fasta_file)),' '), "[", 3)
reads_count <- as.numeric(sapply(strsplit(as.character(reads),'='), "[", 2))

writeXStringSet(fasta_file[reads_count >= cutoff], output_file_name)
