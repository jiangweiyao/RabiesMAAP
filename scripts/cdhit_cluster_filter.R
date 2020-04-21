#!/usr/bin/env Rscript
## Usage is cdhit_cluster_filter.R <cutoff> <fasta_file> <clstr_file> <output_file>

library(Biostrings)

args = commandArgs(trailingOnly=TRUE)

cutoff = as.numeric(args[1])
fasta_file = args[2]
clstr_file = args[3]
output_file = args[4]

clstr <- read.csv(clstr_file, sep = "\t", row.names = NULL, header = FALSE, stringsAsFactors = FALSE)
clstr2 <- clstr
n = nrow(clstr)
x = 0
numbers_only <- function(x) !grepl("\\D", x)
for (row in c(1:n)) {
    if (numbers_only(clstr2[row,1]) == TRUE) {
        clstr2[row,1] <- x}
    else {NULL}
    x <- clstr2[row,1]
}

clstr2 <- subset(clstr2, V2 !="")

reads <- sapply(strsplit(clstr2$V2,' '), "[", 2)
reads1 <- substr(reads, 2, nchar(reads)-3)
alignment <- sapply(strsplit(clstr2$V2,' '), "[", 4)

clstr_final <- data.frame(Cluster = clstr2$V1, Reads = reads1, Alignment = alignment)
clstr_summary <- subset(clstr_final, is.na(Alignment))

clstr_table <- merge(clstr_summary, data.frame(table(clstr_final$Cluster)), by.x = "Cluster", by.y = "Var1")
clstr_cutoff <- subset(clstr_table, Freq >= cutoff)

write.csv(clstr_table, paste0(clstr_file, ".csv"), row.names = FALSE) 

fasta <- readDNAStringSet(fasta_file)
fasta_name <- sapply(strsplit(as.character(names(fasta)),' '), "[", 1)

writeXStringSet(fasta[fasta_name %in% clstr_cutoff$Reads], output_file)

