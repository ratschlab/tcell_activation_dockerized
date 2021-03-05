### this is a preprocessing script to generate multiple filters
library("data.table")

filter_counts <- function(in_table, outdir, filter_amount, outname){
    
    in_table_rp_sum = rowSums(in_table[,grep("RP|FP", colnames(in_table))])
    in_table_rna_TF = apply(in_table[,grep("RNA", colnames(in_table))], 1, function(x) all(x >= 50))
    
    idx_pass = intersect(which(in_table_rp_sum >= filter_amount), which(in_table_rna_TF))
    
    outfile = paste(outdir, outname, "_filtered_table_", filter_amount, ".txt", sep="")
    write.table(in_table[idx_pass,], outfile, row.names=F, quote=F, sep="\t")
    
}

read_in_table <- function(table_file){
    table_df = data.frame(fread(table_file))
    return(table_df)
}

#args = commandArgs(trailingOnly=TRUE)
#table_file = args[1]
#ribodiff_file = args[2]
#outdir = args[3]
outdir = "~/Documents/projects/guido/xpress_tables/"
filter_level = 15


table_file = "~/Documents/projects/guido/xpress_tables/transcript_xpressFP_usualRNA.tsv"
outname = "transcript_xpressFP_usualRNA"
in_table = read_in_table(table_file)
filter_counts(in_table, outdir, filter_level, outname)


table_file = "~/Documents/projects/guido/xpress_tables/cds_xpressFP_previousRNA.tsv"
outname = "cds_xpressFP_previousRNA"
in_table = read_in_table(table_file)
filter_counts(in_table, outdir, filter_level, outname)


table_file = "~/Documents/projects/guido/xpress_tables/transcript_xpressFP_xpressRNA.tsv"
outname = "transcript_xpressFP_xpressRNA"
in_table = read_in_table(table_file)
filter_counts(in_table, outdir, filter_level, outname)


table_file = "~/Documents/projects/guido/xpress_tables/cds_xpressFP_xpressRNA.tsv"
outname = "cds_xpressFP_xpressRNA"
in_table = read_in_table(table_file)
filter_counts(in_table, outdir, filter_level, outname)
