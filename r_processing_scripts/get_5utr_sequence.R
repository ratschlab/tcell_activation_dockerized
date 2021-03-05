
library(data.table)
library(ggplot2)
library(biomaRt)
library(gridExtra)

source("./r_processing_scripts/utils.R")

write_fasta <- function(seq_df, fasta_filename, max_num){
    
    #write into fasta format
    seq_df = seq_df[order(seq_df$padj, decreasing=F), ]
    all_genes_table = seq_df[,"gene_name"]
    unique_genes = names(which(table(all_genes_table) == 1))
    
    utr5_fasta_filename = paste(fasta_filename, "_utr5.fasta", sep="")
    utr3_fasta_filename = paste(fasta_filename, "_utr3.fasta", sep="")
    
    
    # utr5 first
    for(row_idx in 1:max_num){
        if(!is.na(seq_df[row_idx,"utr5_length"]) & 
           seq_df[row_idx,"utr5_length"] < 50000 & 
           seq_df[row_idx,"utr5_length"] > 8 & 
           (seq_df[row_idx,"gene_name"] %in% unique_genes)){
            
            curr_str = paste(">", seq_df[row_idx,"gene_name"], "\n", seq_df[row_idx,"utr5"], sep="")
            if(row_idx == 1){
                write(curr_str,file=utr5_fasta_filename,append=FALSE)
            }else{
                write(curr_str,file=utr5_fasta_filename,append=TRUE)
            }
        }
    }
    
    # utr3 now
    for(row_idx in 1:max_num){
        if(!is.na(seq_df[row_idx,"utr3_length"]) & 
           seq_df[row_idx,"utr3_length"] < 50000 & 
           seq_df[row_idx,"utr3_length"] > 8 & 
           (seq_df[row_idx,"gene_name"] %in% unique_genes)){
            
            curr_str = paste(">", seq_df[row_idx,"gene_name"], "\n", seq_df[row_idx,"utr3"], sep="")
            if(row_idx == 1){
                write(curr_str,file=utr3_fasta_filename,append=FALSE)
            }else{
                write(curr_str,file=utr3_fasta_filename,append=TRUE)
            }
        }
    }
    
}


calc_cutoff_hypergeo <- function(total_genes, ires_genes, ranked_genes, samp_size){
    
    # number of white balls in the urn
    pop_success = length(unique(ires_genes))
    
    # number of black balls in the urn
    total_space = unique(c(total_genes), ires_genes)
    
    pop_fail = sum(! (total_space) %in% (ires_genes))
    
    print("length of all genes")
    print(length(total_space))
    print("number of genes passing")
    print(pop_success)
    print("number of genes failing")
    print(pop_fail)
    
    
    samp_success = sum(ranked_genes[c(1:samp_size)] %in% ires_genes)
    print("IRES genes in enrich")
    print(ranked_genes[which(ranked_genes[c(1:samp_size)] %in% ires_genes)] )
    curr_padj = phyper(samp_success, pop_success, pop_fail, samp_size, lower.tail = FALSE)
    
    enrichment = (samp_success/samp_size) / (pop_success/length(total_genes))
    
    min_pval = c(enrichment, curr_padj, samp_success, samp_size)
    print(min_pval)
    print(c(samp_success, pop_success, pop_fail, samp_size))
    names(min_pval) = c("enrichment", "pval", "samp_success", "samp_size")
    
    min_pval = c(min_pval, ranked_genes[which(ranked_genes[c(1:samp_size)] %in% ires_genes)])
    
    return(min_pval)
    
}

plot_utr_length <- function(seq_df, outname){
    
    # run regression
    seq_df$log2_TE_FC = seq_df$log2FC_TE.DrugTreated.vs.Control.
    a = summary(lm(utr5_length ~ log2_TE_FC, data=seq_df))
    corr_pval_utr5 = round(a$coefficients[2,4], 2)
    
    a = summary(lm(utr3_length ~ log2_TE_FC, data=seq_df))
    corr_pval_utr3 = round(a$coefficients[2,4],2)
    
    
    pdf(paste(outname, "_utr_len.pdf", sep=""))
    gg5 = ggplot(seq_df, aes(x=utr5_length, y=-1*log10(padj))) +
        geom_point() +
        geom_smooth(method='lm',formula=y~x, se=T) +
        annotate("text", x = 2500, y = 2, label = paste("p-val: ", corr_pval_utr5, sep="")) +
        theme_bw() + ggtitle("Correlation between 5'UTR and log2FC TE")

    gg3 = ggplot(seq_df, aes(x=utr3_length, y=-1*log10(padj))) +
        geom_point() +
        geom_smooth(method='lm',formula=y~x, se=T)+
        annotate("text", x = 2500, y = 2, label = paste("p-val: ", corr_pval_utr3, sep="")) +
        theme_bw() + ggtitle("Correlation between 3'UTR and log2FC TE")
    
    gg_all = grid.arrange(gg5, gg3)
    
    print(gg_all)
    dev.off()
}


get_seq_table <- function(res_table){
    
    # get the sequences
    bm <- useMart("ENSEMBL_MART_ENSEMBL", host="www.ensembl.org")
    bm <- useDataset("mmusculus_gene_ensembl", mart=bm)
    
    utr3_tab = getSequence(seqType='3utr',mart=bm,type="ensembl_gene_id",id=res_table$ensembl_gene_id[which(res_table$padj < 0.9999)])
    utr5_tab = getSequence(seqType='5utr',mart=bm,type="ensembl_gene_id",id=res_table$ensembl_gene_id[which(res_table$padj < 0.9999)])
    
    # format the sequences
    colnames(utr3_tab)[1] = "utr3"
    utr3_tab$utr3_length = nchar(utr3_tab$utr3, type = "chars", allowNA = FALSE, keepNA = NA)
    utr3_tab = utr3_tab[order(utr3_tab$utr3_length, decreasing=T),]
    utr3_tab = utr3_tab[which(utr3_tab$utr3 != "Sequence unavailable"),]
    utr3_tab = utr3_tab[which(!duplicated(utr3_tab$ensembl_gene_id)),]
    
    colnames(utr5_tab)[1] = "utr5"
    utr5_tab$utr5_length = nchar(utr5_tab$utr5, type = "chars", allowNA = FALSE, keepNA = NA)
    utr5_tab = utr5_tab[order(utr5_tab$utr5_length, decreasing=T),]
    utr5_tab = utr5_tab[which(utr5_tab$utr5 != "Sequence unavailable"),]
    utr5_tab = utr5_tab[which(!duplicated(utr5_tab$ensembl_gene_id)),]
    
    # merge them together
    merge_tab = merge(res_table, utr5_tab, by="ensembl_gene_id", all.x=TRUE)
    merge_tab = merge(merge_tab, utr3_tab, by="ensembl_gene_id", all.x=TRUE)
    return(merge_tab)
}

#args = commandArgs(trailingOnly=TRUE)
#table_file = args[1]
#ribodiff_file = args[2]
#outdir = args[3]
ribodiff_file = "~/Documents/projects/guido/xpress_tables/transcript_xpressFP_usualRNA_filtered_table_15/result_15_85.txt"
ribodiff_file = "~/Documents/projects/guido/xpress_tables/cds_xpressFP_xpressRNA_filtered_table_15/result_15_85.txt"
outdir = "~/Documents/projects/guido/xpress_tables/cds_xpressFP_xpressRNA_filtered_table_15/"
translate_file = "~/Documents/projects/guido/mm10_gene_trans_name.txt"
ires_ref_file = "~/Documents/projects/guido/annotation/IRES_mouse_human.txt"


# translate the gene names 
translate_df = data.frame(fread(translate_file))
colnames(translate_df) = c("gene_id", "gene_transcript", "gene_name")
translate_df = unique(translate_df[,c("gene_id", "gene_name")])

# read in the ribodiff results
res_table = data.frame(fread(ribodiff_file))
colnames(res_table)[1] = "gene_id"
ribodiff_df_trans = translate_table(res_table, translate_df)
colnames(ribodiff_df_trans)[1] = "ensembl_gene_id"

# now get the 3' and 5' sequences
annot_res_table = get_seq_table(ribodiff_df_trans)

annot_res_table = annot_res_table[order(annot_res_table$padj, decreasing=F),]
write.table(annot_res_table, paste(outdir, "/res_table_utrs.txt", sep=""), quote=F, sep="\t", row.names=F)


# make fasta files for up, down, background to do MEME and TOMTOM
te_up = annot_res_table[annot_res_table$log2FC_TE.DrugTreated.vs.Control. > 0,]
te_down = annot_res_table[annot_res_table$log2FC_TE.DrugTreated.vs.Control. < 0,]

te_bg = te_up[251:nrow(te_up),]
te_bg = rbind(te_bg, te_down[251:nrow(te_down),])


outname_up = paste(outdir, "te_up_padj_top250", sep="")
write_fasta(te_up[1:250,], outname_up, max(nrow(te_up[te_up$padj < 0.1,]), 250))
plot_utr_length(te_up, outname_up)

outname_down = paste(outdir, "te_down_padj_top250", sep="")
write_fasta(te_down[1:250,], outname_down, max(nrow(te_down[te_down$padj < 0.1,]), 250))
plot_utr_length(te_down, outname_down)

outname_bg = paste(outdir, "te_background_padj_bg500", sep="")
write_fasta(te_bg, outname_bg, min(nrow(te_bg), 500))


# now calculate enrichment for IRES
ires_genes = c(unlist(data.frame(fread(ires_ref_file, header=F))))

total_genes = annot_res_table[!is.na(annot_res_table$padj), "gene_name"]
ranked_genes_up = te_up[1:50,"gene_name"]
calc_cutoff_hypergeo(total_genes, ires_genes, ranked_genes=ranked_genes_up, samp_size=length(ranked_genes_up))

ranked_genes_down = te_down[1:250,"gene_name"]
calc_cutoff_hypergeo(total_genes, ires_genes, ranked_genes=ranked_genes_down, samp_size=length(ranked_genes_down))


