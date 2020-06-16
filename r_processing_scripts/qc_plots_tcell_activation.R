
require(ggplot2)
require(data.table)
require(ggrepel)
require(gridExtra)

translate_ENS_to_HGNC <- function(genes_to_trans, translate_df, bind_col, new_col){
    
    new_df = merge(unique(translate_df[,c(bind_col, new_col)]), unique(genes_to_trans), by=bind_col)
    return(new_df)
    
    
}

translate_table <-function(ribodiff_df, translate_df){
    genes_to_trans = ribodiff_df
    print(head(genes_to_trans))
    print(dim(genes_to_trans))
 
    trans_genes = translate_ENS_to_HGNC(genes_to_trans, translate_df, bind_col="gene_id", new_col="gene_name")
    trans_genes = trans_genes[ order(trans_genes$padj), ]
    
    return(trans_genes)
}


normalize_counts <- function(in_df, lib_size){
    # remove gene_names
    in_matr = in_df[,2:ncol(in_df)]
    
    # order the libsize file correctly
    row.names(lib_size) = lib_size$Samples
    lib_size = lib_size[colnames(in_matr),]
    
    # normalize
    in_matr = t(t(in_matr) / lib_size$library_size)
    
    norm_df = data.frame(gene_id=in_df$gene_id, in_matr)
    return(norm_df)
}

get_melted_merged_table <- function(raw_counts_df, norm_counts_df, count_type_df){
    
    raw_counts_df_melt = melt(raw_counts_df)
    colnames(raw_counts_df_melt) = c("gene", "Samples", "raw_count")
    norm_counts_df_melt = melt(norm_counts_df)
    colnames(norm_counts_df_melt) = c("gene", "Samples", "norm_count")
    full_df = merge(raw_counts_df_melt, norm_counts_df_melt)
    full_df = merge(full_df, count_type_df)
    
    return(full_df)
    
}

plot_boxplot <- function(full_df, col_interest, title){
    
    plot_df = full_df
    colnames(plot_df)[colnames(plot_df) == col_interest] = "count"
    
    gg = ggplot(plot_df, aes(x=Samples, y=log10(count+1), fill=count_type)) +
        geom_boxplot() +
        theme_bw() + scale_fill_brewer(palette="Dark2") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(fill = "Count Type") +
        ggtitle(title)
    
    return(gg)
}


plot_corr <- function(full_df, use_normalized, title){
    
    # to do dcast choose the value variable
    cast_var = "norm_count"
    if(!use_normalized){
        cast_var = "raw_count"
    }
    corr_raw_df = dcast(full_df, sample_ids + gene + condition ~ count_type, value.var="raw_count")                          
    corr_raw_df = corr_raw_df[corr_raw_df$sample_ids != "Sample_2hr_1",]                         
    
    get_r2 <- function(in_df){
        a = summary(lm(RNA~RP, in_df))
        return(a$adj.r.squared)
    }
    raw_r2 = by(corr_raw_df, list(corr_raw_df$sample_ids), get_r2)
    raw_r2_df <- data.frame(raw_r2 = paste("Adj. R-Sqr.=", round(as.numeric(unclass(raw_r2)), 3)), sample_ids = names(raw_r2))
    
    gg_raw = ggplot(corr_raw_df, aes(x=RP, y=RNA)) + geom_point() + geom_smooth(method = lm) +
        facet_wrap(~ sample_ids) + theme_bw() + scale_color_brewer(palette="Dark2") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(title) +
        geom_text(data = raw_r2_df, aes(x = 90, y = 100000, label = raw_r2))
    
    return(gg_raw)
    
}

plot_pca <- function(in_df, count_type, cond_vec, title_text){
    pca_res = prcomp(t(in_df[,count_type_df$Samples[count_type_df$count_type == count_type]]), scale. = TRUE)
    pca_raw_RP = data.frame(pca_res$x)
    pca_raw_RP$condition = cond_vec
    pca_raw_RP$Sample = row.names(pca_raw_RP)
    
    percentage <- round(pca_res$sdev / sum(pca_res$sdev) * 100, 2)
    percentage <- paste( colnames(pca_raw_RP), "(", paste( as.character(percentage), "%", ")", sep="") )
    
    gg_pca_raw_RP = ggplot(pca_raw_RP,aes(x=PC1,y=PC2,color=condition, label=Sample))
    gg_pca_raw_RP = gg_pca_raw_RP+geom_point()  +
        theme_bw() + xlab(percentage[1]) + ylab(percentage[2]) +
        ggtitle(title_text) + geom_text_repel()
    
    return(gg_pca_raw_RP)
    
}
plot_volcano <- function(norm_counts_df, gene_names_df, ribodiff_df){
    
    ribodiff_df_trans = translate_table(ribodiff_df, translate_df)
    
    mean_expr_RP_T0 = rowMeans(norm_counts_df[,c("Sample_T0_1_RP",   "Sample_T0_2_RP",   "Sample_T0_3_RP")])
    mean_expr_RP_2hr = rowMeans(norm_counts_df[,c("Sample_2hr_2_RP",  "Sample_2hr_3_RP")])
    mean_expr_RP = data.frame(gene_id = norm_counts_df$gene, T0_mean = mean_expr_RP_T0, hr2_mean = mean_expr_RP_2hr)
    mean_expr_RP = merge(mean_expr_RP, ribodiff_df_trans)
    
    # add plotting columns 
    mean_expr_RP$hist_color = "Nonsignificant"
    mean_expr_RP$hist_color[mean_expr_RP$padj < 0.1 & mean_expr_RP$log2FC_TE.DrugTreated.vs.Control. > 0] = "Pos. Significant"
    mean_expr_RP$hist_color[mean_expr_RP$padj < 0.1 & mean_expr_RP$log2FC_TE.DrugTreated.vs.Control. < 0] = "Neg. Significant"
    mean_expr_RP$hist_color = factor(mean_expr_RP$hist_color, levels=c("Neg. Significant", "Nonsignificant", "Pos. Significant"))
    
    mean_expr_RP$gene_name[mean_expr_RP$padj > 0.1] = ""
    
    
    gg_T0_expr = ggplot(mean_expr_RP, aes(x=log2FC_TE.DrugTreated.vs.Control., y=-log10(padj), color=T0_mean>5, label=gene_name)) + 
        geom_point() + geom_hline(yintercept=1, color="red") + geom_label_repel() +
        theme_bw() + scale_color_brewer(palette="Dark2") + labs(color = "Mean Normalized\nRP counts T0 > 5") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlim(c(-9,9)) + ylim(c(0, 3.5)) +
        ggtitle("RiboDiff Result") + xlab("TE log2FC T0 vs 2hr")
    
    gg_2hr_expr = ggplot(mean_expr_RP, aes(x=log2FC_TE.DrugTreated.vs.Control., y=-log10(padj), color=hr2_mean>5, label=gene_name)) + 
        geom_point() + geom_hline(yintercept=1, color="red") + geom_label_repel() +
        theme_bw() + scale_color_brewer(palette="Dark2") + labs(color = "Mean Normalized\nRP counts 2HR > 5") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlim(c(-9,9)) + ylim(c(0, 3.5)) +
        ggtitle("RiboDiff Result") + xlab("TE log2FC T0 vs 2hr")

    gg_expr = ggplot(mean_expr_RP, aes(x=log2FC_TE.DrugTreated.vs.Control., y=-log10(padj), color=-log10(padj), label=gene_name)) + 
        geom_point() + geom_hline(yintercept=1, color="red") + geom_label_repel() +
        theme_bw() + scale_color_gradient2(midpoint=1, low="black", mid="blue", high="red") + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlim(c(-9,9)) + ylim(c(0, 3.5)) +
        ggtitle("RiboDiff Result") + xlab("TE log2FC T0 vs 2hr")
    
     gg_log2FC_expr = ggplot(mean_expr_RP, aes(x=log2FC_TE.DrugTreated.vs.Control., fill=hist_color)) + 
        geom_histogram(position = "identity", bins = 100) +
        theme_bw() + scale_fill_manual(values=c("#56B4E9", "#999999", "#E69F00")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
         ggtitle("RiboDiff Result") + xlab("TE log2FC T0 vs 2hr")
     
     gg_all = grid.arrange(gg_T0_expr, gg_2hr_expr, gg_expr, gg_log2FC_expr, ncol=2)
    
    return(gg_all)
    
}
    
#args = commandArgs(trailingOnly=TRUE)
#table_file = args[1]
#ribodiff_file = args[2]
#outdir = args[3]
table_file = "~/Documents/projects/guido/total_table.txt"
ribodiff_file = "~/Documents/projects/guido/result_15_85.txt"
translate_file = "~/Documents/projects/guido/mm10_gene_trans_name.txt"
outdir = "~/Documents/projects/guido/new_results/"

raw_counts_df = data.frame(fread(table_file))
colnames(raw_counts_df)[1] = "gene_id"

translate_df = data.frame(fread(translate_file))
colnames(translate_df) = c("gene_id", "gene_transcript", "gene_name")
translate_df = unique(translate_df[,c("gene_id", "gene_name")])

ribodiff_df = data.frame(fread(ribodiff_file))
colnames(ribodiff_df)[1] = "gene_id"
ribodiff_df = ribodiff_df[,1:7]
ribodiff_df = ribodiff_df[ribodiff_df$pval != "NaN",]


# only get counts for genes tested
raw_counts_df = raw_counts_df[raw_counts_df$gene %in% ribodiff_df$gene,]


# count_type definition
count_type_df = data.frame(Samples = c("Sample_T0_1_RP",   "Sample_T0_2_RP",   "Sample_T0_3_RP",   
                                       "Sample_2hr_2_RP",  "Sample_2hr_3_RP",  "Sample_T0_1_RNA",  
                                       "Sample_T0_2_RNA",  "Sample_T0_3_RNA",  "Sample_2hr_1_RNA",
                                       "Sample_2hr_2_RNA", "Sample_2hr_3_RNA"),
                           sample_ids = c("Sample_T0_1",   "Sample_T0_2",   "Sample_T0_3",   
                                          "Sample_2hr_2",  "Sample_2hr_3",  "Sample_T0_1",  
                                          "Sample_T0_2",  "Sample_T0_3",  "Sample_2hr_1",
                                          "Sample_2hr_2", "Sample_2hr_3"),
                           count_type = c(rep("RP",5), rep("RNA", 6)),
                           condition = c(rep("T0", 3), rep("2hr", 2), rep("T0", 3), rep("2hr", 3)))

# normalize
library_size = data.frame(Samples = c("Sample_T0_1_RP",   "Sample_T0_2_RP",   "Sample_T0_3_RP",   
                                      "Sample_2hr_2_RP",  "Sample_2hr_3_RP",  "Sample_T0_1_RNA",  
                                      "Sample_T0_2_RNA",  "Sample_T0_3_RNA",  "Sample_2hr_1_RNA",
                                      "Sample_2hr_2_RNA", "Sample_2hr_3_RNA"),
                          library_size = c(#0.653,  0.758,  0.758,  0.85,   0.871,
                                           #1.149,  1.132,  1.172,  0.822,  0.957,  0.894
                                            0.571,  1,     1,     1.286,  1.493,
                                            1.463,  0.964,  1.461,  0.649,  0.854,  1.036))
norm_counts_df = normalize_counts(raw_counts_df, library_size)

######################                       
# now make plots
######################                       

# plot PCA 
gg_pca_raw_rp = plot_pca(raw_counts_df, count_type="RP", cond_vec=c(rep("T0", 3), rep("2hr", 2)), title_text="PCA Raw RP Counts")
gg_pca_norm_rp = plot_pca(norm_counts_df, count_type="RP", cond_vec=c(rep("T0", 3), rep("2hr", 2)), title_text="PCA Normalized RP Counts")
gg_pca_raw_rna = plot_pca(raw_counts_df, count_type="RNA", cond_vec=c(rep("T0", 3), rep("2hr", 3)), title_text="PCA Raw RNA Counts")
gg_pca_norm_rna = plot_pca(norm_counts_df, count_type="RNA", cond_vec=c(rep("T0", 3), rep("2hr", 3)), title_text="PCA Normalized RNA Counts")
gg_all_pca = grid.arrange(gg_pca_raw_rp, gg_pca_norm_rp, gg_pca_raw_rna, gg_pca_norm_rna, ncol=2)

# plot expression boxplots    
full_df = get_melted_merged_table(raw_counts_df, norm_counts_df, count_type_df)
ribo_rna_boxplot_raw = plot_boxplot(full_df, col_interest="raw_count", title="Raw Counts")
ribo_rna_boxplot_norm = plot_boxplot(full_df, col_interest="norm_count", title="Normalized Counts")
gg_all_boxplot = grid.arrange(ribo_rna_boxplot_raw, ribo_rna_boxplot_norm, ncol=1)

# plot correlations between sample types
gg_corr_raw = plot_corr(full_df, use_normalized=F, title="Raw Counts Correlation")
gg_corr_norm = plot_corr(full_df, use_normalized=T, title="Normalized Counts Correlation")


# plot ribodiff expression results
gg_volcano = plot_volcano(norm_counts_df, gene_names_df, ribodiff_df)


# write them out
outfile = paste(outdir, "qc_plots_pca.pdf", sep="")
ggsave(filename = outfile, plot = gg_all_pca, width = 10, height=10)

outfile = paste(outdir, "qc_plots_boxplot.pdf", sep="")
ggsave(filename = outfile, plot = gg_all_boxplot, width = 10, height=10)


outfile = paste(outdir, "qc_plots_RNA_RP_corr.pdf", sep="")
ggsave(filename = outfile, plot = gg_all_boxplot, width = 10, height=10)

