
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
