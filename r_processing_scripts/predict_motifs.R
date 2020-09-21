
## run prediction of motif across all 5'UTRs
library(data.table)

ires_ref_file = "~/Documents/projects/guido/annotation/IRES_mouse_human.txt"
utr_infile = "~/Documents/projects/guido/xpress_out/filtered_15/res_table_utrs.txt"
motif_res_file = "~/Documents/projects/guido/xpress_out/filtered_15/res_table_motifs.tsv"

all_utrs = data.frame(fread(utr_infile))

gq_pat = "GGG[G]{0,7}[ACTG]{1,7}GGG[G]{0,10}[ACTG]{1,7}GGG[G]{0,10}[ACTG]{1,7}GGG[G]{0,10}"
use_fuzzy = F
gq_idx = grep(gq_pat, all_utrs$utr5)

# 	C yyyyyyyyyyyyyy[0,11,0]G MUST be at TSS
# http://regrna.mbc.nctu.edu.tw/php/entry.php?ID=R0011
top_pat = "^C[CT]{1,15}G"
top_idx = grep(top_pat, all_utrs$utr5)

pyr_pat = "[CT]{9,15}"
pyr_idx = grep(pyr_pat, all_utrs$utr5)


## get IRES genes
ires_in = data.frame(fread(ires_ref_file, header=F))
ires_in = rbind(ires_in, "Lars2")
ires_idx = which(all_utrs$gene_name %in% unlist(ires_in))

all_utrs$has_gq = FALSE
all_utrs$has_gq[gq_idx] = TRUE

all_utrs$has_top = FALSE
all_utrs$has_top[top_idx] = TRUE

all_utrs$has_ires = FALSE
all_utrs$has_ires[ires_idx] = TRUE

all_utrs$has_pyramidine = FALSE
all_utrs$has_pyramidine[pyr_idx] = TRUE

write.table(all_utrs, motif_res_file, sep="\t", quote=F, row.names=F)
