### this is to run the ribodiff analysis

### original

source /Users/natalie/Documents/projects/guido/.venv/bin/activate

expr_design=~/Documents/projects/guido/experimental_setup.txt
count_file=~/Documents/projects/guido/xpress_out/filtered_table_30.txt
outfile=~/Documents/projects/guido/xpress_out/filtered_30/result_30_85.txt
mm10_names=~/Documents/projects/guido/mm10_gene_trans_name.txt
python ~/Documents/projects/RiboDiff/scripts/TE.py -e ${expr_design} -c ${count_file} -o ${outfile} -p 1 -d 0 -s 1

#/cbio/grlab/share/software/R/R-3.1.1/bin/Rscript add_genenames.R mouse ${outfile}

### updated


expr_design=~/Documents/projects/guido/experimental_setup.txt
count_file=~/Documents/projects/guido/xpress_tables/cds_xpressFP_previousRNA_filtered_table_15.txt
outfile=~/Documents/projects/guido/xpress_tables/cds_xpressFP_previousRNA_filtered_table_15/result_15_85.txt
mm10_names=~/Documents/projects/guido/mm10_gene_trans_name.txt
python ~/Documents/projects/RiboDiff/scripts/TE.py -e ${expr_design} -c ${count_file} -o ${outfile} -p 1 -d 0 -s 1

#Library size:
#['Sample_T0_1_FP' 'Sample_T0_2_FP' 'Sample_T0_3_FP' 'Sample_2hr_2_FP'
# 'Sample_2hr_3_FP']
#[ 0.7  0.9  1.   1.4  1.5]
#['Sample_T0_1_RNA' 'Sample_T0_2_RNA' 'Sample_T0_3_RNA' 'Sample_2hr_1_RNA'
# 'Sample_2hr_2_RNA' 'Sample_2hr_3_RNA']
#[ 1.544  0.989  1.538  0.63   0.842  1.011]


expr_design=~/Documents/projects/guido/experimental_setup.txt
count_file=~/Documents/projects/guido/xpress_tables/cds_xpressFP_xpressRNA_filtered_table_15.txt
outfile=~/Documents/projects/guido/xpress_tables/cds_xpressFP_xpressRNA_filtered_table_15/result_15_85.txt
mm10_names=~/Documents/projects/guido/mm10_gene_trans_name.txt
python ~/Documents/projects/RiboDiff/scripts/TE.py -e ${expr_design} -c ${count_file} -o ${outfile} -p 1 -d 0 -s 1
#Library size:
#['Sample_T0_1_FP' 'Sample_T0_2_FP' 'Sample_T0_3_FP' 'Sample_2hr_2_FP'
# 'Sample_2hr_3_FP']
#[ 0.7  0.9  1.   1.4  1.5]
#['Sample_T0_1_RNA' 'Sample_T0_2_RNA' 'Sample_T0_3_RNA' 'Sample_2hr_1_RNA'
# 'Sample_2hr_2_RNA' 'Sample_2hr_3_RNA']
#[ 1.576  1.006  1.558  0.621  0.83   0.994]


expr_design=~/Documents/projects/guido/experimental_setup.txt
count_file=~/Documents/projects/guido/xpress_tables/transcript_xpressFP_usualRNA_filtered_table_15.txt
outfile=~/Documents/projects/guido/xpress_tables/transcript_xpressFP_usualRNA_filtered_table_15/result_15_85.txt
mm10_names=~/Documents/projects/guido/mm10_gene_trans_name.txt
python ~/Documents/projects/RiboDiff/scripts/TE.py -e ${expr_design} -c ${count_file} -o ${outfile} -p 1 -d 0 -s 1
#Library size:
#['Sample_T0_1_FP' 'Sample_T0_2_FP' 'Sample_T0_3_FP' 'Sample_2hr_2_FP'
# 'Sample_2hr_3_FP']
#[ 0.667  0.167  1.     1.25   1.5  ]
#['Sample_T0_1_RNA' 'Sample_T0_2_RNA' 'Sample_T0_3_RNA' 'Sample_2hr_1_RNA'
# 'Sample_2hr_2_RNA' 'Sample_2hr_3_RNA']
#[ 1.52   0.974  1.51   0.642  0.855  1.026]


expr_design=~/Documents/projects/guido/experimental_setup.txt
count_file=~/Documents/projects/guido/xpress_tables/transcript_xpressFP_xpressRNA_filtered_table_15.txt
outfile=~/Documents/projects/guido/xpress_tables/transcript_xpressFP_xpressRNA_filtered_table_15/result_15_85.txt
mm10_names=~/Documents/projects/guido/mm10_gene_trans_name.txt
python ~/Documents/projects/RiboDiff/scripts/TE.py -e ${expr_design} -c ${count_file} -o ${outfile} -p 1 -d 0 -s 1
#Library size:
#['Sample_T0_1_FP' 'Sample_T0_2_FP' 'Sample_T0_3_FP' 'Sample_2hr_2_FP'
# 'Sample_2hr_3_FP']
#[ 0.7  0.9  1.   1.4  1.5]
#['Sample_T0_1_RNA' 'Sample_T0_2_RNA' 'Sample_T0_3_RNA' 'Sample_2hr_1_RNA'
# 'Sample_2hr_2_RNA' 'Sample_2hr_3_RNA']
#[ 1.576  1.006  1.558  0.621  0.83   0.994]
