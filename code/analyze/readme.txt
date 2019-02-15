##########
order of scripts:
##########

predict_main
=====================
inputs
- refs.txt
- strains.txt
- predict_args.txt
- alignments

outputs
- postions_tag.txt.gz
- introgressed_blocks_par/unk_tag.txt
- hmm_init_tag.txt
- hmm_tag.txt
- probs_tag.txt.gz

gene_predictions_chrm_main
=====================
inputs
- chromosome
- refs.txt
- strains.txt
- predict_args.txt
- introgressed_blocks_par_tag.txt
- genbank gene feature files for reference
- alignments

intermediate
- ../reference_chrX_genes.txt
- strain_genes_dic.pkl
- gene_strains_dic.pkl
- chrms_completed.pkl

outputs
- introgressed_blocks_chrX_par_tag_summary.txt
- genes_for_each_region_chrX_tag.txt
- regions_for_each_strain_chrX_tag.txt
- genes_for_each_strain_chrX_tag.txt
- strains_for_each_gene_chrX_tag.txt
- regions/rX.maf.gz
- regions/rX_annotated.txt.gz

annotate_regions_main
=====================
inputs
- chromosome
- refs.txt
- strains.txt
- predict_args.txt
- introgressed_blocks_par_tag.txt
- genbank gene feature files for reference
- ../reference_chrX_genes.txt
- alignments
- masked chromosome sequences

outputs
- site_summaries/strain_chrX_site_summary.txt.gz

mask_regions
=====================
inputs
- tag
- chromosome
- introgressed_blocks_par_tag_summary.txt
- site_summaries/strain_chrX_site_summary.txt.gz
- regions/rX.maf.gz

outputs
- regions/rX_masked.maf.gz

summary_plus_main
=====================
inputs
- tag
- introgressed_blocks_chrX_par_tag_summary.txt
- genes_for_each_region_chrX_tag.txt
- regions/rX_masked.maf.gz (or regions rX.maf.gz)

outputs
- introgressed_blocks_par_tag_summary_plus.txt

filter
=====================
inputs
- tag
- introgressed_blocks_par_tag_summary_plus.txt

outputs
- introgressed_blocks_filtered_par_tag_summary_plus.txt

gene_summary_main
=====================
inputs
- ../../data/S288c_verified_orfs.tsv
- introgressed_blocks_filtered_par_tag_summary_plus.txt

outputs
- gene_filtered_summary.txt


##########
plotting code:
##########

plot_evaluate
- stuff useful for choosing filtering thresholds etc

plot_summary
- general summaries across strains

plot_explore
- ... how is this different from plot_summary?

plot_genome
- plot introgression across individual chromosome with a bar for each strain

plot_genome_composite
- plot introgression across whole genome, aggregated across strains

plot_region
- plot one specific introgressed region (one strain)

plot_compare.R
- plot some things to compare genes I found to those Strope et al found
