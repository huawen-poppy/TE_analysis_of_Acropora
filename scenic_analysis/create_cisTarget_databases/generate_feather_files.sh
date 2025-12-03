conda activate create_cistarget_databases

create_cistarget_databases_dir="/Acropora_sc_analysis/scenic_analysis/create_cisTarget_databases"

fasta_filename="/Acropora_sc_analysis/scenic_analysis/motif_analysis/filtered_renamed_ahem_promoter_up2000.fasta"
# previously used: promoter_region_used_gene_name_filter.fasta

original_species_fasta_filename="/Acropora_sc_analysis/scenic_analysis/motif_analysis/filtered_renamed_ahem_promoter_up2000.fasta"

motifs_dir="/Acropora_sc_analysis/scenic_analysis/create_cisTarget_databases/motifs_cb_format/motifs_cb_format"

motifs_list_filename="/Acropora_sc_analysis/scenic_analysis/create_cisTarget_databases/motifs_cb_format/motifs.lst"

db_prefix="ahem"

nbr_threads=60

python ${create_cistarget_databases_dir}/create_cistarget_motif_databases.py \
    -f "${fasta_filename}" \
    -F "${original_species_fasta_filename}" \
    -M "${motifs_dir}" \
    -m "${motifs_list_filename}" \
    -o "${db_prefix}" \
    -t "${nbr_threads}"

#cp CC7.motifs_vs_regions.rankings.feather test.CC7.motifs_vs_regions.rankings.feather

#db_prefix="test"

#output_dir="/ibex/scratch/projects/c2101/AipSC_analysis_old_genome/scenic_analysis/create_cisTarget_databases/"

#python ${create_cistarget_databases_dir}/create_cross_species_motifs_rankings_db.py \
#    -i "${db_prefix}" \
#    -o "${output_dir}"
