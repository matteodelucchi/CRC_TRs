""" 
Pipeline for the Analysis of Protein TRs Associated with CRC

Requirements:
    TRAL python package
"""
from download_prot_from_uniproturl import retrieve_protein_seq

from download_prot_from_gene import extract_gene_names_from_prognostics_list
from download_prot_from_gene import retrieve_protein

from TR_in_multiple_Protein import return_duplicate_proteins
from TR_in_multiple_Protein import find_protein_repeats

from concatenate_multiple_TRs import concatenate_TR_results

##############################################################
### Download Protein Sequences

### Download Wnt-Proteins:
# Wnt-proteins:
url_wnt = "https://www.uniprot.org/uniprot/?query=wnt%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"
# name of the file with all protein sequences
output_file_path_wnt = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/wnt_proteins_CRC_sp.fasta"
retrieve_protein_seq(output_file_path_wnt, url_wnt)

### Download CRC favorable and unfavorable Proteins by gene name:
gene_file_name_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/prognostic_colorectal_genes_unfavorable.tsv"
gene_file_name_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/prognostic_colorectal_genes_favorable.tsv"
# name of the file with all protein sequences
output_file_path_unfavorable = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/unfavorable_proteins_CRC_sp.fasta"
output_file_path_favorable = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/favorable_proteins_CRC_sp.fasta"
# Get Gene names
unfavorable_genes = extract_gene_names_from_prognostics_list(
    gene_file_name_unfav)
favorable_genes = extract_gene_names_from_prognostics_list(
    gene_file_name_fav)
# Get Protein sequences
retrieve_protein(output_file_path_unfavorable, unfavorable_genes)
retrieve_protein(output_file_path_favorable, favorable_genes)

##############################################################
### Filter for Tandem Repeats
# Defining Paths and Parameters
sequences_file_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/unfavorable_proteins_CRC_sp.fasta"
sequences_file_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/favorable_proteins_CRC_sp.fasta"
sequences_file_wnt = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/wnt_proteins_CRC_sp.fasta"

# Output paths
result_dir_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_unfavorable_proteins_CRC_sp"
result_dir_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_favorable_proteins_CRC_sp"
result_dir_wnt = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_Wnt_proteins_CRC_sp"

# Thresholds for filtering
#pvalue_threshold = 0.05
#divergence_threshold = 0.1
#n_threshold = 2.5  # minimun repeat unit count
#l_threshold = 3  # maximum repeat unit length

# If duplicated proteins, print them:
print(return_duplicate_proteins(sequences_file=sequences_file_unfav))
print(return_duplicate_proteins(sequences_file=sequences_file_fav))
print(return_duplicate_proteins(sequences_file=sequences_file_wnt))

# Call TRAL
find_protein_repeats(sequences_file=sequences_file_unfav,
                     result_dir=result_dir_unfav)
find_protein_repeats(sequences_file=sequences_file_fav,
                     result_dir=result_dir_fav)
find_protein_repeats(sequences_file=sequences_file_wnt,
                     result_dir=result_dir_wnt)

##############################################################
### Concatenate Multiple TR in one file
result_dir_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_favorable_proteins_CRC_sp"
output_tsv_file_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_favorable_proteins_CRC_sp.tsv"
result_dir_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_unfavorable_proteins_CRC_sp"
output_tsv_file_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_unfavorable_proteins_CRC_sp.tsv"
result_dir_wnt = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_Wnt_proteins_CRC_sp"
output_tsv_file_wnt = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_Wnt_proteins_CRC_sp.tsv"

concatenate_TR_results(result_dir_fav, out_file=output_tsv_file_fav)
concatenate_TR_results(result_dir_unfav, out_file=output_tsv_file_unfav)
concatenate_TR_results(result_dir_wnt, out_file=output_tsv_file_wnt)
