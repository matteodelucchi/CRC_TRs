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

# ##############################################################
# ### Download Protein Sequences

# ### Download Pathway Proteins:
# # name of the file with all protein sequences
# output_file_path_wnt = "./data/wnt_proteins_CRC_sp.fasta"
# output_file_path_nfkappab = "./data/nfkappab_proteins_CRC_sp.fasta"
# output_file_path_pi3k_akt_pten = "./data/pi3k_akt_pten_proteins_CRC_sp.fasta"
# output_file_path_ras_raf = "./data/ras_raf_proteins_CRC_sp.fasta"
# output_file_path_gsk3beta = "./data/gsk3beta_proteins_CRC_sp.fasta"

# # Url from which the sequences are requested:
# url_wnt = "https://www.uniprot.org/uniprot/?query=wnt%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"
# url_nfkappab = "https://www.uniprot.org/uniprot/?query=nf%20kappa%20b%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"
# url_pi3k_akt_pten = "https://www.uniprot.org/uniprot/?query=PI3K%20AKT%20PTEN%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"
# url_ras_raf = "https://www.uniprot.org/uniprot/?query=ras%20raf%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"
# url_gsk3beta = "https://www.uniprot.org/uniprot/?query=gsk3b%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"

# # Retrieve Protein Sequences
# retrieve_protein_seq(output_file_path_wnt, url_wnt)
# retrieve_protein_seq(output_file_path_nfkappab, url_nfkappab)
# retrieve_protein_seq(output_file_path_pi3k_akt_pten, url_pi3k_akt_pten)
# retrieve_protein_seq(output_file_path_ras_raf, url_ras_raf)
# retrieve_protein_seq(output_file_path_gsk3beta, url_gsk3beta)


# ### Download CRC favorable and unfavorable Proteins by gene name:
# gene_file_name_unfav = "./data/prognostic_colorectal_genes_unfavorable.tsv"
# gene_file_name_fav = "./data/prognostic_colorectal_genes_favorable.tsv"

# # name of the file with all protein sequences
# output_file_path_unfavorable = "./data/unfavorable_proteins_CRC_sp.fasta"
# output_file_path_favorable = "./data/favorable_proteins_CRC_sp.fasta"

# # Get Gene names
# unfavorable_genes = extract_gene_names_from_prognostics_list(
#     gene_file_name_unfav)
# favorable_genes = extract_gene_names_from_prognostics_list(
#     gene_file_name_fav)

# # Get Protein sequences
# retrieve_protein(output_file_path_unfavorable, unfavorable_genes)
# retrieve_protein(output_file_path_favorable, favorable_genes)

##############################################################
### Filter for Tandem Repeats
# Defining Paths and Parameters
sequences_file_unfav = "./data/unfavorable_proteins_CRC_sp.fasta"
sequences_file_fav = "./data/favorable_proteins_CRC_sp.fasta"
sequences_file_wnt = "./data/wnt_proteins_CRC_sp.fasta"
sequences_file_nfkappab = "./data/nfkappab_proteins_CRC_sp.fasta"
sequences_file_pi3k_akt_pten = "./data/pi3k_akt_pten_proteins_CRC_sp.fasta"
sequences_file_ras_raf = "./data/ras_raf_proteins_CRC_sp.fasta"
sequences_file_gsk3beta = "./data/gsk3beta_proteins_CRC_sp.fasta"

# Output paths
result_dir_unfav = "./results/TRs_unfavorable_proteins_CRC_sp_l1000"
result_dir_fav = "./results/TRs_favorable_proteins_CRC_sp_l1000"
result_dir_wnt = "./results/TRs_Wnt_proteins_CRC_sp_l1000"
result_dir_nfkappab = "./results/TRs_nfkappab_proteins_CRC_sp_l1000"
result_dir_pi3k_akt_pten = "./results/TRs_pi3k_akt_pten_proteins_CRC_sp_l1000"
result_dir_ras_raf = "./results/TRs_ras_raf_proteins_CRC_sp_l1000"
result_dir_gsk3beta = "./results/TRs_gsk3beta_proteins_CRC_sp_l1000"

# Thresholds for filtering
#pvalue_threshold = 0.05
#divergence_threshold = 0.1
#n_threshold = 2.5  # minimun repeat unit count
#l_threshold = 3  # maximum repeat unit length

# If duplicated proteins, print them:
print(return_duplicate_proteins(sequences_file=sequences_file_unfav))
print(return_duplicate_proteins(sequences_file=sequences_file_fav))
print(return_duplicate_proteins(sequences_file=sequences_file_wnt))
print(return_duplicate_proteins(sequences_file=sequences_file_nfkappab))
print(return_duplicate_proteins(sequences_file=sequences_file_pi3k_akt_pten))
print(return_duplicate_proteins(sequences_file=sequences_file_ras_raf))
print(return_duplicate_proteins(sequences_file=sequences_file_gsk3beta))

# Call TRAL
# find_protein_repeats(sequences_file=sequences_file_unfav,
#                      result_dir=result_dir_unfav, l_threshold=1000)
# find_protein_repeats(sequences_file=sequences_file_fav,
#                      result_dir=result_dir_fav, l_threshold=1000)
# find_protein_repeats(sequences_file=sequences_file_wnt,
#                      result_dir=result_dir_wnt, l_threshold=1000)
find_protein_repeats(sequences_file=sequences_file_nfkappab,
                     result_dir=result_dir_nfkappab, l_threshold=1000)
find_protein_repeats(sequences_file=sequences_file_pi3k_akt_pten,
                     result_dir=result_dir_pi3k_akt_pten, l_threshold=1000)
find_protein_repeats(sequences_file=sequences_file_ras_raf,
                     result_dir=result_dir_ras_raf, l_threshold=1000)
find_protein_repeats(sequences_file=sequences_file_gsk3beta,
                     result_dir=result_dir_gsk3beta, l_threshold=1000)

##############################################################
### Concatenate Multiple TR in one file
result_dir_fav = "./results/TRs_favorable_proteins_CRC_sp_l1000"
output_tsv_file_fav = "./results/TRs_favorable_proteins_CRC_sp_l1000.tsv"
result_dir_unfav = "./results/TRs_unfavorable_proteins_CRC_sp_l1000"
output_tsv_file_unfav = "./results/TRs_unfavorable_proteins_CRC_sp_l1000.tsv"
result_dir_wnt = "./results/TRs_Wnt_proteins_CRC_sp_l1000"
output_tsv_file_wnt = "./results/TRs_Wnt_proteins_CRC_sp_l1000.tsv"
result_dir_nfkappab = "./results/TRs_nfkappab_proteins_CRC_sp_l1000"
output_tsv_file_nfkappab = "./results/TRs_nfkappab_proteins_CRC_sp_l1000.tsv"
result_dir_pi3k_akt_pten = "./results/TRs_pi3k_akt_pten_proteins_CRC_sp_l1000"
output_tsv_file_pi3k_akt_pten = "./results/TRs_pi3k_akt_pten_proteins_CRC_sp_l1000.tsv"
result_dir_ras_raf = "./results/TRs_ras_raf_proteins_CRC_sp_l1000"
output_tsv_file_ras_raf = "./results/TRs_ras_raf_proteins_CRC_sp_l1000.tsv"
result_dir_gsk3beta = "./results/TRs_gsk3beta_proteins_CRC_sp_l1000"
output_tsv_file_gsk3beta = "./results/TRs_gsk3beta_proteins_CRC_sp_l1000.tsv"

concatenate_TR_results(result_dir_fav, out_file=output_tsv_file_fav)
concatenate_TR_results(result_dir_unfav, out_file=output_tsv_file_unfav)
concatenate_TR_results(result_dir_wnt, out_file=output_tsv_file_wnt)
concatenate_TR_results(result_dir_nfkappab, out_file=output_tsv_file_nfkappab)
concatenate_TR_results(result_dir_pi3k_akt_pten, out_file=output_tsv_file_pi3k_akt_pten)
concatenate_TR_results(result_dir_ras_raf, out_file=output_tsv_file_ras_raf)
concatenate_TR_results(result_dir_gsk3beta, out_file=output_tsv_file_gsk3beta)
