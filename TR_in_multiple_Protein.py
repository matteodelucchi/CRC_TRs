
import sys
import os
import pickle

from Bio import SeqIO

from pyfaidx import Fasta

import logging
import logging.config

from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.hmm import hmm

def return_duplicate_proteins(sequences_file):
    """
    Checks if there are duplicated protein sequences in the fasta file.
    If so, it outputs a prompt message.

    Args:
        sequence_file (str): whole path to file with >1 protein fasta-sequence

    Output:
        List of duplicate protein IDs
    """
    # Check for Duplicates
    seen = {}
    dupes = []

    for prot in SeqIO.parse(sequences_file, "fasta"): 
        protein = prot.id.split("|")[1]

        if protein not in seen:
            seen[protein] = 1
        else:
            if seen[protein] == 1:
                dupes.append(protein)
            seen[protein] += 1
        #if not dupes:
            #print("there are duplicates. Check them manually and evtl. remove them.")
    return dupes


def find_protein_repeats(sequences_file, result_dir, pvalue_threshold = 0.05, divergence_threshold = 0.1, n_threshold = 2.5, l_threshold = 3):
    """
    Finds tandem repeats (TRs) in the protein sequences provided in 'sequence_file'.
    Filters the TRs according to the thresholds.
    Saves the TRs in the 'result_dir' each individually with the protein identifier as filename.

    Args:
        sequence_file (str): whole path to file with >=1 protein fasta-sequence
        result_dir (str): path where the individual TRAL results are supposed to be stored.
        pvalue_threshold (int): p-value threshold for filtering
        divergence_threshold (int): divergence threshold for filtering
        n_threshold (int): minimun repeat unit count for filtering
        l_threshold (int): maximum repeat unit length for filtering
    
    Output:
        For each protein sequence, a single file with its TRs is produced and stored in 'result_dir'.
    """
    logging.config.fileConfig(config_file("logging.ini"))
    log = logging.getLogger('root')

    # define this as in the config
    CONFIG_GENERAL = configuration.Configuration.instance().config
    CONFIG = CONFIG_GENERAL["repeat_list"]
    score = CONFIG["model"]

    ##########################################################################
    # From .fasta to Class Sequence sequences
    proteins = Fasta(sequences_file)

    all_denovo_repeats = 0
    all_filtered_repeats = 0
    
    for pyfaidx in proteins:
        seq_name = pyfaidx.name.split("|")[1]

        # name is protein identifier
        seq = sequence.Sequence(seq=str(pyfaidx), name=seq_name)

        log.debug("Work on sequence {}".format(seq_name))
        ##########################################################################
        # Getting TRs

        denovo_list = seq.detect(denovo=True, realignment="proPIP")
        for TR in denovo_list.repeats:
            TR.calculate_pvalues()

        ##########################################################################
        # Filtering TRs

        # add number of denovo found repeats
        all_denovo_repeats += len(denovo_list.repeats)

        # filtering for pvalue
        denovo_list = denovo_list.filter(
            "pvalue",
            score,
            pvalue_threshold)

        # filtering for divergence
        denovo_list = denovo_list.filter(
            "divergence",
            score,
            divergence_threshold)

        # filtering for number of repeat units
        denovo_list = denovo_list.filter(
            "attribute",
            "n_effective",
            "min",
            n_threshold)

        # filtering for length of repeat units
        denovo_list = denovo_list.filter(
            "attribute",
            "l_effective",
            "max",
            l_threshold)

        ##########################################################################
        # Building HMM with hmmbuild

        # # De novo TRs were remastered with HMM
        denovo_hmm = [hmm.HMM.create(input_format='repeat', repeat=iTR)
                    for iTR in denovo_list.repeats]  # only possible with hmmbuild
        denovo_list_remastered = seq.detect(lHMM=denovo_hmm)

        ##########################################################################
        # Clustering

        # De novo TRs were clustered for overlap (common ancestry). Only best =
        # lowest p-Value and lowest divergence were retained.
        denovo_list_remastered = denovo_list.filter(
            "none_overlapping", ["common_ancestry"], {"pvalue": score, "divergence": score})

        ##########################################################################
        # Save Tandem Repeats
    
        # Create output directory if not already exists.
        try:
            if not os.path.isdir(result_dir):
                os.makedirs(result_dir)
        except:
            raise Exception(
                "Could not create path to result directory: {}".format(
                    os.path.dirname(result_dir)))
        
        # create filename
        output_pickle_file = os.path.join(result_dir, seq_name + ".pkl")
        output_tsv_file = os.path.join(result_dir, seq_name + ".tsv")
        # output_fasta_file = os.path.join(result_dir, seq_name + ".fasta")

        # save TR-file
        denovo_list_remastered.write(output_format="pickle", file=output_pickle_file)
        denovo_list_remastered.write(output_format="tsv", file=output_tsv_file)
        # denovo_list_remastered.write(output_format='fasta', file=output_fasta_file)



        all_filtered_repeats += len(denovo_list_remastered.repeats)
        print("\n***", seq_name, "***")
        print("denovo repeats:", len(denovo_list.repeats))
        print("repeats after filtering and clustering:",
            len(denovo_list_remastered.repeats))

        for i in range(len(denovo_list_remastered.repeats)):
            print(denovo_list_remastered.repeats[i])

    return print("\nThere where {} repeats found de novo.".format(all_denovo_repeats), "After filtering and clustering there where only {} repeats left.\n".format(
        all_filtered_repeats))

if __name__ == "__main__":
    ##########################################################################
    # Defining Paths and Parameters

    # AA reference
    # Input paths
    #working_dir = "/home/matteo/polybox/MSc_ACLS/master_thesis/data/proteom_reference/pickles"
    sequences_file_unfav = "./data/unfavorable_proteins_CRC_sp.fasta"
    sequences_file_fav = "./data/favorable_proteins_CRC_sp.fasta"

    # Output paths
    result_dir_unfav = "./results/TRs_unfavorable_proteins_CRC_sp"
    result_dir_fav = "./results/TRs_favorable_proteins_CRC_sp"

    # Thresholds for filtering
    #pvalue_threshold = 0.05
    #divergence_threshold = 0.1
    #n_threshold = 2.5  # minimun repeat unit count
    #l_threshold = 3  # maximum repeat unit length
    
    ##########################################################################
    # If duplicated proteins, print them:
    print(return_duplicate_proteins(sequences_file = sequences_file_unfav))
    print(return_duplicate_proteins(sequences_file = sequences_file_fav))

    ##########################################################################
    # Call TRAL
    find_protein_repeats(sequences_file = sequences_file_unfav, result_dir=result_dir_unfav)
    find_protein_repeats(sequences_file = sequences_file_fav, result_dir=result_dir_fav)
