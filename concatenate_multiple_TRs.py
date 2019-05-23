#!/usr/bin/env python3
import os

def concatenate_TR_results(result_dir, out_file):
    """
    Creates one large file with TRs and their protein ID.
    Takes the results from 'TR_in_multiple_Protein.py' as input.
    
    Inputs:
        result_dir: The directory with files containing TRs for each protein. The filename should be the protein identifier.
        out_file: The path and filename of the file with all results.

    Output:
        All TRs are stored in one single .tsv and .pickle file labeled with their protein identifier facilitating further analysis.
    """
    # list all .tsv files in results_dir
    filenames = []
    for file in os.listdir(result_dir):
        if file.endswith(".tsv"):
            filenames.append(file)

    # open each file, write line by line to outfile
    with open(out_file, 'w') as outfile:
        for fname in filenames:
            with open(os.path.join(result_dir, fname)) as infile:
                # open each input-file
                for line in infile:
                    if str('begin') in line:
                        # if the line is a header, then add 'ID'
                        newline = 'ID\t'+line
                    else:
                        # if the line is a TR, add the protein-id (which is the filename) at the start
                        newline = str(fname.split('.')[0])+'\t'+line
                    outfile.write(newline)

    # remove all lines which are headings but not the first line
    with open(out_file, 'r') as outfile:
        lines = outfile.readlines()
    
    TRs = []
    TRs.append(lines[0])  # set the first line as heading
    for line in lines[1:]:
        if not str('begin') in line:
            # add only lines which are not headings
            TRs.append(line) 
    
    with open(out_file, 'w') as outfile:
        for line in TRs:
            outfile.write(line)

if __name__ == "__main__":
    result_dir_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_favorable_proteins_CRC"
    output_tsv_file_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_favorable_proteins_CRC.tsv"
    #output_pickle_file = "/home/matteo/polybox/MSc_ACLS/master_thesis/results/TRs_unfavorable_proteins_CRC.pkl"
    result_dir_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_unfavorable_proteins_CRC"
    output_tsv_file_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/results/TRs_unfavorable_proteins_CRC.tsv"

    concatenate_TR_results(result_dir_fav, out_file=output_tsv_file_fav)
    concatenate_TR_results(result_dir_unfav, out_file=output_tsv_file_unfav)