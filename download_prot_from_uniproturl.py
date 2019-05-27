#!/usr/bin/env python3
import urllib.request
import csv



def retrieve_protein_seq(output_file_path, url):
    """
    Retrieve for the given url as REST API request its proteins in fasta format.
    Append them in a file and store it as .fasta file.
    
    Args: 
        url (str): a url leading to a uniprot sequence list in fasta file format. Access by REST API.
        file_path (str): path/to/file/with/filename.fasta
        
    Output:
        .fasta file with protein sequences"""

    proteins = ""
    try:
        with urllib.request.urlopen(url) as response, open(output_file_path, 'w') as out_file:
            proteins += response.read().decode('utf-8')
            out_file.write(proteins)
        out_file.close() 
    except:
        print("Couldn't retrieve protein sequences.\n Connection issue or no protein sequences for this REST API request.")
            
    return print("Protein sequences saved in: {}".format(output_file_path))

if __name__ == "__main__":
    ##########################################################################
    ######### Set Arguments
    # name of the file with all protein sequences
    output_file_path_wnt = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/wnt_proteins_CRC_sp.fasta"
    # Url from which the sequences are requested:
    url = "https://www.uniprot.org/uniprot/?query=wnt%20AND%20reviewed:yes%20AND%20organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20proteome:up000005640&format=fasta&sort=score"

    ##########################################################################
    ######### Retrieve Protein Sequences
    retrieve_protein_seq(output_file_path_wnt, url)