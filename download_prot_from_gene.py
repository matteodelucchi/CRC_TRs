#!/usr/bin/env python3
import urllib.request
import csv



def retrieve_protein(output_file_path, genes):
    """
    Access Swissprot knowledgebase for Homo Sapiens.  Retrieve for each provided gene its proteins in fasta format.
    Append them in a file and store it as .fasta file.
    
    Args: 
        genes (list): defines for which genes the sequences should be taken
        file_path (str): path/to/file/with/filename.fasta
        
    Output:
        .fasta file with all protein sequences of the provided genes"""

    proteins = ""
    for gene in genes:
        """ The first URL retrieves all human proteins from swissprot of the provided genes.
            The second URL querys for the human reference proteome where reviewed and unreviewed (!) entries are collected. 
                Proteomes: Some proteomes have been (manually and algorithmically) selected as reference proteomes. 
                They cover well-studied model organisms and other organisms of interest for biomedical research and phylogeny.
        """
        #url = "https://www.uniprot.org/uniprot/?query=gene:{}&format=fasta&fil=organism:%22Homo%20sapiens%20(Human)%20[9606]%22%20AND%20reviewed:yes".format(gene)
        url = "https://www.uniprot.org/uniprot/?query=organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+proteome%3Aup000005640+AND+gene%3A{}&format=fasta".format(gene)

        print("Retrieving proteins for gene: {}".format(gene))
        try:
            with urllib.request.urlopen(url) as response, open(output_file_path, 'w') as out_file:
                proteins += response.read().decode('utf-8')
                out_file.write(proteins)
            out_file.close() 
        except TypeError: 
            print("No valid gene name provided. Can't iterate.")
        except:
            print("Couldn't retrieve proteins of gene: {} \n Connection issue or no proteins in Swissprot for this gene.".format(gene))
            pass
    return print("Protein sequences saved in: {}".format(output_file_path))


def extract_gene_names_from_prognostics_list(gene_file_name):
    try:
        csvfile = open(gene_file_name, 'rt')
        csvReader = csv.reader(csvfile, delimiter="\t")
        
        genes_names = list()
        for row in csvReader:
            genes_names.append(row[0])

        del genes_names[0] #remove header: "Gene"
        return genes_names

    except IOError:
        print("File not found")
        return 0


if __name__ == "__main__":
    ##########################################################################
    ######### Set Paths
    # get all genes of which the protein sequences is requested either from here as .tsv:
    # https://www.proteinatlas.org/humanproteome/pathology/colorectal+cancer#colorectal%20canceradditionalinfo
    # or enter them manually:
    #unfavorable_genes = ['POFUT2',  'LRCH4', 'CDX2'] 
    gene_file_name_unfav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/prognostic_colorectal_genes_unfavorable.tsv"
    gene_file_name_fav = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/prognostic_colorectal_genes_favorable.tsv"
    # name of the file with all protein sequences
    output_file_path_unfavorable = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/unfavorable_proteins_CRC.fasta"
    output_file_path_favorable = "/home/matteo/polybox/MSc_ACLS/master_thesis/CRC_TRs/data/favorable_proteins_CRC.fasta"

    ##########################################################################
    ######### Get Gene names
    unfavorable_genes = extract_gene_names_from_prognostics_list(gene_file_name_unfav)
    favorable_genes = extract_gene_names_from_prognostics_list(gene_file_name_fav)

    ##########################################################################
    ######### Get Protein sequences
    retrieve_protein(output_file_path_unfavorable, unfavorable_genes)
    retrieve_protein(output_file_path_favorable, favorable_genes)
