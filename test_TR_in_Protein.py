from tral.paths import config_file, PACKAGE_DIRECTORY
from tral import configuration
from tral.sequence import sequence
from tral.hmm import hmm
from pyfaidx import Fasta

# dowload proteinsequence from: https://www.uniprot.org/uniprot/P98179

##########################################################################
# From .fasta to Class Sequence sequences
proteins = Fasta(
    "/home/matteo/polybox/MSc_ACLS/master_thesis/data/test_protein.fasta")

all_denovo_repeats = 0
all_filtered_repeats = 0

for pyfaidx in proteins:
    seq_name = pyfaidx.name.split("|")[1]
    # name is protein identifier
    seq = sequence.Sequence(seq=str(pyfaidx), name=seq_name)

# Saving this sequences as binary files:
# with open(sequence_pkl, 'wb') as f:
#    pickle.dump(seq, f)

##########################################################################
# Getting TRs
denovo_list = seq.detect(denovo=True)
# TODO: reinstall tral with all phylo files!
for TR in denovo_list.repeats:
    TR.calculate_pvalues()

print(TR)

# Saving this sequences as binary files:
# with open(TRs_pkl, 'wb') as f:
#    pickle.dump(denovo_list, f)

##########################################################################
# Filtering TRs
# define this as in the config
CONFIG_GENERAL = configuration.Configuration.instance().config
CONFIG = CONFIG_GENERAL["repeat_list"]
score = CONFIG["model"]

# Thresholds for filtering
pvalue_threshold = 0.05
divergence_threshold = 0.1
n_threshold = 2.5  # minimun repeat unit count
l_threshold = 3  # maximum repeat unit length


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

# De novo TRs were remastered with HMM
denovo_hmm = [hmm.HMM.create(input_format='repeat', repeat=iTR) 
for iTR in denovo_list.repeats]  # only possible with hmmbuild

denovo_list_remastered = seq.detect(lHMM=denovo_hmm)

##########################################################################
# Clustering

# De novo TRs were clustered for overlap (common ancestry). Only best =
# lowest p-Value and lowest divergence were retained.
denovo_list_remastered = denovo_list.filter(
    "none_overlapping", ["common_ancestry"], {
        "pvalue": score, "divergence": score})

##########################################################################
# Save Tandem Repeats
output_pickle_file = "/home/matteo/polybox/MSc_ACLS/master_thesis/data/test_denovo_list.pkl"
output_tsv_file = "/home/matteo/polybox/MSc_ACLS/master_thesis/data/test_denovo_list.tsv"

denovo_list_remastered.write(
    output_format="pickle", file=output_pickle_file)
    
denovo_list_remastered.write(output_format="tsv", file=output_tsv_file)

# function to save as fasta has to be integrated
all_filtered_repeats += len(denovo_list_remastered.repeats) # add number of clustered repeats
print("\n***", seq_name, "***")
print("denovo repeats:", len(denovo_list.repeats))
print("repeats after filtering and clustering:", len(denovo_list_remastered.repeats))

for i in range(len(denovo_list_remastered.repeats)):
    print(denovo_list_remastered.repeats[i])

print("\nThere where {} repeats found de novo.".format(all_denovo_repeats))
print("After filtering and clustering there where only {} repeats left.\n".format(all_filtered_repeats))
