import pandas as pd
from Bio import SeqIO
import subprocess
import os
import numpy as np
import matplotlib.pyplot as plt
import sys

# NTE50
genome_size = int(sys.argv[1])  # Genome size
perc_TEs = float(sys.argv[2].replace(",", "."))  # TE proportion in the genome (in %)
per_base_tes = genome_size * (perc_TEs/100)
cores = 32

ref_annot = sys.argv[3]  # the RepeatMasker annotation

annot_gff = pd.read_table(ref_annot, sep='\t',
                         names=['score', 'divergence', 'deletion', 'insertion', 'query', 'start', 'end', 'length',
                                'sense', 'family', 'classification', 'start_TE', 'end_TE', 'TE_left', 'ID',
                                'num_assembled', '%_ref'])
annot_lengths = []
ref_len = {}

for x in range(annot_gff.shape[0]):
    annot_lengths.append(int(annot_gff.loc[x, 'length']))

annot_lengths_sorted = sorted(annot_lengths, reverse=True)
perc_cov_ref = []
acum_sum_bp_ref = 0
acum_sum_bp = 0
medium_pos = -1
for i in range(len(annot_lengths_sorted)):
    acum_sum_bp += annot_lengths_sorted[i]
    perc_cov_ref.append((acum_sum_bp*100)/genome_size)
    if acum_sum_bp >= (per_base_tes / 2) and medium_pos == -1:
        medium_pos = i
        acum_sum_bp_ref = acum_sum_bp

print(str(medium_pos) + "\t" + str(annot_lengths_sorted[medium_pos]))
annot_lengths_ref_sorted = annot_lengths_sorted