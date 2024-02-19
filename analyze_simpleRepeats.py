#!/usr/env/bin python
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import sys

dicc_sr_pos = {}
trf_otuput = sys.argv[1]  # output from TRF program
fasta_input = sys.argv[2]  # Fasta file with sequences to test
with open(trf_otuput, 'r') as inFile:
    nbInSeq = 0
    for line in inFile:
        row = line.split(" ")
        if len(row) > 1 and "Sequence:" in row[0]:
            nbInSeq += 1
            seqName = row[1].replace('\n', '')
        if len(row) >= 14 and "Sequence:" not in row[0]:
            start = row[0]
            end = row[1]
            if seqName in dicc_sr_pos.keys():
                dicc_sr_pos[seqName].append([int(start), int(end)])
            else:
                dicc_sr_pos[seqName] = [[int(start), int(end)]]

# collapse all overlapping satellite matches
for sat in dicc_sr_pos.keys():
    pos_list = dicc_sr_pos[sat]
    new_list = [pos_list[0]]
    for i in range(1, len(pos_list)):
        # if there is an overlap
        overlap = False
        for j in range(len(new_list)):
            if (new_list[j][0] <= pos_list[i][0] <= new_list[j][1] or new_list[j][0] <= pos_list[i][1] <=
                new_list[j][1]) \
                    or (pos_list[i][0] <= new_list[j][0] <= pos_list[i][1] or pos_list[i][0] <= new_list[j][1] <=
                        pos_list[i][1]):
                new_start = new_list[j][0] if new_list[j][0] < pos_list[i][0] else pos_list[i][0]
                new_end = new_list[j][1] if new_list[j][1] > pos_list[i][1] else pos_list[i][1]
                new_list[j] = [new_start, new_end]
                overlap = True
        if not overlap:
            new_list.append([pos_list[i][0], pos_list[i][1]])
    dicc_sr_pos[sat] = new_list

dicc_sr = {}
for sat in dicc_sr_pos.keys():
    sum_bp = 0
    pos_list = dicc_sr_pos[sat]
    for i in range(len(pos_list)):
        sum_bp += pos_list[i][1] - pos_list[i][0] + 1
    dicc_sr[sat] = sum_bp

# Remove those TEs with SSR > given threshold or those that match with genes
final_seqs = []
list_len = []
for te in SeqIO.parse(fasta_input, "fasta"):
    if te.id in dicc_sr.keys():
        te_len = len(str(te.seq))
        list_len.append(te_len)
        lenSR = dicc_sr[te.id]
        if ((lenSR * 100) / te_len) >= 60:  # It hasn't less than the given threshold of SSR in its sequence
            print("remove: " + te.id + " " + str((lenSR * 100) / te_len))
        else:
            final_seqs.append(te)

plt.hist([x for x in list_len if x < 2000], bins=40, color='blue', edgecolor='black', alpha=0.7)

plt.title('Histogram')
plt.xlabel('length')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()
