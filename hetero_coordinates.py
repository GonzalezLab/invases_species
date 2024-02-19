import sys
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import multiprocessing
import time
import subprocess

def RM2pandas_lite(RM_curated):
    RMdf = pd.read_csv(RM_curated, delimiter="\t", header=None)

    return RMdf[[4, 5, 6]]


def get_coordinates_parallel(RM_dataframe, genome_file, windows_size, cores):
    start_time = time.time()
    dicc_abundancy = {}
    dicc_sizes = {}

    for seq in SeqIO.parse(genome_file, "fasta"):
        len_chr = len(seq.seq)
        dicc_sizes[seq.id] = len_chr
        if len_chr > 10 * windows_size:
            dicc_abundancy[seq.id] = [0] * len_chr
    end_time = time.time()
    print("Dics created and initialized... [" + str(end_time - start_time) + " seconds]")

    start_time = time.time()
    n = RM_dataframe.shape[0]
    seqs_per_procs = int(n / cores)
    remain = n % cores
    ini_per_thread = []
    end_per_thread = []
    for p in range(cores):
        if p < remain:
            init = p * (seqs_per_procs + 1)
            end = n if init + seqs_per_procs + 1 > n else init + seqs_per_procs + 1
        else:
            init = p * seqs_per_procs + remain
            end = n if init + seqs_per_procs > n else init + seqs_per_procs
        ini_per_thread.append(init)
        end_per_thread.append(end)

    pool = multiprocessing.Pool(processes=cores)
    localresults = [pool.apply_async(get_coordinates, args=[RM_dataframe.loc[ini_per_thread[x]:end_per_thread[x]-1, :].reset_index(),
                                                            dicc_abundancy]) for x in range(cores)]

    local_dfs = [p.get() for p in localresults]
    dicc_abundancy = local_dfs[0]
    local_dfs[0] = None
    for i in range(1, len(local_dfs)):
        print("Process " + str(i) + " done!!")
        for key in local_dfs[i].keys():
            local_dicc = local_dfs[i]
            dicc_abundancy[key] = [x + y for x, y in zip(dicc_abundancy[key], local_dicc[key])]
        local_dfs[i] = None

    pool.close()
    end_time = time.time()
    print("All processed finished... [" + str(end_time - start_time) + " seconds]")

    start_time = time.time()
    for chr_name, cover_i in dicc_abundancy.items():
        abundancy = [sum(cover_i[i:i+windows_size])/windows_size for i in range(0, len(cover_i) - windows_size, windows_size)]
        dicc_abundancy[chr_name] = abundancy

    end_time = time.time()
    print("Abundancy summed and saved... [" + str(end_time - start_time) + " seconds]")

    return dicc_abundancy


def get_coordinates(RM_dataframe, dicc_abundancy):
    for i in range(RM_dataframe.shape[0]):
        chr_name = RM_dataframe.at[i, 4]
        if chr_name in dicc_abundancy:
            start, end = RM_dataframe.at[i, 5], RM_dataframe.at[i, 6]
            dicc_abundancy[chr_name][start:end] = [x + 1 for x in dicc_abundancy[chr_name][start:end]]
    return dicc_abundancy


def get_coordinates2(RM_dataframe, genome_file, windows_size):
    start_time = time.time()
    dicc_abundancy = {}
    dicc_sizes = {}

    for seq in SeqIO.parse(genome_file, "fasta"):
        len_chr = len(seq.seq)
        dicc_sizes[seq.id] = len_chr
        dicc_abundancy[seq.id] = [0] * len_chr
    end_time = time.time()
    print("Dics created and initialized... [" + str(end_time - start_time) + " seconds]")

    start_time = time.time()
    for i in range(RM_dataframe.shape[0]):
        chr_name = RM_dataframe.at[i, 4]
        if chr_name in dicc_abundancy:
            start, end = RM_dataframe.at[i, 5], RM_dataframe.at[i, 6]
            dicc_abundancy[chr_name][start:end] = [x + 1 for x in dicc_abundancy[chr_name][start:end]]

    for chr_name, cover_i in dicc_abundancy.items():
        abundancy = [sum(cover_i[i:i+windows_size])/windows_size for i in range(0, len(cover_i) - windows_size, windows_size)]
        dicc_abundancy[chr_name] = abundancy
    end_time = time.time()
    print("All processed finished... [" + str(end_time - start_time) + " seconds]")
    print(dicc_abundancy)
    print(sorted(dicc_sizes.values()))
    for chr_name, cover_i in dicc_abundancy.items():
        plt.plot([x for x in range(len(cover_i))], cover_i, label=''+chr_name)
        plt.show()


def plot(dicc_abundancy, genome_file, windows_size):
    start_time = time.time()
    map_dict = {}
    output = subprocess.run(
        ['minimap2', '-x', 'asm5', '--secondary=no', genome_file, 'D.melanogaster_heterocromatin.fa'],
        stdout=subprocess.PIPE, text=True)
    lines = [line for line in output.stdout.split("\n") if "[M::" not in line and line != '']
    end_time = time.time()
    print("Minimap2 finished... [" + str(end_time - start_time) + " seconds]")
    #lines = [line for line in open("melanogaster_map.paf", "r").readlines()]
    for line in lines:
        chr_name = line.split("\t")[5]
        start = line.split("\t")[7]
        end = line.split("\t")[8]
        if chr_name in map_dict.keys():
            map_dict[chr_name].append((start, end))
        else:
            map_dict[chr_name] = [(start, end)]

    for chr_name, cover_i in dicc_abundancy.items():
        fig, ax = plt.subplots()

        outputfile = open("windows_"+chr_name+".txt", "w")
        outputfile.write("\n".join([str(x) for x in cover_i]))
        # TE abundance
        ax.plot([x for x in range(len(cover_i))], cover_i, label=chr_name, color='black')
        ax.fill_between([x for x in range(len(cover_i))], cover_i, color='lightblue', alpha=0.5)

        if chr_name in map_dict.keys():
            for start, end in map_dict[chr_name]:
                # Pintar recuadros grises al comienzo y al final
                ax.axvspan(round(int(start)/windows_size), round(int(end)/windows_size), color='lightgray', alpha=0.5)

        # Añadir etiquetas y leyenda
        ax.set_xlabel('Chromosome length')
        ax.set_xticks(range(0, len(cover_i), 20))
        ax.set_ylim(0, 2.5)
        ax.set_ylabel('TE abundance')
        ax.set_title(chr_name)
        ax.legend()

        # Mostrar el gráfico
        plt.show()


def get_hetero(genome_file):
    hetero = []
    for chro in SeqIO.parse(genome_file, "fasta"):
        if "Chr2L_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[0:530000], id=chro.id+"_hetero_5p", description="")
            newseq2 = SeqIO.SeqRecord(chro.seq[18870000:], id=chro.id + "_hetero_3p", description="")
            hetero.append(newseq1)
            hetero.append(newseq2)
        elif "Chr2R_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[0:5982495], id=chro.id+"_hetero_5p", description="")
            newseq2 = SeqIO.SeqRecord(chro.seq[24972477:], id=chro.id + "_hetero_3p", description="")
            hetero.append(newseq1)
            hetero.append(newseq2)
        elif "Chr3L_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[0:750000], id=chro.id+"_hetero_5p", description="")
            newseq2 = SeqIO.SeqRecord(chro.seq[19026900:], id=chro.id + "_hetero_3p", description="")
            hetero.append(newseq1)
            hetero.append(newseq2)
        elif "Chr_3R_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[0:6754278], id=chro.id+"_hetero_5p", description="")
            newseq2 = SeqIO.SeqRecord(chro.seq[31614278:], id=chro.id + "_hetero_3p", description="")
            hetero.append(newseq1)
            hetero.append(newseq2)
        elif "ChrX_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[0:1325967], id=chro.id+"_hetero_5p", description="")
            newseq2 = SeqIO.SeqRecord(chro.seq[21338973:], id=chro.id + "_hetero_3p", description="")
            hetero.append(newseq1)
            hetero.append(newseq2)

    SeqIO.write(hetero, "D.melanogaster_heterocromatin.fa", "fasta")


def get_eucrom(genome_file):
    hetero = []
    for chro in SeqIO.parse(genome_file, "fasta"):
        if "Chr2L_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[530000:18870000], id=chro.id+"_euchromatin", description="")
            hetero.append(newseq1)
        elif "Chr2R_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[5982495:24972477], id=chro.id+"_euchromatin", description="")
            hetero.append(newseq1)
        elif "Chr3L_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[750000:19026900], id=chro.id+"_euchromatin", description="")
            hetero.append(newseq1)
        elif "Chr_3R_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[6754278:31614278], id=chro.id+"_euchromatin", description="")
            hetero.append(newseq1)
        elif "ChrX_" in chro.id:
            newseq1 = SeqIO.SeqRecord(chro.seq[1325967:21338973], id=chro.id+"_euchromatin", description="")
            hetero.append(newseq1)

    SeqIO.write(hetero, "D.melanogaster_euchromatin.fa", "fasta")



if __name__ == '__main__':
    RM_file = sys.argv[1]
    genome_file = sys.argv[2]

    windows_size = 50000
    cores = 14
    # get_hetero(genome_file)  # only for D. melanogaster
    # get_eucrom(genome_file)  # only for D. melanogaster
    start_time = time.time()
    RM_dataframe = RM2pandas_lite(RM_file)
    end_time = time.time()
    print("RM_dataframe created [" + str(end_time - start_time) + " seconds]")
    dicc_abundancy = get_coordinates_parallel(RM_dataframe, genome_file, windows_size, cores)
    pd_df_r = pd.DataFrame(columns=["Chr", "Value"])
    for chrm in dicc_abundancy.keys():
        for val in dicc_abundancy[chrm]:
            pd_df_r = pd.concat([pd_df_r, pd.DataFrame({"Chr": chrm, "Value": [val]})], ignore_index=True)
    pd_df_r.to_csv("abundancy_windows.csv", sep=",")
    plot(dicc_abundancy, genome_file, windows_size)
