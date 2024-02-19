from Bio import SeqIO
import subprocess
import tempfile
import shutil
import os
import sys


def clstr2dicc(clstr_file):
    fileopen = open(clstr_file, "r").readlines()
    dicc_clusters = {}
    cluster_name = ""
    sequences = []
    representative = ""
    for line in fileopen:
        if line[0] == ">":  # it's a new cluster
            # special case of the first cluster
            if line.replace("\n", "") != ">Cluster 0":
                dicc_clusters[representative] = members
            cluster_name = line.replace("\n", "")
            members = []
        else:  # it's a member of the cluster
            member_name = line.split("\t")[1].split(" ")[1].replace(">", "").replace("...", "")
            members.append(member_name)
            sequences.append(member_name)
            if line.split("\t")[1].split(" ")[2].replace("\n", "") == "*":
                representative = member_name

    # add the last cluster
    dicc_clusters[representative] = members
    return dicc_clusters, sequences


def consensus_creation(dicc_clst, lib_path, out_path):
    consensus = []
    for rep in dicc_clst.keys():
        members = [te for te in SeqIO.parse(lib_path, "fasta") if te.id in dicc_clst[rep]]
        print("Doing: " + rep + ", Cluster len: " + str(len(members)))
        if len(members) > 1:
            SeqIO.write(members, "cluster.fa", "fasta")
            # Paso 1: Realizar el Multiple Sequence Alignment (MSA) con MAFFT
            mafft_output = tempfile.NamedTemporaryFile(delete=False)
            mafft_cmd = ['mafft', '--auto', '--quiet', '--thread', '32', 'cluster.fa']
            subprocess.run(mafft_cmd, stdout=mafft_output)
            shutil.copy(mafft_output.name, "cluster_MSA.fa")

            # Paso 2: Generar el archivo de consenso usando cons de EMBOSS
            trim_output = subprocess.run(
                ['trimal', '-in', 'cluster_MSA.fa', '-gt', '0.6', '-fasta', '-out',
                 'cluster_MSA.fa.trim'],
                stdout=subprocess.PIPE, stderr=None, text=True, timeout=1500)

            consensus_output = tempfile.NamedTemporaryFile(delete=False)
            cons_cmd = ['cons', '-sequence', "cluster_MSA.fa.trim", '-outseq', "consensus.fa"]
            subprocess.run(cons_cmd, stdout=consensus_output)

            cons_cluster = [x for x in SeqIO.parse("consensus.fa", "fasta")]
            cons_cluster[0].id = rep
            cons_cluster[0].description = ""

            consensus.append(cons_cluster[0])

            # Eliminar archivos temporales
            os.remove(mafft_output.name)
            os.remove(consensus_output.name)
            os.remove("cluster.fa")
            os.remove("cluster_MSA.fa.trim")
            os.remove("cluster_MSA.fa")
            os.remove("consensus.fa")
        else:
            consensus.append(members[0])

    SeqIO.write(consensus, out_path, "fasta")


# Example usage:
clstrFile = sys.argv[1]  # output from CD-HIT
lib_path = sys.argv[2]  # fasta file containing the TE sequences
out_path = sys.argv[3]  # path to the output file

dicc_clst, sequences = clstr2dicc(clstrFile)
consensus_creation(dicc_clst, lib_path, out_path)