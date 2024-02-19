import sys
import re
import pandas as pd
from Bio import SeqIO

def assign_families(short_name, TE_lib, blast_80_80, outfile):
    TE_80_80 = {}
    final_tes = []
    tes_to_MI = []

    # to get the TE classified by the 80-80-80 rule
    lines_80_80 = open(blast_80_80, "r").readlines()
    for line in lines_80_80:
        hit_len = int(line.split("\t")[3])
        if hit_len > 80:
            query_name = line.split("\t")[0]
            if short_name + "_" not in query_name:
                query_name = short_name + "_" + query_name
            family_name = line.split("\t")[1]
            TE_80_80[query_name.split("#")[0]] = family_name.split("#")[0]

    print("number of TE classified by 80-80-80: " + str(len(TE_80_80)))
    # change TE names to no_fam_yet if they're not in 80-80-80 rule
    num_seq = 0
    for TE in SeqIO.parse(TE_lib, "fasta"):
        te_name = TE.id.split("#")[0]
        classification = TE.id.split("#")[1]
        if te_name not in TE_80_80.keys():
            final_tes.append(TE)
        else:
            original_family = TE_80_80[te_name]
            times = len([x for x in final_tes if original_family in x.id])
            if times > 0:
                original_family += "_" + str(times + 1)
            TE.id = short_name+"_"+original_family + "#" + classification
            TE.description = ""
            final_tes.append(TE)
            num_seq += 1

    SeqIO.write(final_tes, outfile, "fasta")


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


def choose_name(representative, list_seqs):
    if "_family-" not in representative:  # If the representative already has a known family, ise it
        return representative
    else:
        # first case: there is a sequence with already assign family
        for seq in list_seqs:
            if "_family-" not in seq:
                return seq

        # second case: there is no sequences with assign family:
        return representative


def clean_names(sequence, species_dicc_reversed):
    new_seq = sequence
    for spe in species_dicc_reversed.keys():
        new_seq = new_seq.replace(species_dicc_reversed[spe]+"_", "")
    return new_seq


def final_names(dicc_clusters, species_dicc_reversed):
    new_fam_num = {}
    outfile = open("changed_names.txt", "w")
    list_seqs = []
    for seq in dicc_clusters.keys():
        rep_name = choose_name(seq, dicc_clusters[seq])
        classification = rep_name.split("#")[1]
        if "_family-" not in rep_name:  # case 1: already known family
            clean_name = clean_names(rep_name.split("#")[0], species_dicc_reversed)

        else:  # case 2: new family !!!!
            superfamily = classification.split("/")[-1]
            if superfamily in new_fam_num.keys():
                new_fam_num[superfamily] = new_fam_num[superfamily] + 1
            else:
                new_fam_num[superfamily] = 1
            clean_name = "NF_" + superfamily + "_" + str(new_fam_num[superfamily])

        for ind_seq in dicc_clusters[seq]:
            outfile.write(ind_seq + "\t" + clean_name + "#" + classification + "\n")
    return "changed_names.txt"


def write_libs(TE_lib, outfile, short_name, filenames_path):
    final_TEs = []
    final_names_lines = [x for x in open(filenames_path, "r").readlines()]
    for te in SeqIO.parse(TE_lib, "fasta"):
        new_name = [x.split("\t")[1].replace("\n", "") for x in final_names_lines if x.split("\t")[0] == te.id]
        if "NF_" not in new_name[0].split("#")[0]:
            patron = re.compile(r'_\d+')
            family_name = re.sub(patron, '', new_name[0].split("#")[0])
            final_name = short_name + "_" + family_name
        else:
            final_name = short_name + "_" + new_name[0].split("#")[0]
        times = len([x for x in final_TEs if final_name in x.id])
        if times > 0:
            final_name += "_" + str(times + 1)
        te.id = final_name + "#" + new_name[0].split("#")[1]
        te.description = ""
        final_TEs.append(te)

    SeqIO.write(final_TEs, outfile, "fasta")


def family2superfamily(sequenceID):
    dicc_names_LINE = {"_R1": "CLASSI/LINE/R1", "_CR1-": "CLASSI/LINE/CR1", "_Jockey": "CLASSI/LINE/JOCKEY",
                  "_Loa-": "CLASSI/LINE/LOA", "_I-": "CLASSI/LINE/I", "_LOA-": "CLASSI/LINE/LOA",
                  "_Kiri-": "CLASSI/LINE/L2", "_bilbo": "CLASSI/LINE/LOA", "_RTE-": "CLASSI/LINE/RTE",
                  "_Cr1a": "CLASSI/LINE/CR1", "TART-A": "CLASSI/LINE/JOCKEY", "_baggins": "CLASSI/LINE/LOA",
                  "_worf": "CLASSI/LINE/CR1", "Outcast-1": "CLASSI/LINE/OUTCAST", "_HeT-A": "CLASSI/LINE/JOCKEY",
                  "_Rt1b": "CLASSI/LINE/R1", "_Rt1a": "CLASSI/LINE/R1", "con26_BS": "CLASSI/LINE/JOCKEY",
                  "_Juan": "CLASSI/LINE/JOCKEY", "con24_X-element": "CLASSI/LINE/JOCKEY"}

    dicc_names_LTR = {"_BEL": "CLASSI/LTR/BELPAO", "_Gypsy": "CLASSI/LTR/GYPSY", "_diver": "CLASSI/LTR/BELPAO",
                      "_invader": "CLASSI/LTR/GYPSY", "_gypsy": "CLASSI/LTR/GYPSY", "_rover": "CLASSI/LTR/BELPAO",
                      "_roo": "CLASSI/LTR/BELPAO", "UnFmclCluster039_RLX-incomp_COR": "CLASSI/LTR/GYPSY",
                      "_Max-element": "CLASSI/LTR/BELPAO", "_Copia-": "CLASSI/LTR/COPIA", "_Idefix": "CLASSI/LTR/GYPSY",
                      "con7_UnFmcl002_RXX-LARD": "CLASSI/LTR/LARD", "_Nobel": "CLASSI/LTR/BELPAO",
                      "con2_gtwin": "CLASSI/LTR/GYPSY", "_FBte0000588": "CLASSI/LTR/GYPSY"}

    dicc_names_TIR = {"_Tc1": "CLASSII/TIR/TC1MARINER", "_Mariner-": "CLASSII/TIR/TC1MARINER", "_hAT-": "CLASSII/TIR/HAT",
                      "_Transib": "CLASSII/TIR/TRANSIB", "_pogo": "CLASSII/TIR/TC1MARINER", "_hopper": "CLASSII/TIR/TRANSIB",
                      "FBte0000026": "CLASSII/TIR/TC1MARINER", "_Minos": "CLASSII/TIR/TC1MARINER", "_Bari": "CLASSII/TIR/TC1MARINER",
                      "Rehavkus": "CLASSII/TIR/MULE", "MARWOLEN": "CLASSII/TIR/TC1MARINER", "_Uhu_DAn": "CLASSII/TIR/TC1MARINER",
                      "Hoana7": "CLASSII/TIR/HAT", "_Homo4": "CLASSII/TIR/HAT"}

    family_name = sequenceID.split("#")[0]
    order = sequenceID.split("#")[1].split("/")[-1]

    if order == "LINE":
        newSuperfam = [dicc_names_LINE[x] for x in dicc_names_LINE if x in family_name]
        if len(newSuperfam) > 0:
            return newSuperfam[0]
        else:
            print("Not found: "+family_name)
            return "NONE"
    elif order == "LTR":
        newSuperfam = [dicc_names_LTR[x] for x in dicc_names_LTR if x in family_name]
        if len(newSuperfam) > 0:
            return newSuperfam[0]
        else:
            print("Not found: "+family_name)
            return "NONE"
    elif order == "TIR":
        newSuperfam = [dicc_names_TIR[x] for x in dicc_names_TIR if x in family_name]
        if len(newSuperfam) > 0:
            return newSuperfam[0]
        else:
            print("Not found: "+family_name)
            return "NONE"
    else:
        return "NONE"


def put_superfam(TE_lib, outfile):
    final_TEs = []
    for te in SeqIO.parse(TE_lib, "fasta"):
        classification = te.id.split("#")[1]
        if len(classification.split("/")) == 2:
            newSuperFam = family2superfamily(te.id)
            newClass = te.id.split("#")[1]
            if newSuperFam != "NONE":
                newClass = newSuperFam
            te.id = te.id.split("#")[0] + "#" + newClass
            te.description = ""

        final_TEs.append(te)

    SeqIO.write(final_TEs, outfile, "fasta")


def check(TE_lib, families_file):
    final_TEs = []
    families_list = open(families_file, "r").readlines()
    TEs = [te for te in SeqIO.parse(TE_lib, "fasta")]
    for fam in families_list:
        fam_i = fam.replace("\n", "")
        members = [te.id.split("#")[1] for te in TEs if "_"+fam_i+"_" in te.id]
        if len(set(members)) > 1:
            print(fam_i+" has more than one classification")
            print(set(members))
        else:
            print(fam_i + " it's ok")
            print(set(members))
    """for te in SeqIO.parse(TE_lib, "fasta"):
        fam_name = te.id.split("#")[0].split("_")[1:-1]
        fam_name = "_".join(fam_name)
        print(fam_name)"""
    """    if "NF_" in te.id and len(te.id.split("#")[0].split("_")) < 5:
        te.id = te.id.split("#")[0] + "_1#" + te.id.split("#")[1]
        te.description = ""
    final_TEs.append(te)
SeqIO.write(final_TEs, outfile, "fasta")"""


species_dicc_reversed = {'D.melanogaster': 'Dmel', 'D.mauritiana': 'Dmau', 'D.simulans': 'Dsim', 'D.sechellia': 'Dsec',
                        'D.erecta': 'Dere', 'D.orena': 'Dore', 'D.teissieri': 'Dte2', 'D.santomea': 'Dsan',
                        'D.yakuba': 'Dyak', 'D.ambigua': 'Damb', 'D.tristis': 'Dtrs', 'D.obscura': 'Dobs',
                        'D.subsilvestris': 'Dsus', 'D.bifasciata': 'Dbif', 'D.guanche': 'Dgua', 'D.subobscura': 'Dsuo',
                        'D.pseudoobscura': 'Dpso', 'D.persimilis': 'Dper', 'D.miranda': 'Dmir', 'D.lowei': 'Dlow',
                        'D.varians': 'Dvar', 'D.vallismaia': 'Dval', 'D.merina': 'Dmer', 'D.ercepeae': 'Derc',
                        'D.atripex': 'Datr', 'D.monieri': 'Dmon', 'D.ananassae': 'Dana', 'D.p_nigrens': 'Dpni',
                        'D.psepseudoananassae': 'Dpse', 'D.bipectinata': 'Dbip', 'D.parabipectinata': 'Dpbi',
                        'D.malmalerkotliana': 'Dmal', 'D.m_pallens': 'Dpal', 'D.basisetae': 'Dbas',
                        'D.crucigera': 'Dcru', 'D.engyochracea': 'Deng', 'D.grimshawi': 'Dgri', 'D.murphyi': 'Dmur',
                        'D.paucipunta': 'Dpau', 'D.pullipes': 'Dpul', 'D.sproati': 'Dspr', 'D.villosipedis': 'Dvil',
                        'D.albomicans': 'Dalb', 'D.kepulauana': 'Dkep', 'D.nasuta': 'Dnas', 'D.niveifrons': 'Dniv',
                        'D.pallidifrons': 'Dpald', 'D.sulfalbostrigata': 'Dsal', 'D.sulfbilimbata': 'Dsbi',
                        'D.sulfurigaster': 'Dsul', 'D.anceps': 'Danc', 'D.antonietae': 'Dant', 'D.arizonae': 'Dari',
                        'D.borborema': 'Dbor', 'D.hematofila': 'Dhem', 'D.leonis': 'Dleo', 'D.mayaguana': 'Dmay',
                        'D.meridiana': 'Dmera', 'D.meridionalis': 'Dmero', 'D.mettleri': 'Dmet', 'D.mojavensis': 'Dmoj',
                        'D.mulleri': 'Dmul', 'D.nigricuria': 'Dnig', 'D.pegasa': 'Dpeg', 'D.stalkeri': 'Dsta'}


operation = sys.argv[1]
if operation == "assign_fams":
    spe_name = sys.argv[2]
    TE_lib = sys.argv[3]
    blast_80_80 = sys.argv[4]
    outfile = sys.argv[5]
    short_name = species_dicc_reversed[spe_name]
    assign_families(short_name, TE_lib, blast_80_80, outfile)
elif operation == "create_names":
    cd_hit_out = sys.argv[2]
    dicc_clusters, sequences = clstr2dicc(cd_hit_out)
    filenames_path = final_names(dicc_clusters, species_dicc_reversed)
elif operation == "rename":
    spe_name = sys.argv[2]
    TE_lib = sys.argv[3]
    outfile = sys.argv[4]
    short_name = species_dicc_reversed[spe_name]
    write_libs(TE_lib, outfile, short_name, "changed_names.txt")
elif operation == "superfamily":
    TE_lib = sys.argv[2]
    outfile = sys.argv[3]
    put_superfam(TE_lib, outfile)
elif operation == "check":
    TE_lib = sys.argv[2]
    families_file = sys.argv[3]
    check(TE_lib, families_file)





