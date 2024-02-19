#/usr/bin/env python3
import numpy as np
from Bio import SeqIO
import sys
import io
import os
import pandas as pd
import subprocess

from numpy import dtype


def extract_seq_nor(seq_list_file, fasta_file, already_there_file, outfile):
    seq_list = [x.replace('\n', '').split("#")[0] for x in
                open(seq_list_file).readlines()]  # list file (one id per line)
    kept_seqs = []
    remove_seqs = []
    already_kept = [str(te.id).split("#")[0] for te in SeqIO.parse(already_there_file, "fasta")]

    for seq in SeqIO.parse(fasta_file, "fasta"):
        real_TE_name = str(seq.id).split("#")[0].replace("\n", "")
        if real_TE_name not in seq_list and real_TE_name not in already_kept:
            kept_seqs.append(seq)
        else:
            remove_seqs.append(seq)

    print("Number of kept sequences: " + str(len(kept_seqs)))
    SeqIO.write(kept_seqs, outfile, "fasta")


def extract_seqs(seq_list_file, fasta_file, outfile):
    seq_list = [x.replace('\n', '').split("#")[0] for x in open(seq_list_file).readlines()]
    kept_seqs = []
    remove_seqs = []

    for seq in SeqIO.parse(fasta_file, "fasta"):
        if str(seq.id).split("#")[0].replace("\n", "") in seq_list:
            kept_seqs.append(seq)
        else:
            remove_seqs.append(seq)

    print("Number of kept sequences: " + str(len(kept_seqs)))
    SeqIO.write(kept_seqs, outfile, "fasta")


def extract_seqs_from_list(seq_list_input, fasta_file, outfile):
    seq_list = [x.replace('\n', '').split("#")[0] for x in seq_list_input]
    kept_seqs = []
    remove_seqs = []

    for seq in SeqIO.parse(fasta_file, "fasta"):
        if str(seq.id).split("#")[0].replace("\n", "") in seq_list:
            kept_seqs.append(seq)
        else:
            remove_seqs.append(seq)

    print("Number of kept sequences: " + str(len(kept_seqs)))
    SeqIO.write(kept_seqs, outfile, "fasta")


def assign_families(spe_name, TE_lib, blast_80_80, outfile):
    spe_name
    TE_80_80 = {}
    final_tes = []
    tes_to_MI = []

    # to get the TE classified by the 80-80-80 rule
    lines_80_80 = open(blast_80_80, "r").readlines()
    for line in lines_80_80:
        hit_len = int(line.split("\t")[3])
        if hit_len > 80:
            query_name = line.split("\t")[0]
            family_name = line.split("\t")[1]
            TE_80_80[query_name.split("#")[0]] = family_name.split("#")[0]

    print("number of TE classified by 80-80-80: " + str(len(TE_80_80)))

    num_seq = 0
    for TE in SeqIO.parse(TE_lib, "fasta"):
        te_name = TE.id.split("#")[0]
        classification = TE.id.split("#")[1]
        if te_name not in TE_80_80.keys():
            tes_to_MI.append(TE)
            num_seq += 1
        else:
            original_family = TE_80_80[te_name]
            times = len([x for x in final_tes if original_family in x.id])
            if times > 0:
                original_family += "_" + str(times + 1)
            TE.id = original_family + "#" + classification
            TE.description = ""
            final_tes.append(TE)

    SeqIO.write(final_tes, outfile + "_808080.fa", "fasta")
    SeqIO.write(tes_to_MI, outfile + "_non808080.fa", "fasta")


def run_cross_match(fasta_file, lib_path, matrix, toolkit_path):
    output = subprocess.run(
        [toolkit_path + '/cross_match', '-bandwidth', '100', '-gap_ext', '-5', '-gap_init', '-25', '-minscore', '200', '-matrix', matrix, fasta_file, lib_path],
        stdout=subprocess.PIPE, text=True)

    lines = [line.replace(" ", ";").replace(";;", ";").replace(";;", ";").replace(";;", ";") for line in output.stdout.split("\n") if "#CLASS" in line and "Residues:" not in line]
    dataframe = pd.DataFrame(index=range(len(lines)), columns=["score", "div", "inser", "delet", "model", "qstart", "qend", "qleft", "target", "tstart", "tend", "tleft"])
    index = 0
    for line in lines:
        columns = [row for row in line.split(";") if row != '' and row != 'C']
        dataframe.iloc[index] = {"score": columns[0], "div": columns[1], "inser": columns[2], "delet": columns[3],
                            "model": columns[4], "qstart": columns[5], "qend": columns[6], "qleft": columns[7],
                            "target": columns[8], "tstart": columns[9], "tend": columns[10], "tleft": columns[11]}
        index += 1

    return dataframe


def synonymous(order):
    if order == "TIR":
        return "DNA"
    if order == "MITE":
        return "DNA"
    if order == "HELITRON":
        return "RC"
    if order == "CRYPTON":
        return "DNA"
    if order == "DIRS":
        return "LTR"


def order_match(cm_dataframe, tolerance, score):
    min_score = score * (1-tolerance)
    model_list = list(set(cm_dataframe.loc[:, 'model']))
    kept_models = []
    for model in model_list:
        matches = cm_dataframe.loc[cm_dataframe['model'] == model].reset_index()
        if len(model.split("#")[1].split("/")) > 1:
            order = model.split("#")[1].split("/")[1]
            for i in range(matches.shape[0]):
                if min_score < int(matches.loc[i, 'score']):
                    order_hit = matches.loc[i, 'target'].split("#")[1].split("/")[0].split(" ")[0]
                    if order == order_hit or synonymous(order) == order_hit:
                        kept_models.append(model)
                        break

    return kept_models


def rename_short(library, short_name, outfile):
    seq_list = []
    kept_seqs = []

    for seq in SeqIO.parse(library, "fasta"):
        if len(str(seq.seq)) > 100:
            kept_seqs.append(seq.id.split("#")[0])
            seqname = short_name + "_" + str(seq.id)
            seq.id = seqname
            seq.description = ""
            seq_list.append(seq)

    SeqIO.write(seq_list, outfile, "fasta")


def count_domains_by_order(profiles, order):
    right_doms = 0
    other_doms = 0
    if order == "LINE":
        right_doms = len(
            [x for x in profiles.split(",") if '_RT_' in x or '_EN_' in x or '_RNaseH_' in x or '_GAG_' in x])
        other_doms = len([x for x in profiles.split(",") if '_AP_' in x or '_INT_' in x or '_ENV_' in x or '_Tase_' in x
                          or '_HEL_' in x or '_RPA_' in x or '_REP_' in x or '_OTU_' in x or '_SET_' in x or '_Prp' in x
                          or '_ATPase_' in x or '_PhageINT_' in x])
    elif order == "LTR":
        right_doms = len([x for x in profiles.split(",") if
                        '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or '_RNaseH_' in x or '_ENV_' in x])
        other_doms = len([x for x in profiles.split(",") if '_EN_' in x or '_Tase_' in x or '_HEL_' in x or '_RPA_' in x
                          or '_REP_' in x or '_OTU_' in x or '_SET_' in x or '_Prp' in x or '_ATPase_' in x
                          or '_PhageINT_' in x])

    elif order == "DIRS":
        right_doms = len([x for x in profiles.split(",") if
                        '_GAG_' in x or '_RT_' in x or '_RNaseH_' in x or '_PhageINT_' in x])
        other_doms = len([x for x in profiles.split(",") if '_AP_' in x or '_INT_' in x or '_ENV_' in x or '_EN_' in x
                          or '_Tase_' in x or '_HEL_' in x or '_RPA_' in x or '_REP_' in x or '_OTU_' in x or '_SET_'
                          in x or '_Prp' in x or '_ATPase_' in x])

    elif order == "TIR":
        right_doms = len([x for x in profiles.split(",") if '_Tase_' in x])
        other_doms = len([x for x in profiles.split(",") if '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x
                          or '_RNaseH_' in x or '_ENV_' in x or '_EN_' in x or '_HEL_' in x or '_RPA_' in x or '_REP_'
                          in x or '_OTU_' in x or '_SET_' in x or '_Prp' in x or '_ATPase_' in x or '_PhageINT_' in x])

    elif order == "HELITRON":
        right_doms = len([x for x in profiles.split(",") if '_HEL_' in x or '_EN_' in x or '_RPA_' in x or '_REP_' in x
                          or '_OTU_' in x or '_SET_' in x])
        other_doms = len([x for x in profiles.split(",") if
                          '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or '_RNaseH_' in x or '_ENV_' in x
                          or '_Tase_' in x or '_Prp' in x or '_ATPase_' in x or '_PhageINT_' in x])

    elif order == "MAVERICK":
        right_doms = len([x for x in profiles.split(",") if '_Prp' in x or '_ATPase_' in x or '_INT_' in x or '_AP_' in x])
        other_doms = len([x for x in profiles.split(",") if
                          '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x or '_RNaseH_' in x or '_ENV_' in x
                          or '_EN_' in x or '_Tase_' in x or '_HEL_' in x or '_RPA_' in x or '_REP_' in x or '_OTU_' in x
                          or '_SET_' in x or '_PhageINT_' in x])
    elif order == "CRYPTON":
        right_doms = len([x for x in profiles.split(",") if '_PhageINT_' in x])
        other_doms = len([x for x in profiles.split(",") if '_GAG_' in x or '_AP_' in x or '_INT_' in x or '_RT_' in x
                          or '_RNaseH_' in x or '_ENV_' in x or '_EN_' in x or '_HEL_' in x or '_RPA_' in x or '_REP_'
                          in x or '_OTU_' in x or '_SET_' in x or '_Prp' in x or '_ATPase_' in x or '_Tase_' in x])

    return right_doms, other_doms


def stand_by_filtering(fasta_file, features_table_path, flf_path, species, toolkit_path):
    features_table = pd.read_csv(features_table_path, sep='\t')
    flf = pd.read_csv(flf_path, sep='\t')
    to_keep = []
    to_remove = []

    for te in SeqIO.parse(fasta_file, "fasta"):
        TE_name = str(te.id.split("#")[0])
        features_table_selected = features_table.loc[features_table["Seq_name"] == TE_name].reset_index()
        results = flf.loc[flf["TE"] == TE_name].reset_index()
        flf_i = results.at[0, "fullLgthFrags"]
        keep = True

        # case 1:
        if ";  )" in features_table_selected.at[0, "struct"]:  # no domains neither terminal repeats
            if flf_i < 3:
                print("Removed " + TE_name + " due to case 1")
                keep = False

        # case 2
        else:
            codings = features_table_selected.at[0, "coding"].replace("coding=(", " ").replace(")", "")
            codigs = codings.split(";")
            profiles = [cod for cod in codigs if "profiles:" in cod]
            order = features_table_selected.at[0, "order"]
            if len(profiles) > 0 and len(profiles[0].split(",")) >= 1:
                right_doms, other_doms = count_domains_by_order(profiles[0], order)
                bad_TR = False
                if order in ["LTR", "TRIM", "LARD", "PLE"] and (
                        "termTIR" in features_table_selected.at[0, "struct"] or "polyAtail" in
                        features_table_selected.at[0, "struct"]):
                    bad_TR = True
                elif order in ["TIR", "MITE", "MAVERICK"] and (
                        "termLTR" in features_table_selected.at[0, "struct"] or "polyAtail" in
                        features_table_selected.at[0, "struct"]):
                    bad_TR = True
                elif order in ["LINE", "SINE"] and (
                        "termTIR" in features_table_selected.at[0, "struct"] or "termLTR" in features_table_selected.at[
                    0, "struct"]):
                    bad_TR = True

                if (other_doms > 0 or bad_TR is True) and flf_i < 3:
                    print("Removed " + TE_name + " due to case 2")
                    keep = False

                # Case 3
                else:
                    # Step 0 extract the TE seq
                    SeqIO.write([te], "TE_seq.tmp", "fasta")

                    # Step 1 extract the ORF from the TE seq
                    output = subprocess.run(
                        ['getorf', '-sequence', "TE_seq.tmp", '--outseq', "TE.orfs", '-minsize', "400"],
                        stdout=subprocess.PIPE, text=True)

                    # Step2 blastp the ORF against the TE+Aid dataset
                    command = f"blastp -query TE.orfs -db "+toolkit_path+"/TE-Aid-master/db/RepeatPeps.lib -outfmt 6 | sort -k1,1 -k12,12nr | sort -u -k1,1 | sed 's/#/--/g'"
                    output = subprocess.check_output(command, shell=True).decode()
                    if len(output) > 0:
                        blastresult = pd.read_csv(io.StringIO(output), sep='\t', header=None)

                        if blastresult.shape[0] > 10 and flf_i < 3:
                            print("Removed " + TE_name + " due to case 3")
                            keep = False

        if keep:
            to_keep.append(te)
        else:
            to_remove.append(te)

    SeqIO.write(to_keep, species + "_standBy_kept.fa", "fasta")
    SeqIO.write(to_remove, species + "_standBy_v3.fa", "fasta")


def cut_pipeline(fasta_file, seqname, pos_start, pos_end, genome, te_aid_path):
    extract_seqs_from_list([seqname], fasta_file, "extracted_seq.fa")
    cut_seq = [te for te in SeqIO.parse("extracted_seq.fa", "fasta")]
    cut_seq = cut_seq[0]
    cut_seq.seq = cut_seq.seq[pos_start:pos_end]
    cut_seq.description = ""
    SeqIO.write([cut_seq], seqname.split("#")[0] + ".fa", "fasta")
    os.remove("extracted_seq.fa")
    output = subprocess.run(
        [te_aid_path + '/TE-Aid', '-q', seqname.split("#")[0] + ".fa", '-g', genome, '-o', seqname.split("#")[0] + ".fa_teAid"],
        stdout=subprocess.PIPE, text=True)

    # Create a file with the coverage:
    # run blastn command and read output into pandas dataframe
    command = f"blastn -query {seqname.split('#')[0] + '.fa'} -db {genome} -evalue 10e-8 -num_threads 1 -outfmt '6 qstart qend'"
    output = subprocess.check_output(command, shell=True).decode()

    if len(output) > 0:
        blast = pd.read_csv(io.StringIO(output), sep='\t', header=None)
        # create matrix of zeros with dimensions based on length of blast dataframe and cons_len variable
        coverage = np.zeros((len(blast), len(cut_seq.seq)), dtype=np.bool_)

        # iterate over rows of blast dataframe and set corresponding values in coverage matrix
        for i, row in blast.iterrows():
            start = int(row[0]) - 1
            end = int(row[1]) - 1
            if start <= end:
                coverage[i][start:end + 1] = 1
            else:
                coverage[i][end:start + 1] = 1

        # calculate column sums of coverage matrix and write to file
        coverage_sum = np.sum(coverage, axis=0, dtype=np.int32).T
        np.savetxt("TE_coverage.tab", coverage_sum, delimiter=",", fmt='%s')


if __name__ == '__main__':
    # variables for 70-70-70 filtering
    toolkit_path = os.path.dirname(os.path.realpath(__file__))
    matrix = toolkit_path + "/25p41g.matrix"
    lib_path = toolkit_path + "/your_own_Repbase_lib.fa"
    cross_match_path = toolkit_path + "/pathToCross_match/"
    te_aid_path = toolkit_path + "/pathToTE-Aid/"
    tolerance = 0.10
    score = 300

    option = sys.argv[1]

    if option == "extractseq_nor":
        if len(sys.argv) == 6:
            seq_list_file = sys.argv[2]  # sequences to not include
            fasta_file = sys.argv[3]  # fasta file with all the sequences
            already_there_file = sys.argv[4]  # fasta file with already included sequences
            outfile = sys.argv[5]  # Path/name of the output file
            extract_seq_nor(seq_list_file, fasta_file, already_there_file, outfile)
        else:
            print("ERROR in parameters, usage: python3 "+sys.argv[0]+" extractseq_nor seq_list_to_exclude.txt file_with_all_sequences.fasta file_with_already_included_seqs.fasta path/to/output.fasta")
            sys.exit(0)
    elif option == "extractseq":
        if len(sys.argv) == 5:
            seq_list_file = sys.argv[2]  # List of sequences to extract (without ">")
            fasta_file = sys.argv[3]  # Fasta file containing the sequences to extract
            outfile = sys.argv[4]  # Path/name of the output file
            extract_seqs(seq_list_file, fasta_file, outfile)
        else:
            print("ERROR in parameters, usage: python3 "+sys.argv[0]+" extractseq seq_list_to_extract.txt file_with_all_sequences.fasta path/to/output.fasta")
            sys.exit(0)
    elif option == "rename":
        if len(sys.argv) == 5:
            library = sys.argv[2]  # Fasta file containing the sequences to rename
            short_name = sys.argv[3]  # Short name of the species (i.e. Dmel for D. melanogaster)
            outfile = sys.argv[4]  # Path/name of the output file
            rename_short(library, short_name, outfile)
        else:
            print("ERROR in parameters, usage: python3 "+sys.argv[0]+" rename file_with_sequences_to_rename.fasta short_species_name path/to/output.fasta")
            sys.exit(0)
    elif option == "assign_families":
        if len(sys.argv) == 6:
            spe_name = sys.argv[2]  # Short name of the species (i.e. Dmel for D. melanogaster)
            TE_lib = sys.argv[3]  # Fasta file containing the sequences to assign families
            blast_80_80 = sys.argv[4]  # result of the BLASTn (tabular format, --outfmt 6)
            outfile = sys.argv[5]  # Path/name of the output file
            assign_families(spe_name, TE_lib, blast_80_80, outfile)
        else:
            print("ERROR in parameters, usage: python3 "+sys.argv[0]+" assign_families short_species_name file_with_sequences_tobe_assigned.fasta tabular_blastn_output path/to/output.fasta")
            sys.exit(0)
    elif option == "707070_filtering":
        if len(sys.argv) == 4:
            if not os.path.exists(lib_path):
                print("Installation error: Please set up the variable 'lib_path' correctly")
                sys.exit(0)
            if not os.path.exists(cross_match_path):
                print("Installation error: Please set up the variable 'cross_match_path' correctly")
                sys.exit(0)
            fasta_file = sys.argv[2]  # Fasta file containing the sequences to filter
            output_path = sys.argv[3]  # Path/name of the output file

            results = run_cross_match(fasta_file, lib_path, matrix, cross_match_path)
            kept_models = order_match(results, tolerance, score)
            extract_seqs_from_list(kept_models, fasta_file, output_path)
            print("Total number of model kept: " + str(len(kept_models)))
        else:
            print("ERROR in parameters, usage: python3 "+sys.argv[0]+" 707070_filtering file_with_sequences_tobe_filtered.fasta path/to/output.fasta")
            sys.exit(0)
    elif option == "stand_by_filtering":
        if len(sys.argv) == 6:
            fasta_file = sys.argv[2]  # Fasta file containing the sequences to filter
            features_table_path = sys.argv[3]  # Tabular file containing structural features from MCHelper (denovoLibTEs_PC.classif)
            flf_path = sys.argv[4]  # Tabular file containing number of full-length copies from MCHelper (fullLengthFrag.txt)
            species = sys.argv[5]  # Short name of the species (i.e. Dmel for D. melanogaster)
            toolkit_path = os.path.dirname(os.path.realpath(__file__))
            stand_by_filtering(fasta_file, features_table_path, flf_path, species, toolkit_path)
        else:
            print("ERROR in parameters, usage: python3 "+sys.argv[0]+" stand_by_filtering file_with_sequences_tobe_filtered.fasta /path/to/MCHelper_output/denovoLibTEs_PC.classif /path/to/MCHelper_output/fullLengthFrag.txt short_species_name path/to/output.fasta")
            sys.exit(0)
    elif option == "cut_pipeline":
        if len(sys.argv) == 7:
            if not os.path.exists(te_aid_path):
                print("Installation error: Please set up the variable 'te_aid_path' correctly")
                sys.exit(0)
            fasta_file = sys.argv[2]  # Fasta file containing the sequences to filter
            seqname = sys.argv[3]  # Tabular file containing structural features from MCHelper (denovoLibTEs_PC.classif)
            pos_start = int(sys.argv[4])  # Tabular file containing number of full-length copies from MCHelper (fullLengthFrag.txt)
            pos_end = int(sys.argv[5])  # Short name of the species (i.e. Dmel for D. melanogaster)
            genome = sys.argv[6]
            cut_pipeline(fasta_file, seqname, pos_start, pos_end, genome, te_aid_path)
        else:
            print("ERROR in parameters, usage: python3 "+sys.argv[0]+" cut_pipeline file_with_the_sequence_tobe_cut.fasta seq_id pos_start pos_end genome")
            sys.exit(0)
    else:
        print("The option "+option+" is not valid. Options: extractseq, extractseq_nor, assign_families, 707070_filtering, and stand_by_filtering")
        sys.exit(0)
