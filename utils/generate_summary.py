#! /usr/bin/env python3


import sys
import glob
import os
import re


result_dir = sys.argv[1]
fasta_files = glob.glob("{}/60_fasta/*.fasta".format(result_dir))

fasta_files = [fasta_file for fasta_file in fasta_files if "all" not in os.path.basename(fasta_file)]


transcript_ids = {}
transcript_ids_cooksCutoff_FALSE = {}


def add_pfam_id(transcript_ids, transcript_id, pfam_id):

    if transcript_id not in transcript_ids:

        transcript_ids[transcript_id] = [pfam_id]
    
    else:

        transcript_ids[transcript_id].append(pfam_id)

    return transcript_ids


def get_expressions(transcript_id):

    with open("{}/30_count_matrix/RSEM.isoform.TMM.EXPR.matrix".format(result_dir)) as f:

        for line in f:

            if transcript_id in line:

                line = line.rstrip("\n")

                cols = line.split("\t")

                return cols[1:]


def get_pvalue(transcript_id):

    with open("{}/40_DEGseq2/DEGseq2_isoform_result_cooksCutoff_FALSE/RSEM.isoform.counts.matrix.N_vs_V.DESeq2.DE_results.cooksCutoff_FALSE".format(result_dir)) as f:

        for line in f:

            if transcript_id in line:

                line = line.rstrip("\n")

                cols = line.split("\t")

                pvalue = cols[9]

                return pvalue





for fasta_file in fasta_files:

    pfam_id = os.path.basename(fasta_file).split(".")[0]

    with open(fasta_file) as f:

        for line in f:

            if re.match(">", line):

                transcript_id = line.lstrip(">")
                transcript_id = transcript_id.rstrip("\n")

                if re.match(">", line):

                    if "cooksCutoff_FALSE" in fasta_file:
                         
                        transcript_ids_cooksCutoff_FALSE = add_pfam_id(transcript_ids_cooksCutoff_FALSE, \
                                                                        transcript_id, \
                                                                        pfam_id)

                    else:

                        transcript_ids = add_pfam_id(transcript_ids, \
                                                      transcript_id, \
                                                      pfam_id)


for k, v in transcript_ids_cooksCutoff_FALSE.items():

    pfam_ids = ",".join(set(v))
    exps = "\t".join(get_expressions(k))
    pvalue = get_pvalue(k)

    if k in transcript_ids:

        print(k, "-", pfam_ids, ",".join(transcript_ids[k]), pvalue, exps, sep="\t")

    else:

        print(k, "+", pfam_ids, "", pvalue, exps, sep="\t")

