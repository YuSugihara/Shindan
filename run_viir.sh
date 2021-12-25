#!/usr/bin/env bash -e


######################### Parameters #####################################
OUT_DIR=$1
FASTQ_LIST=$2
PVALUE=$3
PFAM_ID_LIST=$4
N_THREADS=$5
MAX_MEMORY=$6
SS_LIB_TYPE=$7
ADAPTER_FASTA=$8
##########################################################################


SCRIPT_DIR=$(cd $(dirname $0); pwd)

mkdir -p ${OUT_DIR}

cd ${OUT_DIR}

OUT_DIR=`pwd`

mkdir -p ${OUT_DIR}/00_fastq

cd ${OUT_DIR}/00_fastq

if [ ${ADAPTER_FASTA} = "Default_adapter" ]
then

    wget https://raw.githubusercontent.com/YuSugihara/ViiR/master/adapters.fasta \
         -O ${OUT_DIR}/00_fastq/adapter.fasta

    ADAPTER_FASTA=${OUT_DIR}/00_fastq/adapter.fasta

fi


FASTQ_CNT=0
TRINITY_LEFT=""
TRINITY_RIGHT=""


while read LINE || [ -n "${LINE}" ]
do

    COLS=(${LINE})

    SAMPLE_TYPE=${COLS[0]}
    FASTQ1=${SCRIPT_DIR}/${COLS[1]}
    FASTQ2=${SCRIPT_DIR}/${COLS[2]}


    mkdir -p ${OUT_DIR}/00_fastq/${SAMPLE_TYPE}_${FASTQ_CNT}

    PREFIX=${OUT_DIR}/00_fastq/${SAMPLE_TYPE}_${FASTQ_CNT}/${SAMPLE_TYPE}_${FASTQ_CNT}

    trimmomatic PE -threads ${N_THREADS} -phred33 \
    ${FASTQ1} \
    ${FASTQ2} \
    ${PREFIX}.1.trimmed.fastq.gz \
    ${PREFIX}.1.unpaired.trimmed.fastq.gz \
    ${PREFIX}.2.trimmed.fastq.gz \
    ${PREFIX}.2.unpaired.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTER_FASTA}:2:30:10 \
    LEADING:20 \
    TRAILING:20 \
    SLIDINGWINDOW:4:15 \
    MINLEN:75

    echo "${SAMPLE_TYPE}    \
          ${SAMPLE_TYPE}${FASTQ_CNT}    \
          ${PREFIX}.1.trimmed.fastq.gz    \
          ${PREFIX}.2.trimmed.fastq.gz" >> ${OUT_DIR}/00_fastq/fastq_list.txt


    if [ ${FASTQ_CNT} = 0 ]
    then

        TRINITY_LEFT="${PREFIX}.1.trimmed.fastq.gz"
        TRINITY_RIGHT="${PREFIX}.2.trimmed.fastq.gz"

    else

        TRINITY_LEFT="${TRINITY_LEFT},${PREFIX}.1.trimmed.fastq.gz"
        TRINITY_RIGHT="${TRINITY_RIGHT},${PREFIX}.2.trimmed.fastq.gz"

    fi


    FASTQ_CNT=$((FASTQ_CNT+1))


done < ${SCRIPT_DIR}/${FASTQ_LIST}


mkdir -p ${OUT_DIR}/10_trinity


if [ ${SS_LIB_TYPE} = "No" ]
then

    Trinity --seqType fq \
            --max_memory ${MAX_MEMORY} \
            --left ${TRINITY_LEFT} \
            --right ${TRINITY_RIGHT} \
            --output ${OUT_DIR}/10_trinity/trinity_assembly \
            --CPU ${N_THREADS} \
            --full_cleanup

else

    Trinity --seqType fq \
            --max_memory ${MAX_MEMORY} \
            --left ${TRINITY_LEFT} \
            --right ${TRINITY_RIGHT} \
            --output ${OUT_DIR}/10_trinity/trinity_assembly \
            --SS_lib_type ${SS_LIB_TYPE} \
            --CPU ${N_THREADS} \
            --full_cleanup

fi


mkdir -p ${OUT_DIR}/20_estimate_abundance

cd ${OUT_DIR}/20_estimate_abundance


align_and_estimate_abundance.pl --transcripts ${OUT_DIR}/10_trinity/trinity_assembly.Trinity.fasta \
                                --gene_trans_map ${OUT_DIR}/10_trinity/trinity_assembly.Trinity.fasta.gene_trans_map \
                                --seqType fq \
                                --samples_file ${OUT_DIR}/00_fastq/fastq_list_for_DEG.txt \
                                --SS_lib_type FR \
                                --est_method RSEM \
                                --aln_method bowtie \
                                --coordsort_bam \
                                --trinity_mode \
                                --prep_reference \
                                --thread_count ${N_THREADS}



mkdir -p ${OUT_DIR}/30_count_matrix

cd ${OUT_DIR}/30_count_matrix


abundance_estimates_to_matrix.pl --est_method RSEM \
                                 --gene_trans_map ${OUT_DIR}/10_trinity/trinity_assembly.Trinity.fasta.gene_trans_map \
                                 --name_sample_by_basedir \
                                 --out_prefix RSEM \
                                 ${OUT_DIR}/20_estimate_abundance/*/RSEM.isoforms.results



mkdir -p ${OUT_DIR}/40_DEGseq2

cd ${OUT_DIR}/40_DEGseq2

cut -f 1,2 ${OUT_DIR}/00_fastq/fastq_list.txt > fastq_list_for_DESeq2.txt


run_DE_analysis.pl --matrix ${OUT_DIR}/30_count_matrix/RSEM.gene.counts.matrix \
                   --method DESeq2 \
                   --samples_file fastq_list_for_DESeq2.txt \
                   --output ${OUT_DIR}/40_DEGseq2/DEGseq2_gene_result


run_DE_analysis.pl --matrix ${OUT_DIR}/30_count_matrix/RSEM.isoform.counts.matrix \
                   --method DESeq2 \
                   --samples_file fastq_list_for_DESeq2.txt \
                   --output ${OUT_DIR}/40_DEGseq2/DEGseq2_isoform_result


function make_cooksCutoff_FALSE() {

    DATA_TYPE=$1

    DIR_PREFIX=${OUT_DIR}/40_DEGseq2/DEGseq2_${DATA_TYPE}_result
    FILE_PREFIX=RSEM.${DATA_TYPE}.counts.matrix.V_vs_N.DESeq2

    mkdir -p ${DIR_PREFIX}_cooksCutoff_FALSE

    cat ${DIR_PREFIX}/${FILE_PREFIX}.Rscript | \
    sed "s/contrast/contrast\,\ cooksCutoff\=FALSE/g" \
    > ${DIR_PREFIX}_cooksCutoff_FALSE/${FILE_PREFIX}.cooksCutoff_FALSE.Rscript

    Rscript ${DIR_PREFIX}_cooksCutoff_FALSE/${FILE_PREFIX}.cooksCutoff_FALSE.Rscript

    mv ${DIR_PREFIX}_cooksCutoff_FALSE/${FILE_PREFIX}.DE_results \
       ${DIR_PREFIX}_cooksCutoff_FALSE/${FILE_PREFIX}.DE_results.cooksCutoff_FALSE

}


function get_prefix() {

    DATA_TYPE=$1
    COOKSCUTOFF=$2

    if [ ${COOKSCUTOFF} = "TRUE" ]
    then

        DIR_PREFIX=${OUT_DIR}/40_DEGseq2/DEGseq2_${DATA_TYPE}_result
        FILE_PREFIX=RSEM.${DATA_TYPE}.counts.matrix.V_vs_N.DESeq2

    else

        DIR_PREFIX=${OUT_DIR}/40_DEGseq2/DEGseq2_${DATA_TYPE}_result_cooksCutoff_FALSE
        FILE_PREFIX=RSEM.${DATA_TYPE}.counts.matrix.V_vs_N.DESeq2.cooksCutoff_FALSE

    fi

    PREFIX=${DIR_PREFIX}/${FILE_PREFIX}

    echo ${PREFIX}

}


function get_significant_list() {

    DATA_TYPE=$1
    COOKSCUTOFF=$2

    PREFIX=`get_prefix ${DATA_TYPE} ${COOKSCUTOFF}`

    cat ${PREFIX} | \
    awk '{if ($10 < '${PVALUE}') print $0}' | \
    awk '{if ($7 > 0) print $0}' | \
    cut -f 1 \
    > ${PREFIX}.significant_${DATA_TYPE}s

}


function get_significant_fasta() {

    COOKSCUTOFF=$1


    PREFIX=`get_prefix isoform ${COOKSCUTOFF}`

    samtools faidx -r ${PREFIX}.significant_${DATA_TYPE}s
                      ${OUT_DIR}/10_trinity/trinity_assembly.Trinity.fasta \
                    > ${PREFIX}.significant_${DATA_TYPE}s.fasta

    esl-translate ${PREFIX}.significant_${DATA_TYPE}s.fasta \
    sed "s/\ source\=/-/g" > ${PREFIX}.significant_${DATA_TYPE}s.AA.fasta

}


function get_hmmscan_list() {

    HMMSCAN_RESULT=$1

    cat ${HMMSCAN_RESULT} | \
    grep -v "#" | \
    awk '{print $3}' | \
    ut -f 2 -d "-" | \
    sort | \
    uniq

}


function get_hmmscan_fasta() {

    ISOFORM_LIST=$1
    COOKSCUTOFF=$2

    if [ `cat ${ISOFORM_LIST} | wc -l` -gt 0 ]
    then

        if [ ${COOKSCUTOFF} = "TRUE" ]
        then

            samtools faidx -r ${ISOFORM_LIST}
                            ${OUT_DIR}/10_trinity/trinity_assembly.Trinity.fasta \
                            > ${OUT_DIR}/60_fasta/${PFAM_ID}.fasta

        else

            samtools faidx -r ${ISOFORM_LIST}
                            ${OUT_DIR}/10_trinity/trinity_assembly.Trinity.fasta \
                            > ${OUT_DIR}/60_fasta/${PFAM_ID}.cooksCutoff_FALSE.fasta

        fi

    fi

}


export make_cooksCutoff_FALSE
export get_prefix
export get_significant_list
export get_significant_fasta
export get_hmmscan_list
export get_hmmscan_fasta


make_cooksCutoff_FALSE "gene" 
make_cooksCutoff_FALSE "isoform"


get_significant_list "gene"    "TRUE"
get_significant_list "gene"    "FALSE"
get_significant_list "isoform" "TRUE"
get_significant_list "isoform" "FALSE"


get_significant_fasta "TRUE"
get_significant_fasta "FALSE"


mkdir -p ${OUT_DIR}/50_hmmer
mkdir -p ${OUT_DIR}/60_fasta

cd ${OUT_DIR}/50_hmmer


if [ ${PFAM_ID_LIST} = "Default_list" ]
then

    wget https://raw.githubusercontent.com/YuSugihara/ViiR/master/Pfam_IDs_list.txt

    PFAM_ID_LIST=${OUT_DIR}/50_hmmer/Pfam_IDs_list.txt

fi


while read PFAM_ID || [ -n "${PFAM_ID}" ]
do

    wget https://pfam.xfam.org/family/${PFAM_ID}/hmm \
         -O ${PFAM_ID}.hmm


    hmmpress ${PFAM_ID}.hmm


    PREFIX=`get_prefix isoform TRUE`

    hmmscan  --cpu ${N_THREADS} \
             --tblout ${PFAM_ID}_hmmscan.txt \
             ${PFAM_ID}.hmm \
             ${PREFIX}.significant_isoforms.AA.fasta \
             1> /dev/null


    get_hmmscan_list ${PFAM_ID}_hmmscan.txt \
                   > ${PFAM_ID}_hmmscan.isoform_list.txt

    get_hmmscan_fasta ${PFAM_ID}_hmmscan.isoform_list.txt "TRUE"


    PREFIX=`get_prefix isoform FALSE`

    hmmscan  --cpu ${N_THREADS} \
             --tblout ${PFAM_ID}_cooksCutoff_FALSE_hmmscan.txt \
             ${PFAM_ID}.hmm \
             ${PREFIX}.significant_isoforms.AA.fasta \
             1> /dev/null

    get_hmmscan_list ${PFAM_ID}_cooksCutoff_FALSE_hmmscan.txt \
                   > ${PFAM_ID}_cooksCutoff_FALSE_hmmscan.isoform_list.txt

    get_hmmscan_fasta ${PFAM_ID}_cooksCutoff_FALSE_hmmscan.isoform_list.txt "FALSE"

done < ${PFAM_ID_LIST}


mv ${SCRIPT_DIR}/run_viir.sh ${OUT_DIR}
