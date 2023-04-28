#!/bin/bash

### exit when error
set -e

### Local Variables


TO_RUN=$1
FROM_IDX=$2
TO_IDX=$3
BASE_UTILS=$4

if [[ "${FROM_IDX}" == "" ]] || [[ "${TO_IDX}" == "" ]] ; then
    echo "Usage: ${0} SRA|TCGA from_idx to_idx"
    exit 0
fi

if [[ "${TO_RUN}" == "SRA" ]] || [[ "${TO_RUN}" == "TCGA" ]]; then 
    echo "Running ${TO_RUN} from "${FROM_IDX}" to "${TO_IDX}""
else
    echo "Usage: ${0} SRA|TCGA from_idx to_idx"
    exit 0
fi




BASE_DIR=$(realpath ./ )

OUT_DIR="${BASE_DIR}/output"


ngc="${BASE_UTILS}/tokens/dbGAP.ngc"
token_tcga="${BASE_UTILS}/tokens/TCGA.txt"



TCGA_manifest="${BASE_DIR}/metadata/gdc_manifest.TCGA_prostate.txt"
SRA_ids="${BASE_DIR}/metadata/SRR_ids.txt"  

### References 

GENCODE_VERSION="39"

reference_gtf="${BASE_UTILS}/gencode_v${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf"
reference_genome="${BASE_UTILS}/gencode_v${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa"

### Process resources

pipeline_version="3.6"
threads="16"
max_mem="80.GB"

MAX_SRA_TRY=3
## run only fastqc
skip_processes="--skip_stringtie --skip_preseq --deseq2_vst false --skip_bigwig  --skip_dupradar  --skip_qualimap  --skip_rseqc --skip_biotype_qc --skip_deseq2_qc --skip_multiqc --skip_markduplicates "

### Some files check


if [ ! -f ${reference_gtf} ]; then
    wget -O ${reference_gtf}.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/gencode.v${GENCODE_VERSION}.annotation.gtf.gz 
    gzip -d ${reference_gtf}.gz
fi



if [ ! -f ${reference_genome} ]; then
    wget -O ${reference_genome}.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_VERSION}/GRCh38.primary_assembly.genome.fa.gz
    gzip -d ${reference_genome}.gz
fi

### The STAR references will be generated on the first run if not exsist in the reference gencode directory
STAR_REF="${BASE_UTILS}/gencode_v${GENCODE_VERSION}/STAR/"


### Temporary fastq location ( each time remove processed raw data )
fastq_dir="${BASE_DIR}/fastq"
mkdir -p ${fastq_dir} 
## useful vars

timestamp=$(date +%d_%m_%y_%H_%M_%S)
tmp_dir="${BASE_DIR}/tmp_${timestamp}"

export NXF_WORK="${BASE_DIR}/work_${timestamp}"

mkdir -p ${tmp_dir}

mkdir -p ${OUT_DIR}
######################## SRA
LINE_N=0
if [[ "${TO_RUN}" == "SRA" ]] ; then
while read line; do
    if [[ "${line}" != "" ]]; then
        LINE_N=$((LINE_N + 1 ))
        SRR=${line}
        if [[ ! -d ${OUT_DIR}/${SRR}/ ]] && [ $LINE_N -ge ${FROM_IDX} ] && [ $LINE_N -lt ${TO_IDX} ] ; then
            echo "${LINE_N} : ${SRR}"
            if [ -f ${fastq_dir}/${SRR}_1.fastq.gz ]; then
                fq1=${fastq_dir}/${SRR}_1.fastq.gz                    
                if [ -f ${fastq_dir}/${SRR}_2.fastq.gz ]; then
                   fq2=${fastq_dir}/${SRR}_2.fastq.gz
                else
                   fq2=""
                fi
            else
                if [ -f ${fastq_dir}/${SRR}.fastq.gz ] ; then 
                    fq1=${fastq_dir}/${SRR}.fastq.gz
                    fq2=""
                else
                    TEST_SRA=0
                    while [ ${TEST_SRA} -lt ${MAX_SRA_TRY} ] ; do
                        echo "Try ${TEST_SRA} "
                        rm -fr ${fastq_dir}/${SRR}*
                        prefetch -O ${fastq_dir}/ --ngc ${ngc} ${SRR} && fasterq-dump --split-3 --outdir ${fastq_dir}/ ${fastq_dir}/${SRR} && TEST_SRA=${MAX_SRA_TRY} ||  echo "ERROR FOR ${SRR}" && TEST_SRA=$((TEST_SRA+1))
                    done
                    if [ -f ${fastq_dir}/${SRR}_1.fastq ]; then
                        gzip ${fastq_dir}/${SRR}_1.fastq
                        fq1=${fastq_dir}/${SRR}_1.fastq.gz
                        gzip ${fastq_dir}/${SRR}_2.fastq
                        fq2=${fastq_dir}/${SRR}_2.fastq.gz
                    else
                        if [ -f ${fastq_dir}/${SRR}.fastq ]; then
                            gzip ${fastq_dir}/${SRR}.fastq
                            fq1=${fastq_dir}/${SRR}.fastq.gz 
                            fq2=""
                        else   
                            fq1=""
                            fq2=""
                        fi
                    fi
                fi
            fi
            if [[ "${fq1}" == "" ]] ;then
                echo "${SRR}" >> ./missing_SRA.txt
            else
                echo -e "sample,fastq_1,fastq_2,strandedness\n${SRR},${fq1},${fq2},unstranded" > ${tmp_dir}/${SRR}.csv
                if [ -d ${STAR_REF} ]; then
                    parameters="--star_index ${STAR_REF}"
                else
                    parameters="--save_reference"
                fi
                nextflow run nf-core/rnaseq -profile singularity -r ${pipeline_version} --aligner star_salmon --gencode --max_cpus ${threads} --max_memory ${max_mem} --input  ${tmp_dir}/${SRR}.csv  ${skip_processes} \
                    --outdir ${OUT_DIR}/${SRR}/  --fasta ${reference_genome}  --gtf ${reference_gtf}   ${parameters}
                if [ ! -d ${STAR_REF} ]; then
                    mv ${OUT_DIR}/${SRR}/genome/index/star/ ${STAR_REF}
                fi
                ### cleanup
                rm -fr ${fastq_dir}/${SRR}*
                mv ${OUT_DIR}/${SRR}/star_salmon/salmon_tx2gene.tsv ./
                rm -fr ${OUT_DIR}/${SRR}/genome/ ${OUT_DIR}/${SRR}/star_salmon/*.bam  ${OUT_DIR}/${SRR}/star_salmon/*.bai  ${OUT_DIR}/${SRR}/star_salmon/*.rds ${OUT_DIR}/${SRR}/star_salmon/*scaled.tsv 
                if (( LINE_N % 10 == 9 )) ; then
                    rm -fr ${NXF_WORK}
                fi
           fi
        fi
    fi
    
done < ${SRA_ids}

rm -fr ${tmp_dir}
exit 0
fi




########################    TCGA

header="id	filename	md5	size	state"

while read line ; do
     if [[ "${line}" != "" && ${line} != "${header}" ]] ; then
        LINE_N=$((LINE_N + 1 ))
        SRR=$(echo "${line}" | awk -F '\t' '{print $1}')
        if [[ ! -d ${OUT_DIR}/${SRR}/ ]] && [ $LINE_N -ge ${FROM_IDX} ] && [ $LINE_N -lt ${TO_IDX} ] ; then
            fname=$(echo "${line}" | awk -F '\t' '{print $2}')
            echo "LINE: ${LINE_N} ${fname}"            
            echo -e "${header}\n${line}" > ${tmp_dir}/${SRR}.manifest.txt
            gdc-client download -t ${token_tcga} -d ${fastq_dir} -m ${tmp_dir}/${SRR}.manifest.txt
            bam_file="${fastq_dir}/${SRR}/${fname}"
            if [ -f ${bam_file} ]; then
                if [ -f ${fastq_dir}/${SRR}_1.fastq.gz ]; then
                    fq1=${fastq_dir}/${SRR}_1.fastq.gz                    
                    if [ -f ${fastq_dir}/${SRR}_2.fastq.gz ]; then
                        fq2=${fastq_dir}/${SRR}_2.fastq.gz
                    else
                        fq2=""
                    fi
                else
                    samtools sort -@ ${threads} -n ${bam_file} | samtools fastq -@ ${threads} -s /dev/null -0 /dev/null -1 ${fastq_dir}/${SRR}_1.fastq -2 ${fastq_dir}/${SRR}_2.fastq
                    rm -fr ${fastq_dir}/${SRR}/
                    if [ -f ${fastq_dir}/${SRR}_1.fastq ] && [ -f ${fastq_dir}/${SRR}_2.fastq ]; then
                        gzip ${fastq_dir}/${SRR}_1.fastq
                        fq1=${fastq_dir}/${SRR}_1.fastq.gz
                        gzip ${fastq_dir}/${SRR}_2.fastq
                        fq2=${fastq_dir}/${SRR}_2.fastq.gz
                    else
                        gzip ${fastq_dir}/${SRR}_1.fastq
                        fq1=${fastq_dir}/${SRR}_1.fastq.gz
                        fq2=""
                    fi
                fi
            else
               echo "Error in downloading $SRR"
               echo $SRR >> ./missing_TCGA.txt
            fi
            echo -e "sample,fastq_1,fastq_2,strandedness\n${SRR},${fq1},${fq2},unstranded" > ${tmp_dir}/${SRR}.csv
            if [ -d ${STAR_REF} ]; then
                parameters="--star_index ${STAR_REF}"
            else
                parameters="--save_reference"
            fi
            nextflow run nf-core/rnaseq -profile singularity -r ${pipeline_version} --aligner star_salmon --gencode --max_cpus ${threads} --max_memory ${max_mem} --input  ${tmp_dir}/${SRR}.csv  ${skip_processes} \
                --outdir ${OUT_DIR}/${SRR}/  --fasta ${reference_genome}  --gtf ${reference_gtf}   ${parameters}
            if [ ! -d ${STAR_REF} ]; then
                mv ${OUT_DIR}/${SRR}/genome/index/star/ ${STAR_REF}
            fi
            ### cleanup
            rm ${fastq_dir}/${SRR}*
            mv ${OUT_DIR}/${SRR}/star_salmon/salmon_tx2gene.tsv ./
            rm -fr ${NXF_WORK}/*/*/*.bam
            rm -fr ${OUT_DIR}/${SRR}/genome/ ${OUT_DIR}/${SRR}/star_salmon/*.bam  ${OUT_DIR}/${SRR}/star_salmon/*.bai  ${OUT_DIR}/${SRR}/star_salmon/*.rds ${OUT_DIR}/${SRR}/star_salmon/*scaled.tsv 
        fi
    fi  
done < ${TCGA_manifest}

rm -fr ${tmp_dir}




