# Atlas generation scripts
Scripts to run in order to generate the Prostate Cancer Atlas.

### Requirements
The following requirements must be met:
 - Singularity ( Apptainer ) installed
 - Nextflow
 - fasterq-dump and prefetch ( sra-tools)
 - gdc-client
 - samtools
 - A folder (Utils directory) containing the tokens for restricted access databases under the paths ${BASE_UTILS}/tokens/dbGAP.ngc and ${BASE_UTILS}/tokens/TCGA.txt. The STAR reference will be generated in this folder with the first run and the references from GENCODE will be downloaded here.
 - In the folder where you run the 01_Download_and_process.sh there must be a subfolder called metadata with the information about the data to download, provided in this folder:
    - "${BASE_DIR}/metadata/gdc_manifest.TCGA_prostate.txt"
    - "${BASE_DIR}/metadata/SRR_ids.txt"
 - R with the requirements listed in the given scripts
## 01: Download and process

The first script download and process the raw FASTQ files downloaded from SRA or GDC ( TCGA ).
You can set the limits for each job at lines 52 and 53 with the variables 'threads' and 'max_mem'

### Run the script
To run all the samples, you can use the following commands:
```
./01_Download_and_process.sh SRA 0 100000 ~/Utils/
./01_Download_and_process.sh TCGA 0 100000 ~/Utils/
```
Where 0 is the lower limit, 100000 is the upper ( you can divide the work in batches ) and the ~/Utils/ is the BASE_UTILS folder containing the tokens.
The first process will generate a STAR reference in the Utils folder that will be reused by other processes, the first time wait that at least one of job is completed to run other batches.

### Output
A an ./output folder containing all the results. ( Salmon quantification, fastqc and trimgalore )

## 02: Collect data
Collect the data from the output folder and compress them. 

### Run the script
In the folder where you ran the previous script, and therefore containing the subfolder output, run
```
./02_Collect_data.sh
```

### Output
An export folder containing subfolders named with the sample ID and the compressed salmon's output quant.sf.gz.



