#!/bin/bash
#$ -q 1-day
#$ -o $JOB_NAME.stdout
#$ -e $JOB_NAME.stderr
#$ -M wazim.ismail@gmail.com
#$ -m abe
#$ -pe threaded 16
#$ -l h_vmem=8G

. $HOME/.bash_profile

myhome=/
softwares=${myhome}/softwares
references=${myhome}/references

cr_path=${softwares}/cellranger-6.0.1/bin

ref=${references}/transcriptomes/refdata-gex-GRCh38-2020-A

cd ${outdir}
${cr_path}/cellranger count --id=${id} \
                            --transcriptome=${ref} \
                            --fastqs=${readsdir} \
                            --sample=${sample} \
                            --expect-cells=${ecells} \
                            --localcores=16 \
                            --localmem=115

