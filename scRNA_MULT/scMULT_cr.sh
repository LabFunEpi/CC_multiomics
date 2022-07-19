#!/bin/bash
#$ -q 1-day
#$ -o /research/labs/experpath/maia/m237371/logs/$JOB_NAME.stdout
#$ -e /research/labs/experpath/maia/m237371/logs/$JOB_NAME.stderr
#$ -M wazim.ismail@gmail.com
#$ -m abe
#$ -pe threaded 16
#$ -l h_vmem=8G

. $HOME/.bash_profile

maiahome=/
myhome=${maiahome}/m237371
softwares=${myhome}/softwares
references=${myhome}/references

cr_arc_path=${softwares}/cellranger-arc-2.0.0

ref=${references}/transcriptomes/refdata-cellranger-arc-GRCh38-2020-A-2.0.0

cd $outdir
${cr_arc_path}/cellranger-arc count --id=${id} \
                                    --reference=${ref} \
                                    --libraries=${libraries} \
                                    --localcores=16 \
                                    --localmem=115
