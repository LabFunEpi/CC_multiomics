
crgex=/CC_multiomics/scRNA_MULT/scRNA_cr.sh
rawdata=/rawdata
crout=/cellranger_output
crmult=/CC_multiomics/scRNA_MULT/scMULT_cr.sh

qsub -N PM001 -v id=PM001,sample=PM001,ecells=4000,readsdir=${rawdata}/scRNA/batch3/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N M4666 -v id=M4666,sample=M4666,ecells=4000,readsdir=${rawdata}/scRNA/batch3/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N C8 -v id=C8,sample=cv-C8,ecells=4000,readsdir=${rawdata}/scRNA/batch6/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N M3399 -v id=M3399,sample=M3399,ecells=4000,readsdir=${rawdata}/scRNA/batch3/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N C3 -v id=C3,sample=cv-C3,ecells=4000,readsdir=${rawdata}/scRNA/batch5/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N C1 -v id=C1,sample=cv-C1,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N M5167 -v id=M5167,sample=M5167,ecells=4000,readsdir=${rawdata}/scRNA/batch3/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N C4 -v id=C4,sample=cv-C4,ecells=4000,readsdir=${rawdata}/scRNA/batch5/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J072 -v id=072,sample=CV-072,ecells=3960,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J984 -v id=984,sample=CV-984,ecells=3278,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B36 -v id=B36,sample=B36,ecells=4000,readsdir=${rawdata}/scRNA/batch8/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B23 -v id=B23,sample=cv-B23,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B39 -v id=B39,sample=B39,ecells=4000,readsdir=${rawdata}/scRNA/batch8/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B33 -v id=B33,sample=B33,ecells=4000,readsdir=${rawdata}/scRNA/batch7/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J454 -v id=454,sample=CV-454,ecells=4372,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B31 -v id=B31,sample=B31,ecells=4000,readsdir=${rawdata}/scRNA/batch7/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J789 -v id=789,sample=CV-789,ecells=3245,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B34 -v id=B34,sample=B34,ecells=4000,readsdir=${rawdata}/scRNA/batch8/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B38 -v id=B38,sample=B38,ecells=4000,readsdir=${rawdata}/scRNA/batch8/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J570 -v id=570,sample=CV-570,ecells=3135,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J120 -v id=120,sample=CV-120,ecells=3399,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B32 -v id=B32,sample=B32,ecells=4000,readsdir=${rawdata}/scRNA/batch7/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N B30 -v id=B30,sample=B30,ecells=4000,readsdir=${rawdata}/scRNA/batch7/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J890 -v id=890,sample=CV-890,ecells=3630,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A6 -v id=A6,sample=cv-A6,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A25 -v id=A25,sample=cv-A25,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A42 -v id=A42,sample=A42,ecells=4000,readsdir=${rawdata}/scRNA/batch7/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A16 -v id=A16,sample=cv-A16,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N J108 -v id=108,sample=CV-108,ecells=3245,readsdir=${rawdata}/scRNA/batch2/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A5 -v id=A5,sample=cv-A5,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A12 -v id=A12,sample=cv-A12,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A40 -v id=A40,sample=A40,ecells=4000,readsdir=${rawdata}/scRNA/batch7/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A30 -v id=A30,sample=cv-A30,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}
qsub -N A26 -v id=A26,sample=cv-A26,ecells=4000,readsdir=${rawdata}/scRNA/batch4/,outdir=${crout}/scRNA_cr_6.0.1 ${crgex}

qsub -N C2 -v id=C2,libraries=${crout}/scMULT_v2/C2.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N C5 -v id=C5,libraries=${crout}/scMULT_v2/C5.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N C8 -v id=C8,libraries=${crout}/scMULT_v2/C8.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N C3 -v id=C3,libraries=${crout}/scMULT_v2/C3.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N C1 -v id=C1,libraries=${crout}/scMULT_v2/C1.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N C4 -v id=C4,libraries=${crout}/scMULT_v2/C4.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N B36 -v id=B36,libraries=${crout}/scMULT_v2/B36.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N B23 -v id=B23,libraries=${crout}/scMULT_v2/B23.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N B39 -v id=B39,libraries=${crout}/scMULT_v2/B39.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N B31 -v id=B31,libraries=${crout}/scMULT_v2/B31.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N B34 -v id=B34,libraries=${crout}/scMULT_v2/B34.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N B30 -v id=B30,libraries=${crout}/scMULT_v2/B30.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N A6 -v id=A6,libraries=${crout}/scMULT_v2/A6.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N A25 -v id=A25,libraries=${crout}/scMULT_v2/A25.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N A16 -v id=A16,libraries=${crout}/scMULT_v2/A16.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N A5 -v id=A5,libraries=${crout}/scMULT_v2/A5.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N A12 -v id=A12,libraries=${crout}/scMULT_v2/A12.csv,outdir=${crout}/scMULT_v2 ${crmult}
qsub -N A40 -v id=A40,libraries=${crout}/scMULT_v2/A40.csv,outdir=${crout}/scMULT_v2 ${crmult}

######################### Renaming filenames with patient IDs used in manuscript ################

cd /rawdata/scRNA/batch3
for i in PM001*; do cp "$i" /GEO/1_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in M4666*; do cp "$i" /GEO/4_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch6
for i in cv-C8*; do cp "$i" /GEO/5_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch3
for i in M3399*; do cp "$i" /GEO/6_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch5
for i in cv-C3*; do cp "$i" /GEO/7_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch4
for i in cv-C1*; do cp "$i" /GEO/8_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch3
for i in M5167*; do cp "$i" /GEO/9_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch5
for i in cv-C4*; do cp "$i" /GEO/12_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch2
for i in CV-072*; do cp "$i" /GEO/13_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in CV-984*; do cp "$i" /GEO/14_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch8
for i in B36*; do cp "$i" /GEO/15_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch4
for i in cv-B23*; do cp "$i" /GEO/16_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch8
for i in B39*; do cp "$i" /GEO/18_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch7
for i in B33*; do cp "$i" /GEO/19_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch2
for i in CV-454*; do cp "$i" /GEO/20_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch7
for i in B31*; do cp "$i" /GEO/21_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch2
for i in CV-789*; do cp "$i" /GEO/23_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch8
for i in B34*; do cp "$i" /GEO/24_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in B38*; do cp "$i" /GEO/25_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch2
for i in CV-570*; do cp "$i" /GEO/27_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in CV-120*; do cp "$i" /GEO/28_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch7
for i in B32*; do cp "$i" /GEO/29_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in B30*; do cp "$i" /GEO/30_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch2
for i in CV-890*; do cp "$i" /GEO/34_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch4
for i in cv-A6*; do cp "$i" /GEO/35_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in cv-A25*; do cp "$i" /GEO/36_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch7
for i in A42*; do cp "$i" /GEO/37_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch4
for i in cv-A16*; do cp "$i" /GEO/38_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch2
for i in CV-108*; do cp "$i" /GEO/40_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch4
for i in cv-A5*; do cp "$i" /GEO/41_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in cv-A12*; do cp "$i" /GEO/42_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch7
for i in A40*; do cp "$i" /GEO/43_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
cd /rawdata/scRNA/batch4
for i in cv-A30*; do cp "$i" /GEO/45_scRNA_"$(echo $i | cut -d"_" -f2-)"; done
for i in cv-A26*; do cp "$i" /GEO/46_scRNA_"$(echo $i | cut -d"_" -f2-)"; done

cd /rawdata/scMULT/RNA
for i in C2_GEX_11*; do cp "$i" /GEO/1_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in C5_GEX_14*; do cp "$i" /GEO/4_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in C8_GEX_17*; do cp "$i" /GEO/5_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in C3_GEX_12*; do cp "$i" /GEO/7_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in C1_GEX_10*; do cp "$i" /GEO/8_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in C4_GEX_13*; do cp "$i" /GEO/12_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in B36_GEX_23*; do cp "$i" /GEO/15_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in B23_GEX_9*; do cp "$i" /GEO/16_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in B39_GEX_25*; do cp "$i" /GEO/18_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in B31_GEX_20*; do cp "$i" /GEO/21_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in B34_GEX_22*; do cp "$i" /GEO/24_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in B30_GEX_19*; do cp "$i" /GEO/30_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in A6_GEX_1*; do cp "$i" /GEO/35_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in A25_GEX_2*; do cp "$i" /GEO/36_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in A16_GEX_4*; do cp "$i" /GEO/38_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in A5_GEX_8*; do cp "$i" /GEO/41_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in A12_GEX_5*; do cp "$i" /GEO/42_GEX_"$(echo $i | cut -d"_" -f4-)"; done
for i in A40_GEX_18*; do cp "$i" /GEO/43_GEX_"$(echo $i | cut -d"_" -f4-)"; done

cd /rawdata/scMULT/ATAC
for i in C2_ATAC_11*; do cp "$i" /GEO/1_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in C5_ATAC_14*; do cp "$i" /GEO/4_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in C8_ATAC_17*; do cp "$i" /GEO/5_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in C3_ATAC_12*; do cp "$i" /GEO/7_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in C1_ATAC_10*; do cp "$i" /GEO/8_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in C4_ATAC_13*; do cp "$i" /GEO/12_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in B36_ATAC_23*; do cp "$i" /GEO/15_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in B23_ATAC_9*; do cp "$i" /GEO/16_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in B39_ATAC_25*; do cp "$i" /GEO/18_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in B31_ATAC_20*; do cp "$i" /GEO/21_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in B34_ATAC_22*; do cp "$i" /GEO/24_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in B30_ATAC_19*; do cp "$i" /GEO/30_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in A6_ATAC_1*; do cp "$i" /GEO/35_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in A25_ATAC_2*; do cp "$i" /GEO/36_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in A16_ATAC_4*; do cp "$i" /GEO/38_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in A5_ATAC_8*; do cp "$i" /GEO/41_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in A12_ATAC_5*; do cp "$i" /GEO/42_ATAC_"$(echo $i | cut -d"_" -f4-)"; done
for i in A40_ATAC_18*; do cp "$i" /GEO/43_ATAC_"$(echo $i | cut -d"_" -f4-)"; done


