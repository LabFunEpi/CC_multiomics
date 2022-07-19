import scrublet as scr
import scipy.io
import numpy as np
import os
import sys

if len(sys.argv) <= 1:
 sys.exit("Usage: python scrub.py <indir>")

mtx = sys.argv[1] + "/matrix.mtx"
tsv = sys.argv[1] + "/features.tsv"

counts_matrix = scipy.io.mmread(mtx).T.tocsc()
genes = np.array(scr.load_genes(tsv, delimiter='\t', column=1))

print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
print('Number of genes in gene list: {}'.format(len(genes)))

scrub = scr.Scrublet(counts_matrix)
doublet_scores, predicted_doublets = scrub.scrub_doublets()

f = open(sys.argv[1] + "/scrublet.tsv", "w")
for i in range(0, len(scrub.doublet_scores_obs_)):
 f.write(str(scrub.doublet_scores_obs_[i]) + "\t" + str(scrub.doublet_scores_sim_[i]) + "\t" + str(predicted_doublets[i]) + "\n")
f.close()

