#!/usr/bin/env python
from glob import glob
import pandas as pd

# set up for reading in the blast data
colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'gene_info']
# cutoffs
pident_cutoff = 95.0
evalue_cutoff = 1e-25

for i, infile in enumerate(sorted(glob('./*_swissprot.txt'))):
    df_tmp = pd.read_csv(infile, sep='\t', names=colnames)
    df_tmp = df_tmp[(df_tmp['pident']>pident_cutoff) & (df_tmp['evalue']<evalue_cutoff)]
    if i ==0:
        df = df_tmp
    else:
        df = pd.concat((df, df_tmp))
df.to_csv('swissprot_best_hits.tsv', sep='\t', index=False)
