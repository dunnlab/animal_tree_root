#!/usr/bin/env python
from glob import glob
import pandas as pd

colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
qseq_colnames = ['mscript1', 'dataset1', 'taxon1', 'taxid1', 'part1']
sseq_colnames = ['mscript2', 'dataset2', 'taxon2', 'taxid2', 'part2']
pident_cutoff = 95.0

for i, infile in enumerate(sorted(glob('blast/results/*_permissive.txt'))):
    # first file
    df_tmp = pd.read_csv(infile, sep='\t', names=colnames)
    df_tmp = df_tmp[ df_tmp['pident']>pident_cutoff ]

    qseq = df_tmp['qseqid'].str.split(':', expand=True)
    qseq.columns = qseq_colnames
    sseq = df_tmp['sseqid'].str.split(':', expand=True)
    sseq.columns = sseq_colnames

    df_tmp = pd.concat((qseq, sseq, df_tmp), axis=1)
    if i ==0:
        df = df_tmp
    else:
        df = pd.concat((df, df_tmp))

parts_graph = df[ ['mscript1', 'dataset1', 'part1', 'mscript2', 'dataset2', 'part2'] ].drop_duplicates()
parts_graph.to_csv('blast/graphs/partitions_graph.tsv', sep='\t', index=False)
taxon_graph = df[ ['mscript1', 'dataset1', 'taxon1', 'mscript2', 'dataset2', 'taxon2'] ].drop_duplicates()
taxon_graph.to_csv('blast/graphs/taxon_graph.tsv', sep='\t', index=False)

cross_mscript = df[(df['taxon1'] == df['taxon2']) & (df['mscript1'] != df['mscript2'])]
cross_mscript.to_csv('blast/graphs/cross_mscript_taxa.tsv', sep='\t', index=False)
