#!/usr/bin/env python
from glob import glob
import pandas as pd
import networkx as nx

# set up for reading in the blast data
colnames = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
qseq_colnames = ['mscript1', 'dataset1', 'taxon1', 'taxid1', 'part1']
sseq_colnames = ['mscript2', 'dataset2', 'taxon2', 'taxid2', 'part2']
# cutoffs
pident_cutoff = 95.0
evalue_cutoff = 1e-25

# read in blast results
for i, infile in enumerate(sorted(glob('blast/results/*_permissive.txt'))):
    df_tmp = pd.read_csv(infile, sep='\t', names=colnames)
    df_tmp = df_tmp[(df_tmp['pident']>pident_cutoff) & (df_tmp['evalue']<evalue_cutoff)]
    if i ==0:
        df = df_tmp
    else:
        df = pd.concat((df, df_tmp))

# keep only lowest evalue for each pair of hits
unique_best_evalue = df.groupby(['qseqid', 'sseqid'])['evalue'].min().reset_index()
# explode the qseq columns
qseq = unique_best_evalue['qseqid'].str.split(':', expand=True)
sseq = unique_best_evalue['sseqid'].str.split(':', expand=True)
qseq.columns = qseq_colnames
sseq.columns = sseq_colnames
unique_best_evalue = pd.concat((qseq, sseq, unique_best_evalue), axis=1)

# only retain hits from different matrices/datasets
inter_mscript = unique_best_evalue[unique_best_evalue['dataset1'] != unique_best_evalue['dataset2']]
# count unique mscript:dataset:partition <-> mscript:dataset:partition 
groupby_columns = ['mscript1', 'dataset1', 'part1', 'mscript2', 'dataset2', 'part2']
parts = inter_mscript.groupby(groupby_columns).size().reset_index(name='count')

# construct graph from hits
parts_graph = nx.Graph()
parts_graph.add_edges_from([(tuple(x[:3]), tuple(x[3:6]), {'weight':int(x[6])}) for x in  parts.values])
# divide graph into separate components, sorted by size
components = nx.connected_components(parts_graph)
out = open('blast/graphs/partition_clusters.tsv', 'w')
print('\t'.join(['manuscript', 'dataset', 'partition_name','cluster_number']), file=out)
for i, component in enumerate(sorted(components, key=len, reverse=True)):
    [print('\t'.join(x+(str(i),)), file=out) for x in component]
out.close()
