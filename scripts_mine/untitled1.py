# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 13:34:41 2022

@author: matil
"""

#in blast/diamond result (ava_df?) rename genes with component number and see if any components map to other components
#If you change the dataset_geneID to component number.
#Are any case from the diamond files (premissive) where different component numbers hit each other?
#!/usr/bin/env python
from glob import glob
# import matplotlib.pyplot as plt
import pandas as pd
from goatools import obo_parser


# read in blast results
def read_blast_results(dirpattern, colnames, pident_cutoff, evalue_cutoff):
    df = None
    for i, infile in enumerate(sorted(glob(dirpattern))):
        tmp_df = pd.read_csv(infile, sep='\t', names=colnames)

        tmp_df = tmp_df[((tmp_df['pident']>pident_cutoff) & 
                         (tmp_df['evalue']<evalue_cutoff))]
        if i ==0:
            df = tmp_df
        else:
            df = pd.concat((df, tmp_df))
    return df


# Read in All vs All results (ava) from diaomond blast runs
# cutoffs: chose to retain enough results to keep at least 1 hit for ~90% of partitions
ava_pident_cutoff = 50
ava_evalue_cutoff = 1e-5
ava_colnames = ['qseqid', 'sseqid', 'pident', 'length', 
                'mismatch', 'gapopen', 'qstart', 'qend', 
                'sstart', 'send', 'evalue', 'bitscore']
ava_df = read_blast_results('../reconciliation/blast/diamond_results/dunn_philippe/*_permissive.txt', 
                            ava_colnames, 
                            ava_pident_cutoff, 
                            ava_evalue_cutoff)

# filter out self-matches
ava_df_sm = ava_df[ava_df['qseqid'] != ava_df['sseqid']]       #Still includes missing partitions

ava_df_sm = ava_df["qseqid"].str.split(":", n = 3, expand = True)
ava_df["qseqid"]= ava_df_sm[0]+ava_df_sm[3]

ava_df_sm = ava_df["sseqid"].str.split(":", n = 3, expand = True)
ava_df["sseqid"]= ava_df_sm[0]+ava_df_sm[3]


partitions_df = pd.read_csv("../files_mine/component_partition.csv", sep="\t")

# to do: compare ids in ava to matrix-partition in partitions and add component numbers to ava as column

#for sseqid and for qseqid 
#if name (xseqid) in ava is in list of names in partitions (matrix_partition) add component nr as new column (xcomp_nr)

for i,name in zip(ava_df.index.values.astype(int), ava_df.sseqid):
    for co_nr , name_list in zip(partitions_df.component_number, partitions_df.matrix_partition):
        if name in name_list:
            ava_df.at[i,'s_comp_nr'] = co_nr

for i,name in zip(ava_df.index.values.astype(int), ava_df.qseqid):
    for co_nr , name_list in zip(partitions_df.component_number, partitions_df.matrix_partition):
        if name in name_list:
            ava_df.at[i,'q_comp_nr'] = co_nr
            
ava_df[['s_comp_nr','q_comp_nr']]

non_match_comp_df = ava_df[ava_df['s_comp_nr'] != ava_df['q_comp_nr']]