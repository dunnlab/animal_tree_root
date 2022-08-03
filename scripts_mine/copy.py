# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 12:06:20 2022

@author: matil
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 10:15:14 2022

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
from time import sleep
import sys

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
            df = pd.concat((df, tmp_df),ignore_index=True)
    return df

# find and set component number
def set_co_nr(name):
    co_nr= partitions_df_expanded[partitions_df_expanded.matrix_partition.isin([name])].component_number
    if len(co_nr) > 0:
        co_nr = int(co_nr)
    else:
        co_nr = -1
    return co_nr

# Add component number to dataframe    
def add_co_nr(i,df,col_name,col_name2,name):
#    tempdf = df.copy()
#    com_nr = set_co_nr(name)
#    tempdf.at[i,col_name]=com_nr
#    tempdf.at[i,col_name2]=name
#    return tempdf
    com_nr = set_co_nr(name)
    df.at[i,col_name]=com_nr
    df.at[i,col_name2]=name
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

ava_df = ava_df[['sseqid','qseqid']]

#Import partitions
partitions_df_expanded = pd.read_csv("../files_mine/component_partition_expansded.csv", sep="\t")

import progressbar
from time import sleep
bar = progressbar.ProgressBar(maxval=len(ava_df), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])

bar.start()
for i,qname,sname in zip(ava_df.index.values.astype(int), ava_df.qseqid,ava_df.sseqid):
    ava_df = add_co_nr(i,ava_df,'scomp_nr','sname',sname)
    #ava_df = add_co_nr(i,ava_df,'qcomp_nr','qname',qname)
    bar.update(i+1)
    sleep(0.1)
bar.finish()
    
    
#ava_df_s = ava_df.apply(lambda row : add_co_nr(i,ava_df,'scomp_nr','sname',ava_df.sseqid), axis = 1)
    