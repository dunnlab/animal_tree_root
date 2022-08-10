# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:12:54 2022

@author: matil
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:44:27 2022

@author: matil
"""

#Idea: add component number for each seq based on the blast query name 
import pandas as pd
from Bio import SeqIO
from compare_comp_nrs.py import add_co_nr,process_df
from glob import glob
import numpy as np
# read in blast results, function copied from blast_partitions_graph.py
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
#add component number
def edit_df(df,partition_df):
    for i,name in zip(df.index.values.astype(int), df.sseqid):
        df = add_co_nr(i,df,'comp_nr','name',name,partition_df)
        print('\rProcessing row:'+str(i)+'/'+str(len(df)),end="")
    return df

#read fasta headers into dataframe
def read_blast_query(dirpattern,colnames):
    l1 = []
    l2 = []
    for i, infile in enumerate(sorted(glob(dirpattern))):
        with open(infile, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                l1.append(record.name)
                l2.append(record.seq)
                #print(record)
#        if i ==0:
#            l = l
#        else:
#            df = pd.concat((l, tmp_l),ignore_index=True)
    df = pd.DataFrame({colnames[0]:l1,colnames[1]:l2})
    return df

#Read blast results, add component nrs only for result not for query
ava_pident_cutoff = 50
ava_evalue_cutoff = 1e-5
ava_colnames = ['qseqid', 'sseqid', 'pident', 'length', 
                'mismatch', 'gapopen', 'qstart', 'qend', 
                'sstart', 'send', 'evalue', 'bitscore']
ava_df = read_blast_results('../add_new_data_test/diamond_results/*_permissive.txt',   #read only Dunn and Philippe) 
                            ava_colnames, 
                            ava_pident_cutoff, 
                            ava_evalue_cutoff)
# filter out self-matches
ava_df_sm = ava_df[ava_df['qseqid'] != ava_df['sseqid']]

#unify names with partition map
ava_df_sm = ava_df["qseqid"].str.split(":", n = 3, expand = True)
ava_df["qseqid"]= ava_df_sm[0]+ava_df_sm[3]
ava_df_sm = ava_df["sseqid"].str.split(":", n = 3, expand = True)
ava_df["sseqid"]= ava_df_sm[0]+ava_df_sm[3]

#keep only columns with names
ava_df = ava_df[['sseqid','qseqid']]

partitions_df_expanded_r = pd.read_csv("../files_mine/part_map_glob.csv", sep="\t")
#ava_df = ava_df.drop([0,3],axis=1)

#process partition map
partitions_df_expanded_r,group_df_r = process_df(partitions_df_expanded_r)

#add component numbers
df_processed = edit_df(ava_df.copy(),partitions_df_expanded_r)
