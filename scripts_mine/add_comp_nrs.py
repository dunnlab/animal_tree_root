# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:44:27 2022

@author: matil
"""

#Idea: add component number for each seq based on the blast query name 
import pandas as pd
from Bio import SeqIO
from compare_comp_nrs import add_co_nr,process_df
from glob import glob
import numpy as np

#add component number
def edit_df(df,partition_df):
    for i,name in zip(df.index.values.astype(int), df.name):
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

#import fasta headers from blast queries to use as name
infile = r'..\reconciliation\blast\queries\dunn_philippe\*.fa'
colnames = ["fasta_header","seq"]
headers_df = read_blast_query(infile,colnames)

#add name column
headers_df = headers_df.join(headers_df.fasta_header.str.split(':', expand = True))[[colnames[0],colnames[1],0,3]]
headers_df["name"] = headers_df[0] + headers_df[3]
#import partition map 
partitions_df_expanded_r = pd.read_csv("../manuscript/part_map_glob.csv", sep="\t")
headers_df = headers_df.drop([0,3],axis=1)

#process partition map
partitions_df_expanded_r,group_df_r = process_df(partitions_df_expanded_r)

#add component numbers
headers_df_processed = edit_df(headers_df.copy(),partitions_df_expanded_r)

#split dataframe on component nr
dfs = dict(tuple(headers_df_processed.groupby('comp_nr')))


#write outfile
#list of lists of turned into biopython record
records = []
for i in dfs:
    sep_rec = []
    for l,j in zip(dfs[i].fasta_header,dfs[i].sequence):
        sep_rec.append(SeqIO.SeqRecord(j, id = l,description="Component_nr:"+str(int(i))))
    records.append([sep_rec,i])
    
    
for i in range(0,len(records)):
    SeqIO.write(records[i][0], "component"+str(records[i][1])+".fasta", "fasta")