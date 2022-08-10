# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 10:44:27 2022

@author: matil
"""
import pandas as pd
from Bio import SeqIO
from compare_comp_nrs import add_co_nr,process_df
from glob import glob

def edit_df(df,partition_df):
    for i,name in zip(df.index.values.astype(int), df.name):
        df = add_co_nr(i,df,'comp_nr','name',name,partition_df)
        print('\rProcessing row:'+str(i)+'/'+str(len(df)),end="")
    return df

def read_blast_query(dirpattern,colname):
    l = []
    for i, infile in enumerate(sorted(glob(dirpattern))):
        with open(infile, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                l.append(record.description)
#        if i ==0:
#            l = l
#        else:
#            df = pd.concat((l, tmp_l),ignore_index=True)
    df = pd.DataFrame(l, columns={colname})
    return df

#imort fasta headers from blast queries to use as name
infile = r'..\reconciliation\blast\queries\dunn_philippe\*.fa'
colname = "fasta_header"
headers_df = read_blast_query(infile,colname)

#headers = []
#with open(infile, "r") as f:
#    for record in SeqIO.parse(f, "fasta"):
#        headers.append(record.description)
#headers_df = pd.DataFrame(headers, columns={"fasta_header"})

#add name column
headers_df = headers_df.join(headers_df.fasta_header.str.split(':| ', expand = True))[["fasta_header",0,3]]
headers_df["name"] = headers_df[0] + headers_df[3]
#import partition map 
partitions_df_expanded_r = pd.read_csv("../manuscript/part_map_glob.csv", sep="\t")
headers_df = headers_df.drop([0,3],axis=1)

#process partition map
partitions_df_expanded_r,group_df_r = process_df(partitions_df_expanded_r)

headers_df_processed = edit_df(headers_df.copy(),partitions_df_expanded_r)