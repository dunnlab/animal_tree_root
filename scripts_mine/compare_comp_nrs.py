# -*- coding: utf-8 -*-
"""
Created on Tue Aug  2 10:15:14 2022

@author: matil
"""

#Idea of script: from blast/diamond result rename sequence ids to: 
#matrix (paper) + partition. Use component map (Table S9 from supplementary
#material AND partition_map_global from R) to find what component maps to each 
#sequence. Make a subset of sequences where at least one sequence doesnt map 
#to any component and  (if any) a subset om sequences that dont match component
#number


#!/usr/bin/env python
from glob import glob
# import matplotlib.pyplot as plt
import pandas as pd
#from goatools import obo_parser
#import sys
#import time

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

#get component number corresponding to name (paper+partition), when no match is found for a name the component nr is set to -1 
def set_co_nr(name,partition_df):
    co_nr= partition_df[partition_df.matrix_partition.isin([name])].component_number
    if len(co_nr) > 0:
        co_nr = int(co_nr)
    else:
        co_nr = -1
    return co_nr

#Add component number based on name as column in dataframe    
def add_co_nr(i,df,col_name,col_name2,name,partition_df):
    com_nr = set_co_nr(name,partition_df)
    df.at[i,col_name]=com_nr
    #df.at[i,col_name2]=name
    return df

#Add component numbers for both search and query in all rows of dataframe, input df must have columns with names qreqid and sseqid
def edit_df(df,partition_df):
    for i,qname,sname in zip(df.index.values.astype(int), df.qseqid,df.sseqid):
        df = add_co_nr(i,df,'scomp_nr','sname',sname,partition_df)
        df = add_co_nr(i,df,'qcomp_nr','qname',qname,partition_df)
        print('\rProcessing row:'+str(i)+'/'+str(len(df)),end="")
    return df

#return dataframe consisting of non matching entries (compares component ids)
def find_non_matches(df):
    non_match_comp_df = df[df['scomp_nr'] != df['qcomp_nr']]
    return non_match_comp_df

#return dataframe consisting of entries missing component number
def find_empty(df):
    empty_df = df[(df['scomp_nr'] == -1) & (df['qcomp_nr'] == -1)]
    return empty_df

#return dataframe consisting of entries missing component number
def find_one_empty(df):
    empty_df = df[(df['scomp_nr'] == -1) | (df['qcomp_nr'] == -1)]
    return empty_df

#process partition df to standardise naming and remove excess columns
def process_df(df):
    col_names = ["component_number","matrix","partition"]
    df = df[col_names]
    df_sep = df.copy()
    df_sep["matrix_partition"] = df[col_names[1]] + df[col_names[2]].astype(str)+","
    df["matrix_partition"] = df[col_names[1]] + df[col_names[2]].astype(str)
    df = df[[col_names[0],"matrix_partition"]] #expanded (original) file
    df_sep=df_sep[[col_names[0],"matrix_partition"]] #non expanded file, to be grouped
    group_df = df_sep.groupby([col_names[0]]).agg({'matrix_partition':'sum'}) #group, non expanded file
    return df,group_df

def main():
    print("Hello World!")
# Read in All vs All results (ava) from diaomond blast runs 
# cutoffs: chose to retain enough results to keep at least 1 hit for ~90% of partitions
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
#ava_df_unique = ava_df.drop_duplicates()

#load partition maps
supp_mat_df = pd.read_excel('../Li_etal_Supplementary_Tables.xlsx', sheet_name="Table S9",)
partitions_df_expanded_r = pd.read_csv("../files_mine/part_map_glob.csv", sep="\t")

#filter out hits with only one match
#supp_mat_df = supp_mat_df[supp_mat_df["nodes_in_component"] !=1] 


partitions_df_expanded,group_df_supp = process_df(supp_mat_df)
partitions_df_expanded_r,group_df_r = process_df(partitions_df_expanded_r)

#process all v all dataframe
ava_df_process = edit_df(ava_df.copy(),partitions_df_expanded)
ava_df_process_r = edit_df(ava_df.copy(),partitions_df_expanded_r)  #see if results differ when using the global map from r vs the supplementary material


#get non matching results
non_matching_df = find_non_matches(ava_df_process)
non_matching_df_unique = non_matching_df.drop_duplicates()  #21 rows, same as below
#get empty results
no_hits_df = find_empty(ava_df_process)
no_hits_df_unique = no_hits_df.drop_duplicates()
#get dataframe with at least one misisng comp nr
no_hits_df_one = find_one_empty(ava_df_process)
no_hits_df_one_unique = no_hits_df_one.drop_duplicates()    #33 rows, missing Moroz2014_3d0129	92.0, Moroz2014_3d0170	144.0, Moroz2014_3d0170	144.0


#get non matching results (r)
non_matching_df_r = find_non_matches(ava_df_process_r)
non_matching_df_unique_r = non_matching_df_r.drop_duplicates()  #18 rows
#get empty results
no_hits_df_r = find_empty(ava_df_process_r)
no_hits_df_unique_r = no_hits_df_r.drop_duplicates()
#get dataframe with at least one misisng comp nr
no_hits_df_one_r = find_one_empty(ava_df_process_r)
no_hits_df_one_unique_r = no_hits_df_one_r.drop_duplicates()    #30 rows

if __name__ == "__main__":
    main()