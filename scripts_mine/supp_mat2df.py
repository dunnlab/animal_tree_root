# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 12:29:58 2022

@author: matil
"""
import pandas as pd
#Idea: 
#take li et all supp s9 and turn into df with comp nr and geneid
#output csv files with partitioning, one expanded and one where sequences 
#are grouped, expanded one is used in compare_comp_nrs

#!FUNCTION ALSO IN compare_comp_nrs.py

#import and transform dataframe
supp_mat_df = pd.read_excel('../Li_etal_Supplementary_Tables.xlsx', sheet_name="Table S9",)

#filter out hits with only one match
#supp_mat_df = supp_mat_df[supp_mat_df["nodes_in_component"] !=1] 

col_names = ["component_number","matrix","partition"]

def process_df(df,col_names):
    df = df[col_names]
    df_sep = df.copy()
    df_sep["matrix_partition"] = df[col_names[1]] + df[col_names[2]].astype(str)+","
    df["matrix_partition"] = df[col_names[1]] + df[col_names[2]].astype(str)
    df = df[[col_names[0],"matrix_partition"]] #expanded (original) file
    df_sep=df_sep[[col_names[0],"matrix_partition"]] #non expanded file, to be grouped
    group_df = df_sep.groupby([col_names[0]]).agg({'matrix_partition':'sum'}) #group, non expanded file
    return df,group_df

supp_mat_df,group_df = process_df(supp_mat_df,col_names)

#save files when filter used
#group_df.to_csv("../files_mine/component_partition_larger_1.csv", sep="\t")
#supp_mat_df.to_csv("../files_mine/component_partition_larger_than_1_expanded.csv", sep="\t")

#save files when filter not used
group_df.to_csv("../files_mine/component_partition.csv", sep="\t")
supp_mat_df.to_csv("../files_mine/component_partition_expanded.csv", sep="\t")