# -*- coding: utf-8 -*-
"""
Created on Tue Aug  9 14:16:08 2022

@author: matil
"""
import pandas as pd

infile = r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\dosReisDonoghueYang_2015\dosReis_etal_superalignment.nex'

partitions = pd.read_csv(infile, sep="\t")
partitions['new'] = "CHARSET"
partitions['new1'] = partitions.protein
partitions['new2'] = "="
partitions['new3'] = partitions.start.astype(str)+'-'+partitions.end.astype(str)+';'

partitions_new = partitions[['new','new1','new2','new3']]

partitions_new.loc[0] = ['BEGIN SETS;',None, None, None]
partitions_new.loc[partitions_new.shape[0]] = ['END SETS;',None, None, None]
                      
partitions_new.to_csv(r'C:\Users\matil\OneDrive\Dokument\GitHub\animal_tree_root_fork\dosReisDonoghueYang_2015\dosReis_etal_superalignment_new.nex',sep="\t",index=False,header=False,)