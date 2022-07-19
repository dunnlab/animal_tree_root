# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:09:23 2022

@author: matil
"""

import pandas

infile = r'C:\Users\matil\OneDrive\Dokument\GitHub\summer_project_2022\files_mine\gene_positions_Philippe_fixed'
data_fixed =pandas.read_csv(infile, sep=' ', header=None, names=["Gene_name", "Position"])

infile2 = r'C:\Users\matil\OneDrive\Dokument\GitHub\summer_project_2022\files_mine\gene_positions_Philippe_raw'
data_raw =pandas.read_csv(infile2, sep=' ', header=None, names=["Gene_name", "Position"])

data_non_overlap = data_fixed[~data_fixed.apply(tuple,1).isin(data_raw.apply(tuple,1))]