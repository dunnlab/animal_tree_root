# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 11:39:34 2022

@author: matil
"""

"""
Pickle2JSON is a simple Python Command Line program for converting Pickle file to JSON file.
Arguments: Only one (1) argument is expected which is the pickle file.
Usage: python pickle2json.py myfile.pkl
Output: The output is a JSON file bearing the same filename containing the JSON document of the converted Pickle file.
"""

# import libraries
import pickle
import json
import sys
import os
import pandas as pd

inf = r"C:/Users/matil/OneDrive/Dokument/GitHub/summer_project_2022/reconciliation/taxonomy_info/known_ids.pickle"
#outfile =r"C:/Users/matil/OneDrive/Dokument/GitHub/summer_project_2022/reconciliation/taxonomy_info/known_tax.json"
# open pickle file
with open(inf, 'rb') as infile:
    objects = pickle.load(infile)

# convert pickle object to json object
json_obj = json.loads(json.dumps(obj, default=str))

# write the json file
with open(
        os.path.splitext(inf)[0] + '.json',
        'w',
        encoding='utf-8'
    ) as outfile:
    json.dump(json_obj, outfile, ensure_ascii=False, indent=4)