#!/usr/bin/env python3
# adapted from :
# https://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython
import re
import os
import sys
import glob
import pickle
from nexus import NexusReader # https://pypi.org/project/python-nexus/
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'casey.dunn@yale.edu'

def get_names(file_name):
    """read in a phy or nex file, return the taxa names found"""
    names = []
    if file_name.endswith('.phy'):
        phy = SeqIO.parse(file_name, 'phylip-relaxed')
        names = [t.name for t in phy]
    elif file_name.endswith('.nex'):
        nex = NexusReader(file_name)
        try:
            names = [t for t in nex.taxa]
        except AttributeError as e:
            # when there isn't a separate taxa block, get names from sequences
            names = [block[0] for block in nex.blocks['data']]
    return names

def get_tax_id(query_name):
    """sanitize query name, search for and return NCBI ID"""
    query_term = re.sub('[ _-]', '+', query_name)
    search = Entrez.esearch(term = query_term, db = "taxonomy", retmode = "xml")
    record = Entrez.read(search)
    try:
        tax_id = record['IdList'][0]
    except IndexError as e:
        print("Query", query_name, "not found.")
        tax_id = 'nan'
    return tax_id

def get_tax_data(tax_id):
    """Fetch taxonomy info with NCBI ID"""
    if tax_id == 'nan':
        return [{'Lineage':'na'}]
    search = Entrez.efetch(id = tax_id, db = "taxonomy", retmode = "xml")
    return Entrez.read(search)

# be kind to API, don't ask for what we already know
if os.path.isfile('known_ids.pickle'):
    known_ids = pickle.load(open('known_ids.pickle'))
else:
    known_ids = {} 

if os.path.isfile('known_tax.pickle'):
    known_tax = pickle.load(open('known_tax.pickle'))
else:
    known_tax = {} 

table_columns = ["original_matrix", "matrix_name", "ncbi_tax_id", 
                 "ncbi_taxonomy", "relabelled_name", "clade_assignment"]
table_rows = []

# glob hack to match nex and phy
for file_name in glob.glob('../raw_data/**/*.[pn][eh][xy]', recursive=True):
    names_list = []
    print("Reading file", file_name)
    try:
        names_list = get_names(file_name)
    except Exception as e:
        print("Error loading", file_name, ":\n", e)
    for name in names_list:
        if name not in known_ids:
            known_ids[name] = get_tax_id(name) 
            known_tax[known_ids[name]] = get_tax_data(known_ids[name])[0]

        table_rows.append([ file_name, name, known_ids[name],
                            known_tax[known_ids[name]]['Lineage'],
                            name, name])

pickle.dump(known_ids, open('known_ids.pickle', 'wb'))
pickle.dump(known_tax, open('known_tax.pickle', 'wb'))

# print to table file
sep = '\t'
table_out = open("taxon_table.tsv", 'w')
print(sep.join(table_columns), file=table_out)
for row in table_rows:
    print(sep.join(row), file=table_out)
