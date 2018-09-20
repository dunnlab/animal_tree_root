#!/usr/bin/env python3
# adapted from :
# https://stackoverflow.com/questions/16504238/attempting-to-obtain-taxonomic-information-from-biopython
import re
import os
import sys
import glob
import urllib
import pickle
from nexus import NexusReader # https://pypi.org/project/python-nexus/
from Bio import SeqIO
from Bio import Entrez
Entrez.email = 'b.evans@yale.edu'

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
    query_name = query_name.rstrip('_.')
    if query_name.endswith('_sp'):
        query_name[:-3]
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
    try:
        search = Entrez.efetch(id = tax_id, db = "taxonomy", retmode = "xml")
        return Entrez.read(search)
    except Exception as e:
        print("taxonomy lookup error: ", e)
    return [{'Lineage':'na'}]

def get_manual_entries():
    fh = open('manual_taxonomy_map.tsv')
    fh.readline() # discard header
    manual = dict()
    for line in fh:
        tmp = line.rstrip().split()
        manual[tmp[1]] = tmp[2]
    return manual

def get_pickles():
    """be kind to API, don't ask for what we already know"""
    if os.path.isfile('known_ids.pickle'):
        k_ids = pickle.load(open('known_ids.pickle', 'rb'))
    else:
        k_ids = {} 

    if os.path.isfile('known_tax.pickle'):
        k_tax = pickle.load(open('known_tax.pickle', 'rb'))
    else:
        k_tax = {}
    return k_ids, k_tax

if __name__ == '__main__':

    known_ids, known_tax = get_pickles()
    manual_entries = get_manual_entries()
    table_columns = ["original_matrix", "matrix_name", "ncbi_tax_id", 
                     "ncbi_taxonomy", "relabelled_name", "clade_assignment"]
    sep = '\t'
    table_out = open("taxon_table.tsv", 'w')
    print(sep.join(table_columns), file=table_out)

    # glob hack to match nex and phy
    for file_name in glob.glob('../considered_data/**/*.[pn][eh][xy]', recursive=True):
        names_list = []
        print("Reading file", file_name)
        try:
            names_list = get_names(file_name)
        except Exception as e:
            print("Error loading", file_name, ":\n", e)
        for original_name in names_list:
            manual_name = ''
            if original_name in manual_entries:
                manual_name = manual_entries[original_name]
                known_ids[original_name] = get_tax_id(manual_name)
            elif original_name not in known_ids:
                known_ids[original_name] = get_tax_id(original_name)

            known_tax[known_ids[original_name]] = get_tax_data(known_ids[original_name])[0]

            print(sep.join([ file_name, original_name, known_ids[original_name],
                                known_tax[known_ids[original_name]]['Lineage'],
                                manual_name, original_name]), file=table_out)

    table_out.close()
    pickle.dump(known_ids, open('known_ids.pickle', 'wb'))
    pickle.dump(known_tax, open('known_tax.pickle', 'wb'))
