#!/usr/bin/env python
import os
import sys
import argparse
from shutil import copy2 as copy
from collections import defaultdict
from os.path import  abspath, basename, dirname, join, splitext
from nexus import NexusReader # https://pypi.org/project/python-nexus/
from Bio import AlignIO

def read_taxon_table(taxon_table_fn):
    names_map = defaultdict(lambda: {})
    
    table = open(taxon_table_fn, 'r')
    header = table.readline().split()
    original_file = header.index('original_matrix')
    old_name = header.index('matrix_name')
    new_name = header.index('relabelled_name')
    for line in table:
        tmp = line.rstrip('\n').split('\t')
        manuscript = basename(dirname(tmp[original_file]))
        names_map[manuscript][tmp[old_name]] = tmp[new_name]
    return names_map

def process_nexus(manuscript_name, names_map, in_filename, out_phylip_fn, out_partition_fn):
    n = NexusReader(in_filename)
    out_phylip = open(out_phylip_fn, 'w')

    # write phylip header
    out_phylip.write('{} {}\n'.format(n.data.ntaxa, n.data.nchar))

    # write character matrix
    for taxon, characters in n.data:
        if taxon in names_map[manuscript_name] and names_map[manuscript_name][taxon] != 'na':
            out_phylip.write(names_map[manuscript_name][taxon]+' ')
        else:
            out_phylip.write(taxon+' ')
        if 'missing' in n.data.format:
            missing =  n.data.format['missing']
            for i, character in enumerate(characters):
                if character == missing:
                    characters[i] = '-'
        out_phylip.write(''.join(characters)+'\n')
    out_partition = open(out_partition_fn, 'w')
    out_partition.write('\n'.join(n.sets.block)+'\n')

def process_phylip(manuscript_name, names_map, in_filename, out_phylip_fn, out_partition_fn):
    
    try:
        in_partition_fn = splitext(in_filename)[0] + '.nex'
        copy(in_partition_fn, out_partition_fn)
    except IOError as e:
        print("Error, couldn't find partition for file: {}".format(in_filename))
        sys.exit(1)

    phylip_out = open(out_phylip_fn,'w')
    alignment = AlignIO.read(in_filename, 'phylip-relaxed')
    for taxon in alignment:
        if taxon.id in names_map[manuscript_name] and names_map[manuscript_name][taxon.id] != 'na':
                taxon.id = names_map[manuscript_name][taxon.id]
    AlignIO.write(alignment, out_phylip_fn, 'phylip-relaxed')


desc = """prepare_phylip.py
Simple script to take a nexus or phylip format file from the manuscripts we are considering for this study and
1) convert it to a standard phylip + nexus partitions format
2) rename original taxa to conform to a global namespace in the process

"""

parser = argparse.ArgumentParser(description=desc,
                                 usage='%(prog)s --in-file file.nex --out-dir out_dir/ --taxon-table path/to/taxon_table.tsv', 
                                 formatter_class=argparse.RawTextHelpFormatter,
                                 prog='prepare_phylip.py')
parser.add_argument('-i', '--in-file',
                    required=True,
                    help='nexus or phylip file. if phylip, also expect a .nex file with same prefix where partitions are')
parser.add_argument('-o', '--out-dir',
                    required=True,
                    help='Directory to put converted data files')
parser.add_argument('-t', '--taxa-table',
                    required=True,
                    help='Table of taxon names and metadata. To be used in renaming samples.')
args = parser.parse_args()

in_filename = args.in_file
out_prefix = args.out_dir
dataset_name = splitext(basename(in_filename))[0]
manuscript_name = basename(dirname(abspath(in_filename)))

out_phylip_fn =    join(out_prefix, '{}_{}.phy'.format(manuscript_name, dataset_name))
out_partition_fn = join(out_prefix, '{}_{}.nex'.format(manuscript_name, dataset_name))
print('processing {}'.format(in_filename))
names_map = read_taxon_table(args.taxa_table)
if in_filename.endswith('.nex'):
    process_nexus(manuscript_name, names_map, in_filename, out_phylip_fn, out_partition_fn)
elif in_filename.endswith('.phy'):
    process_phylip(manuscript_name, names_map, in_filename, out_phylip_fn, out_partition_fn)

