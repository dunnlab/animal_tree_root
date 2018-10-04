#!/usr/bin/env python
import sys
import os
from nexus import NexusReader # https://pypi.org/project/python-nexus/

in_filename = sys.argv[1]
out_phylip = sys.argv[2]

n = NexusReader(in_filename)
out = open(out_phylip, 'w')

# write phylip header
out.write('{} {}\n'.format(n.data.ntaxa, n.data.nchar))

# write character matrix
for taxon, characters in n.data:
    out.write(taxon+' ')
    if 'missing' in n.data.format:
        missing =  n.data.format['missing']
        for i, character in enumerate(characters):
            if character == missing:
                characters[i] = '-'
    out.write(''.join(characters)+'\n')
try:
    partition_fn = os.path.splitext(out_phylip)[0]+'.nex'
    partition_out = open(partition_fn, 'w')
    partition_out.write('\n'.join(n.sets.block)+'\nend;\n')
except Exception as e:
    if os.path.exists(filename):
        partition_out.close()
        os.remove(filename)

