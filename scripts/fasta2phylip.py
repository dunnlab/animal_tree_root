#!/usr/bin/env python
import sys
from Bio import AlignIO

AlignIO.write(AlignIO.parse(sys.argv[1], 'fasta'), sys.argv[2], 'phylip-relaxed')
