# Animal root

A comparison of phylogenetic studies relevant to placing the root of the 
animal phylogeny.

## scripts directory

The `nexus2phylip.py` script (which requires [python-nexus](https://pypi.org/project/python-nexus/)) accepts an input nexus and an output phylip file. It will write a .nex file of the input nexus partitions if they are found. 

The `generate_provenance_table.py` script (which requires [python-nexus](https://pypi.org/project/python-nexus/) and [BioPython](https://biopython.org/)) is meant to be run from the `taxonomy_info` directory. It looks recursively at all `.phy` and `.nex` files in the `raw_data` directory, parses out taxa where it can, and queries NCBI for taxonomy info. NCBI ID and Lineage information are `nan` and `na` respectively when the name -> ID search fails.

