# Animal root

A comparison of phylogenetic studies relevant to placing the root of the 
animal phylogeny.

## scripts directory

To run the python scripts in this directory, you'll need the `animal_root` conda environment. To recreate this environment, install [mini|ana]conda and run:
```
cd scripts
./create_conda_env.sh
```

The `nexus2phylip.py` script (which requires [python-nexus](https://pypi.org/project/python-nexus/)) accepts an input nexus and an output phylip file. It will write a .nex file of the input nexus partitions if they are found. 

The `generate_provenance_table.py` script is meant to be run from the `taxonomy_info` directory. See the `readme.md` there for more info.
