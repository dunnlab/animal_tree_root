#!/bin/bash

module load bioinfo-tools
module load augustus/3.3.3-CGP
module load BUSCO

run_BUSCO.py -c 20 -o busco_results/animal_root_metazoa_busco -m prot -l $BUSCO_LINEAGE_SETS/metazoa_odb10 -i /crex/proj/metazoa_phylo/private/animal_tree_root_fork/add_new_data_test/db_test/animal_root_test.fa
run_BUSCO.py -c 20 -o busco_results/animal_root_eukaryota_busco -m prot -l $BUSCO_LINEAGE_SETS/eukaryota_odb10 -i /crex/proj/metazoa_phylo/private/animal_tree_root_fork/add_new_data_test/db_test/animal_root_test.fa        

