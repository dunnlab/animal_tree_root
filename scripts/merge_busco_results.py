#!/usr/bin/env python
import pandas as pd

metazoa_results_file = 'blast/busco/run_animal_root_metazoa_busco/full_table_animal_root_metazoa_busco.tsv'
busco_ids_file = 'blast/busco/busco_datasets/metazoa_33208_OrthoDB9_orthogroup_info.txt'
outfile = 'blast/graphs/busco_metazoa_results.tsv'

metazoa_results = pd.read_csv(metazoa_results_file, sep='\t', 
                 names = ['Busco_Id', 'Status', 'Sequence', 'Score', 'Length'],
                 comment='#', index_col=0)
busco_ids = pd.read_csv(busco_ids_file, sep='\t', index_col=0)

combined = metazoa_results.join(busco_ids)
combined.to_csv(outfile, sep='\t', index_label="Busco_Id")
