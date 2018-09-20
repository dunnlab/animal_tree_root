## `taxonomy_info`
We want to track the sample names from each respective manuscript globally to ensure we are consistant and transparent in how the matrices compared. Here is where we try to get all the sample names sorted out by searching for them at NCBI.

### Taxon table
The table is generated with the python script `scripts/generate_provenance_table.py`. You'll need to be using the `animal_root` conda environment to run it.

To rerun the taxonomy table generation, run:

```
source activate animal_root
cd taxonomy_info/
# If you want to start from scratch, remove cached results first:
# rm *.pickle taxon_table.tsv provenance_table.log
./gen_prov_table.sh
```
The script looks in the `considered_data` directory for `.phy` and `.nex` files, and stderr and stdout from the script running gets saved to `provenance_table.log`. It will try looking up each name from `manual_taxonomy_map.tsv` and the cached names in the pickles saved in this directory befor making queries to NCBI. The taxonomy table gets saved to `taxon_table.tsv`. 

`taxon_table.tsv` is A tab-separated table where each row corresponds to a sample in the matrices that got read in. NCBI ID and Lineage information are `nan` and `na` respectively when the name -> ID search fails. The output table has the following columns gets generated:

| Column Name | Description |
|---|---|
| `original_matrix` | name of the original matrix file |
| `matrix_name` | name from the original matrix |
| `ncbi_tax_id` | NCBI Taxonomy ID of the  |
| `ncbi_taxonomy` | All the lineage info for this sample |
| `relabelled_name` | If the sample was manually assigned a name, this is it |
| `clade_assignment` | dummy field right now, currently matches `matrix_name` |

### Taxa that need a human eye

Running the script `get_unique_ncbi_taxids_notfound.sh ` will generate a file called `find_our_taxids.txt` with two columns: The manuscript of origin, and the name that couldn't be found with the `generate_provenance_table.py` script. These will have to be looked up / decided on manually and added to the taxon_table.tsv.