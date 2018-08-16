awk '$3 == "nan" {n=split($1,A,"/"); print A[3],$2}' taxon_table.tsv | sort | uniq > find_our_taxids.txt

