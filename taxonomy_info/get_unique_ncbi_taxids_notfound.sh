awk '$4 == "nan" {n=split($1,A,"/"); print A[3],$6}' taxon_table.tsv | sort | uniq > find_our_taxids.txt

