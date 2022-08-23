
while read line;
do	
	perl /proj/metazoa_phylo/private/ShareScripts/Newick2Nexus.pl -i /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/mafft/test/$line/TEST_$line.contree -o /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/trees/nexus/$line.contree.nex
	#mv $line.contree.nex /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/trees/nexus/
	perl /proj/metazoa_phylo/private/ShareScripts/metazoa_phylo_TreeColors.pl -i /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/trees/nexus/$line.contree.nex -o /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/trees/colored/new/$line.colr.treefile -p /proj/metazoa_phylo/private/ShareScripts/taxon_table.ManuallyCorrected -c /proj/metazoa_phylo/private/ShareScripts/taxon_table.ColoringScheme

done < /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/name_list_paralogs.txt

