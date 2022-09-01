
while read line;
do	
	perl /proj/metazoa_phylo/private/ShareScripts/Newick2Nexus.pl -i /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine//LG+C60+F_results/iqtree/$line/*.contree -o /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine//LG+C60+F_results/nexus/$line.contree.nex
	perl /proj/metazoa_phylo/private/ShareScripts/metazoa_phylo_TreeColors.pl -i /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine//LG+C60+F_results/nexus/$line.contree.nex -o /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine//LG+C60+F_results/colored_trees/$line.colr.treefile -p /proj/metazoa_phylo/private/ShareScripts/taxon_table.ManuallyCorrected -c /proj/metazoa_phylo/private/ShareScripts/taxon_table.ColoringScheme

done < /proj/metazoa_phylo/private/animal_tree_root_fork/files_mine/LG+C60+F_results/iqtree/successful_runs.txt
