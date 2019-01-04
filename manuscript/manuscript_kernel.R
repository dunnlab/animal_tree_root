# The computationally intensive analyses for the manuscript 
# presented at https://github.com/caseywdunn/animal_root/tree/master/manuscript .
#
# Executing this code generates manuscript.RData, which contains analysis 
# results. That file is then read by manuscript.rmd for rendering and 
# presentation of the results.
# 
# The code presented here is roughly in the order of the analyses presented 
# in the manuscript, though there are exceptions

## Preliminaries

time_start = Sys.time()

library( tidyverse )
library( magrittr )
library( viridis )
library( ape )
library( igraph )
library( hutan )
library( doParallel )

source( "functions.R" )

# Set system computational parameters
cores = detectCores() - 1
if ( cores < 1 ) {
	cores = 1
}

# Register parallel workers for %dopar%
registerDoParallel( cores )

# Set up constants
focal_matrices = 
	read.table(text = 
		"manuscript matrix
		Borowiec2015	Total1080
		Chang2015 Chang2015
		Dunn2008  Dunn2008
		Hejnol2009 Hejnol2009
		Moroz2014 ED3a
		Nosenko2013 nonribosomal_9187_smatrix
		Nosenko2013 ribosomal_11057_smatrix
		Philippe2009	Philippe2009
		Ryan2013	est.opisthokonta
		Ryan2013  genome.opisthokonta
		Simion2017	supermatrix_97sp_401632pos_1719genes
		Whelan2015	Metazoa_Choano_RCFV_strict
		", 
		header = TRUE, 
		stringsAsFactors = FALSE)

bootstrap_threshold = 90
posterior_prob_threshold = 95

# Load data

papers = read_tsv( "../data_processed/tables/previously_published_manuscripts.tsv" )
datasets = read_tsv( "../data_processed/tables/previously_published_matrices.tsv" )
analyses = 
	read_tsv( "../data_processed/tables/previously_published_analyses.tsv", col_types = 'iccccccccccclnnnncc' ) %>%
	mutate( clade = factor( clade, levels= c( 'Choanimalia', 'Holozoa', 'Opisthokonta' ) ) )

taxonomy_reference = read_tsv("../reconciliation/taxonomy_info/taxon_table.tsv")

partition_map_global = 
	read_tsv("../reconciliation/blast/graphs/diamond_partition_clusters.tsv") %>%
	dplyr::rename(partition=partition_name) %>%
	dplyr::rename(gene=cluster_number) %>%
	mutate( gene = as.character(gene) )

analyses$result = "Unresolved"
analyses$result[ ( analyses$inference == "Bayesian" ) & (analyses$support_porifera_sister >= posterior_prob_threshold) ] = "Porifera-sister"
analyses$result[ ( analyses$inference == "Bayesian" ) & (analyses$support_ctenophora_sister >= posterior_prob_threshold) ] = "Ctenophora-sister"
analyses$result[ ( analyses$inference == "ML" ) &       (analyses$support_porifera_sister >= bootstrap_threshold) ] = "Porifera-sister"
analyses$result[ ( analyses$inference == "ML" ) &       (analyses$support_ctenophora_sister >= bootstrap_threshold) ] = "Ctenophora-sister"
analyses$result = factor( analyses$result )

analyses$model_summary = analyses$model_rate_matrix
long_summary = paste(analyses$model_rate_matrix, analyses$model_equilibrium, sep="+")
analyses$model_summary[ ! is.na(analyses$model_equilibrium) ] = long_summary[ ! is.na(analyses$model_equilibrium) ]
analyses$model_summary[ grepl("WAG", analyses$model_summary) ] = "WAG" # Simplify WAG
analyses$model_summary = factor( analyses$model_summary, levels=c("GTR+CAT", "F81+CAT", "LG", "WAG") )


# Matrix taxon composition

clades = c( "Fungi", "Ichthyosporea", "Filasterea", "Choanoflagellida", "Ctenophora", "Porifera", "Placozoa", "Bilateria", "Cnidaria" )

taxa = 
	taxonomy_reference %>% 
	distinct( relabelled_name, clade_assignment, ncbi_tax_id ) %>%
	dplyr::rename( taxon=relabelled_name, clade=clade_assignment ) %>%
	mutate( clade = factor( clade, levels=clades )  )

matrix_path = "../data_processed/matrices/"
phylip_file_names = list.files(path = matrix_path, pattern = ".+\\.phy$")

parse_phylip =	function( phylip_file ){
	print(phylip_file)
	partition_file = sub( "phy$", "nex", phylip_file )
	
	phylip_path = paste(c(matrix_path, phylip_file), collapse="")
	partition_path = paste(c(matrix_path, partition_file), collapse="")
	
	base_name = sub( "\\.phy$", "", phylip_file )
	elements = strsplit(base_name, "_") %>% unlist()
	manuscript_name = elements[1]
	matrix_name = base_name
	
	# Do some cludgey renaming to make it easy to step into function below for debugging
	alignment_path = phylip_path
	partition_file = partition_path
	taxon_map_global = taxa
	
	sequence_matrix = 
		PartitionedMultipleAlignment( 
			alignment_path, 
			partition_file, 
			taxon_map_global, 
			partition_map_global=partition_map_global, 
			manuscript_name=manuscript_name, 
			matrix_name=matrix_name
		)
	
	return( sequence_matrix )
}

sequence_matrices = foreach( phylip_file=phylip_file_names) %dopar%
	parse_phylip( phylip_file )

# Matrix gene composition

busco_results = 
	read_tsv("../reconciliation/blast/graphs/busco_metazoa_results.tsv") %>% 
	filter( Status != "Missing" ) 

sequences_all_txt = 
	read_tsv("../reconciliation/blast/all_parts_list.txt")

names(sequences_all_txt) = c("full_name")
sequences_all = str_split_fixed( sequences_all_txt$full_name, ":", 4 ) %>% 
	as_tibble()
names( sequences_all ) = c( "matrix", "species", "ncbi_taxon_id", "partition" )

partitions_all = 
	sequences_all %>% 
	group_by( matrix, partition ) %>% 
	summarise( n_sequences =n() )


# Multiple fields are in one colon delimited string. Need to parse them out.
# Busco result example
# "Moroz2014:ED3a:Capitella:51293:0241"  
# manuscript:matrix:species:NCBI_taxon_id:partition

Bs = 
	str_split_fixed( busco_results$Sequence, ":", 4 ) %>% 
	as_tibble()

names( Bs ) = c( "matrix", "species", "ncbi_taxon_id", "partition" )

Bs %<>% mutate( ncbi_taxon_id = as.integer(ncbi_taxon_id) )

busco_results %<>% bind_cols( Bs )

# Find overrepresented partitions, ie partitions with genes that hit more than one Busco
busco_distinct =  busco_results %>% select( matrix, partition, Description ) %>% distinct(  )
busco_overrepresented = busco_distinct %>% group_by( matrix, partition ) %>% summarise( n = n() ) %>% filter( n > 1 )

busco_overrepresented_full = left_join( busco_overrepresented, busco_results )

# Combine and summarize results
partition_to_busco_map = 
	busco_distinct %>% 
	group_by( matrix, partition ) %>% 
	summarise( BUSCO=names(which(table(Description) == max(table(Description)))[1]) )

partition_map_global %<>% left_join( partition_to_busco_map, by=c("matrix", "partition") )


busco_summary = 
	partition_to_busco_map %>% 
	group_by( matrix ) %>%
	summarise( n_busco_partitions = n() )

matrix_summary =
	partitions_all %>%
	group_by( matrix ) %>%
	summarise( n_partitions = n() ) %>%
	left_join(
		busco_summary,
		by=c("matrix")
	)

matrix_summary$manuscript = str_split( matrix_summary$matrix, "_", simplify=TRUE )[,1]

cluster_summary =
	partition_map_global %>% 
	group_by( gene ) %>% 
	summarise( 
		n_partitions = n(), 
		n_matrices= length(unique(matrix)),
		n_with_busco = sum(! is.na(BUSCO)),
		n_unique_busco = length(unique(na.omit(BUSCO)))
	)


# Matrix overlap
matrix_overlap = 
	lapply(sequence_matrices, function(x) lapply(sequence_matrices, function(y) compute_matrix_overlap(x,y))) %>% 
	unlist(recursive=FALSE) %>% 
	bind_rows()

# Remove reciprocal comparisons and comparisons to self
n = nrow(matrix_overlap)
mask = lower.tri(matrix(1:n^2, n, n))
dim(mask) = NULL
matrix_overlap = matrix_overlap[mask,]


# New analyses of published matrices

trees_path = "../trees_new/"
contree_extension = "treefile"
tree_file_names = list.files(path = trees_path, pattern = str_c( "\\.", contree_extension, "$" ))

# Temp fix re issue #8
tree_file_names = tree_file_names[ ! grepl("Philippe2009", tree_file_names) ]

tree_file_names = tree_file_names[ ! grepl("Nosenko2013", tree_file_names) ]

parse_tree = 	function( tree_file ){
	print(tree_file)
	log_file = sub( str_c( contree_extension, "$" ), "log", tree_file )
	bootstrap_file = sub( str_c( contree_extension, "$" ), "ufboot", tree_file )
	
	tree_file_path = paste(c(trees_path, tree_file), collapse="")
	log_file_path = paste(c(trees_path, log_file), collapse="")
	
	if(!file.exists(log_file_path)){
		warning( str_c( "Missing log file: ", log_file_path ) )
		return( NA )
	}
	
	bootstrap_file_path = paste(c(trees_path, bootstrap_file), collapse="")
	if(!file.exists(bootstrap_file_path)){
		warning( str_c( "Missing bootstrap file: ", bootstrap_file_path ) )
		return( NA )
	}
	
	base_name = sub( str_c( "\\.", contree_extension, "$" ), "", tree_file )
	elements = strsplit(base_name, "\\.") %>% unlist()
	matrix_name = elements[1]
	model_name = elements[2]
	
	tree = read.tree( file = tree_file_path )
	
	# Get the clades of the tips
	clades = 
		taxonomy_reference$clade_assignment[ 
			match( tree$tip.label, taxonomy_reference$relabelled_name ) 
			]
	
	# Get the outgroup sampling
	
	outgroup = NA
	
	sampling = "Metazoa"
	if( "Choanoflagellida" %in% clades ){
		sampling = "Choanimalia"
		outgroup = which( clades == "Choanoflagellida" )
	}
	
	if( ("Ichthyosporea" %in% clades) | ("Filasterea" %in% clades) ){
		sampling = "Holozoa"
		outgroup = c( which( clades == "Ichthyosporea" ),  which( clades == "Filasterea" ))
	}
	
	if( "Fungi" %in% clades ){
		sampling = "Opisthokonta"
		outgroup = which( clades == "Fungi" )
	}
	
	if( is.na(outgroup) ){
		stop( "There are no outgroups in the tree, can't test rooting topology" )
	}
	
	# Reroot the tree
	# Use a single taxon from the outgroup since the outgroup is not guaranteed to be monophyletic
	tree = root( tree, outgroup[1] )
	
	# Get the clades of the tips again, since tip order has changed
	clades = 
		taxonomy_reference$clade_assignment[ 
			match( tree$tip.label, taxonomy_reference$relabelled_name ) 
			]
	
	if( sum(is.na(clades)) > 0 ){
		missing_clades = tree$tip.label[ is.na(clades) ]
		stop( paste( c( "Some tips have missing clade names:", missing_clades ), collapse=" " ))
	}
	
	if( sum(clades == "Ctenophora") < 1 ){
		stop( "There are no ctenophores in the tree, can't test rooting topology" )
	} 
	
	if( sum(clades == "Porifera") < 1 ){
		stop( "There are no sponges in the tree, can't test rooting topology" )
	} 
	
	# Get the set of taxa in the group Ctenophora+Placozoa+Bilateria+Cnidaria, the clade that exists under Proifera-sister
	CtPlBiCn = tree$tip.label[ clades %in% c("Ctenophora", "Placozoa", "Bilateria", "Cnidaria") ]
	
	# Get the set of taxa in the group Porifera+Placozoa+Bilateria+Cnidaria, the clade that exists under Ctenophora-sister
	PoPlBiCn = tree$tip.label[ clades %in% c("Porifera", "Placozoa", "Bilateria", "Cnidaria") ]
	
	# Read the bootstrap trees
	bootstrap_trees = read.tree(file = bootstrap_file_path)
	
	Ctenophora_sister = lapply(
		bootstrap_trees,
		function( phy ){
			is_monophyletic( phy, PoPlBiCn )
		}
	) %>% 
		unlist() %>%
		mean()
	
	Porifera_sister = lapply(
		bootstrap_trees,
		function( phy ){
			is_monophyletic( phy, CtPlBiCn )
		}
	) %>% 
		unlist() %>%
		mean()				
	
	# Stuff a few other things into the tree object
	tree$matrix = matrix_name
	tree$model = model_name
	tree$modelfinder = FALSE
	tree$clades = clades
	tree$sampling = sampling
	tree$ctenophora_sister = Ctenophora_sister
	tree$porifera_sister = Porifera_sister
	
	# parse model
	log_lines = read_lines( log_file_path )
	
	best_line = log_lines[ grepl("Best-fit model:", log_lines) ]
	
	if( length(best_line) > 1 ){
		warning("More than one best model found")
	}
	
	if( length(best_line) > 0 ){
		tree$modelfinder = TRUE
		elements = strsplit( best_line[1], " " ) %>% unlist()
		tree$model = elements[3]
	}
	
	return( tree )
}


trees = foreach( tree_file=tree_file_names) %dopar%
	parse_tree( tree_file )

tree_summary = lapply(
	trees,
	function( tree ){
		data.frame(
			matrix = tree$matrix,
			model = tree$model,
			modelfinder = tree$modelfinder,
			sampling = tree$sampling,
			stringsAsFactors = FALSE
		)
	}
) %>% 
bind_rows()


# Partition comparison across matrices

Dpart = 
	read_tsv("../reconciliation/blast/graphs/partitions_graph.tsv") %>% 
	mutate( node1 = paste(matrix1, part1, sep = ':') ) %>% 
	mutate( node2 = paste(matrix2, part2, sep = ':') )

# Hits to self not included anymore, so skip this bit
# Dpart_self = Dpart %>% 
#	 filter( node1==node2 ) %>%
#	 select( matrix1, part1, node1, count )
# ggplot( Dpart_self ) + geom_density(aes(count))


# Retain only nonself edges
Dpart %<>% filter( node1!=node2 )


# Network analysis
# Builds on https://www.jessesadler.com/post/network-analysis-with-r/

D_source = Dpart %>%
	distinct(node1, .keep_all = TRUE ) %>%
	select( node1, matrix1, part1 ) %>%
	dplyr::rename(label = node1, matrix = matrix1, part = part1 )

D_dest = Dpart %>%
	distinct(node2, .keep_all = TRUE ) %>%
	select( node2, matrix2, part2 ) %>%
	dplyr::rename(label = node2, matrix = matrix2, part = part2 )

nodes = 
	bind_rows( D_source, D_dest ) %>% 
	distinct( label, .keep_all = TRUE ) %>%
	arrange( label ) %>% 
	rowid_to_column("id")


sources = Dpart %>%
	distinct(node1) %>%
	dplyr::rename(label = node1)

destinations = Dpart %>%
	distinct(node2) %>%
	dplyr::rename(label = node2)

nodes2 = 
	full_join(sources, destinations, by = "label") %>% 
	arrange( label ) %>% 
	rowid_to_column("id")

edges = Dpart %>% 
	select( node1, node2, count ) %>%
	left_join(nodes, by = c("node1" = "label")) %>% 
	dplyr::rename(from = id) %>% 
	left_join(nodes, by = c("node2" = "label")) %>% 
	dplyr::rename(to = id) %>%
	dplyr::rename(weight = count) %>%
	select( from, to, weight )


gene_igraph = graph_from_data_frame(d = edges, vertices = nodes, directed = FALSE)
gene_clusters = clusters(gene_igraph, mode=c("weak", "strong"))

nodes$cluster = gene_clusters$membership



## Record information about the session
session_info_kernel = sessionInfo()
system_time_kernel = Sys.time()

commit_kernel = 
	system("git log | head -n 1", intern=TRUE) %>% 
	str_replace("commit ", "")

time_stop = Sys.time()
time_run = time_stop - time_start

## Write the results to prepare them for manuscript.rmd

save.image("manuscript.RData")