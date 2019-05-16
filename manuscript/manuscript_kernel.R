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
library( igraph )
library( doParallel )

# setwd("/animal_root/manuscript")

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
analyses_published = 
	read_tsv( "../data_processed/tables/previously_published_analyses.tsv", col_types = 'iccccccccccclnnnnccc' ) %>%
	mutate( clade = factor( clade, levels= c( 'Choanimalia', 'Holozoa', 'Opisthokonta' ) ) )

taxonomy_reference = read_tsv("../reconciliation/taxonomy_info/taxon_table.tsv")

# raw columns:
# component_number, matrix, partition_name, edges, nodes_in_component, component_density, BUSCO_ID, BUSCO_description, 
# SwissProt_accession, SwissProt_description, GO_annotations, ribo_found
partition_map_global = 
	read_tsv("../reconciliation/blast/graphs/partition_components_split_annotated.tsv") %>%
	dplyr::rename(partition=partition_name) %>%
	dplyr::rename(gene=component_number) %>%
	mutate( gene = as.character(gene) )

analyses_published$result = "Unresolved"
analyses_published$result[ ( analyses_published$inference == "Bayesian" ) & (analyses_published$support_porifera_sister >= posterior_prob_threshold) ] = "Porifera-sister"
analyses_published$result[ ( analyses_published$inference == "Bayesian" ) & (analyses_published$support_ctenophora_sister >= posterior_prob_threshold) ] = "Ctenophora-sister"
analyses_published$result[ ( analyses_published$inference == "ML" ) &       (analyses_published$support_porifera_sister >= bootstrap_threshold) ] = "Porifera-sister"
analyses_published$result[ ( analyses_published$inference == "ML" ) &       (analyses_published$support_ctenophora_sister >= bootstrap_threshold) ] = "Ctenophora-sister"
analyses_published$result = factor( analyses_published$result )

analyses_published$model_summary = analyses_published$model_rate_matrix
long_summary = paste(analyses_published$model_rate_matrix, analyses_published$model_equilibrium, sep="+")
analyses_published$model_summary[ ! is.na(analyses_published$model_equilibrium) ] = long_summary[ ! is.na(analyses_published$model_equilibrium) ]
analyses_published$model_summary[ grepl("WAG", analyses_published$model_summary) ] = "WAG" # Simplify WAG
analyses_published$model_summary = factor( analyses_published$model_summary, levels=c("GTR+CAT", "F81+CAT", "LG", "WAG") )


# Matrix taxon composition

clades = c( "Fungi", "Ichthyosporea", "Filasterea", "Choanoflagellida", "Ctenophora", "Porifera", "Placozoa", "Bilateria", "Cnidaria" )

taxa = 
	taxonomy_reference %>% 
	distinct( relabelled_name, clade_assignment, ncbi_tax_id ) %>%
	dplyr::rename( taxon=relabelled_name, clade=clade_assignment ) %>%
	mutate( clade = factor( clade, levels=clades )  )

matrix_path = "../data_processed/matrices/"
phylip_file_names = list.files(path = matrix_path, pattern = ".+\\.phy$", full.names=TRUE)

parse_phylip =	function( phylip_path ){
	
	partition_path = sub( "phy$", "nex", phylip_path )
	base_name = sub( "\\.phy$", "", basename( phylip_path ) )
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

# Make contraint trees for each matrix
constraint_tree_path = "../trees_new/constraint_trees/"

generate_constraint_trees = function(seq_matrix){
	
	clades = unique(seq_matrix@taxon_map)
	
	# Get the set of taxa in the group Ctenophora+Placozoa+Bilateria+Cnidaria, the clade that exists under Proifera-sister
	CtPlBiCn = rownames(seq_matrix)[ seq_matrix@taxon_map %in% c("Ctenophora", "Placozoa", "Bilateria", "Cnidaria") ]
	
	# Get the set of taxa in the group Porifera+Placozoa+Bilateria+Cnidaria, the clade that exists under Ctenophora-sister
	PoPlBiCn = rownames(seq_matrix)[ seq_matrix@taxon_map %in% c("Porifera", "Placozoa", "Bilateria", "Cnidaria") ]
	
	# Write constraint trees for hypothesis testing
	CtPlBiCn_not = rownames(seq_matrix)[ ! rownames(seq_matrix) %in% CtPlBiCn ]
	porifera_sister_constraint_tree = generate_constaint_tree( CtPlBiCn, CtPlBiCn_not )
	write.tree( 
		porifera_sister_constraint_tree, 
		file = paste( constraint_tree_path, seq_matrix@matrix_name, ".porifera_sister_constraint.tree", sep="") 
	)
	
	PoPlBiCn_not = rownames(seq_matrix)[ ! rownames(seq_matrix) %in% PoPlBiCn ]
	ctenophora_sister_constraint_tree = generate_constaint_tree( PoPlBiCn, PoPlBiCn_not )
	write.tree( 
		ctenophora_sister_constraint_tree, 
		file = paste( constraint_tree_path, seq_matrix@matrix_name, ".ctenophora_sister_constraint.tree", sep="") 
	)	
}

lapply(sequence_matrices, generate_constraint_trees)

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

trees_path_iqtree = "../trees_new/iqtree"
iqtree_ext = "\\.treefile$"
file_names_iqtree = list.files( path = trees_path_iqtree, pattern = iqtree_ext, full.names = TRUE )

trees_path_pb = "../trees_new/phylobayes"
pb_tree_ext = "\\.con\\.tre$"
file_names_pb = list.files( path = trees_path_pb, pattern = pb_tree_ext, full.names = TRUE )

# Temp fix re issue #8
file_names_iqtree = file_names_iqtree[ ! grepl("Philippe2009", file_names_iqtree) ]
file_names_iqtree = file_names_iqtree[ ! grepl("Nosenko2013", file_names_iqtree) ]


trees = foreach( tree_file = file_names_iqtree ) %dopar%
	parse_tree_iqtree( tree_file, taxonomy_reference )

analyses_new = lapply(
	trees,
	function( tree ){
		data.frame(
			matrix = tree$matrix,
			model = tree$model,
			modelfinder = tree$modelfinder,
			clade = tree$sampling,
			support_ctenophora_sister = tree$ctenophora_sister * 100,
			support_porifera_sister = tree$porifera_sister * 100,
			stringsAsFactors = FALSE
		)
	}
) %>% 
	bind_rows()

analyses_new$inference = rep( "ML", nrow( analyses_new ) )

# Summarize result
analyses_new$result = "Unresolved"
analyses_new$result[ ( analyses_new$inference == "Bayesian" ) & (analyses_new$support_porifera_sister >= posterior_prob_threshold) ] = "Porifera-sister"
analyses_new$result[ ( analyses_new$inference == "Bayesian" ) & (analyses_new$support_ctenophora_sister >= posterior_prob_threshold) ] = "Ctenophora-sister"
analyses_new$result[ ( analyses_new$inference == "ML" ) &       (analyses_new$support_porifera_sister >= bootstrap_threshold) ] = "Porifera-sister"
analyses_new$result[ ( analyses_new$inference == "ML" ) &       (analyses_new$support_ctenophora_sister >= bootstrap_threshold) ] = "Ctenophora-sister"
analyses_new$result = factor( analyses_new$result )

# Parse model components
analyses_new$model_summary = analyses_new$model
analyses_new$model_summary = factor( analyses_new$model_summary, levels=c("WAG+C60+F+G", "LG+C60+F+G", "poisson_C60", "GTR20", "WAG", "LG") )

# Partition comparison across matrices

n_total_partitions =
	partitions_all %>% 
	group_by(matrix) %>%
	summarize("n_total_partitions"=n_distinct(partition)) 
n_compandBUSCO = 
	partition_map_global %>%
		group_by(matrix) %>%
		summarize(
			"n_components"=n_distinct(gene),
			"n_distinct_BUSCO"=n_distinct(BUSCO_ID)
		)
n_components_with_BUSCO =
	partition_map_global %>%
		filter(BUSCO_ID != "") %>%
		group_by(matrix,gene) %>%
		summarize(n()) %>%
		group_by(matrix) %>%
		tally(name="n_components_with_BUSCO")
n_ribo = 
	partition_map_global %>% 
		group_by(matrix,gene) %>% 
		tally(ribo_found) %>% 
		group_by(matrix) %>% 
		tally()
discarded_parts = 
	read_tsv("../reconciliation/blast/graphs/discarded_nodes.tsv") %>%
	group_by(matrix) %>%
	summarize("n_partitions_discarded"=n_distinct(partition_name))

partition_network_summary = 
	n_total_partitions%>%
	left_join(n_compandBUSCO,          by="matrix") %>%
	left_join(n_components_with_BUSCO, by="matrix") %>%
	left_join(n_ribo,                  by="matrix") %>%
	left_join(discarded_parts,         by="matrix") %>% 
	mutate(n_partitions_discarded = replace_na(n_partitions_discarded, 0))

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