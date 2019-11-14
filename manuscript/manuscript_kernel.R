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

setwd("/animal_root/manuscript")

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
	read_tsv( "../data_processed/tables/previously_published_analyses.tsv", col_types = 'icccccccccccclccccccc' ) %>%
	mutate( clade = factor( clade, levels= c( 'Choanimalia', 'Holozoa', 'Opisthokonta' ) ) )

taxonomy_reference = read_tsv("../reconciliation/taxonomy_info/taxon_table.tsv")

# raw columns:
# component_number, matrix, partition_name, edges, nodes_in_component, component_density, BUSCO_ID, BUSCO_description,
# SwissProt_accession, SwissProt_description, GO_annotations, ribo_found
partition_map_global =
	read_tsv("../reconciliation/blast/graphs/partition_components_split_annotated.tsv") %>%
	dplyr::rename(partition=partition_name) %>%
	mutate( component_number = as.character(component_number) )

analyses_published$result = "Unresolved"
analyses_published$result[ ( analyses_published$inference == "Bayesian" ) & (analyses_published$support_porifera_sister >= posterior_prob_threshold) ] = "Porifera-sister"
analyses_published$result[ ( analyses_published$inference == "Bayesian" ) & (analyses_published$support_ctenophora_sister >= posterior_prob_threshold) ] = "Ctenophora-sister"
analyses_published$result[ ( analyses_published$inference == "ML" ) &       (analyses_published$support_porifera_sister >= bootstrap_threshold) ] = "Porifera-sister"
analyses_published$result[ ( analyses_published$inference == "ML" ) &       (analyses_published$support_ctenophora_sister >= bootstrap_threshold) ] = "Ctenophora-sister"
analyses_published$result = factor( analyses_published$result )

analyses_published$model_combined = factor( analyses_published$model_combined, levels=c("WAG", "LG", "GTR", "data partitioning", "Recoding + GTR", "Recoding + GTR + CAT", "Poisson + CAT", "GTR + CAT" ) )

# Matrix taxon composition

clades = c( "Fungi", "Ichthyosporea", "Filasterea", "Choanoflagellida", "Ctenophora", "Porifera", "Placozoa", "Bilateria", "Cnidaria" )

taxa =
	taxonomy_reference %>%
	distinct( relabelled_name, clade_assignment, ncbi_tax_id ) %>%
	dplyr::rename( taxon=relabelled_name, clade=clade_assignment ) %>%
	mutate( clade = factor( clade, levels=clades )  )

matrix_path = "../data_processed/matrices"
phylip_file_names = list.files(path = matrix_path, pattern = ".+\\.phy$", full.names=TRUE)

sequence_matrices = foreach( phylip_file=phylip_file_names) %dopar%
	parse_phylip( phylip_file )

# Make contraint trees for each matrix
constraint_tree_path = "../trees_new/constraint_trees/"



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
	group_by( component_number ) %>%
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
mask = lower.tri(matrix(nrow=sqrt(n), ncol=sqrt(n)))
dim(mask) = NULL
matrix_overlap = matrix_overlap[mask,]


# New analyses of published matrices

trees_path_iqtree = "../trees_new/iqtree"
iqtree_ext = "\\.treefile$"
file_names_iqtree = list.files( path = trees_path_iqtree, pattern = iqtree_ext, full.names = TRUE )

# Temp fix re issue #8
file_names_iqtree = file_names_iqtree[ ! grepl("Philippe2009", file_names_iqtree) ]
file_names_iqtree = file_names_iqtree[ ! grepl("Nosenko2013", file_names_iqtree) ]

# read iqtrees
trees_iq = foreach( tree_file = file_names_iqtree ) %dopar%
	parse_tree_iqtree( tree_file, taxonomy_reference )

trees_paths_pb = list.dirs("../trees_new/")
trees_paths_pb <- trees_paths_pb[ grepl("_phylobayes", trees_paths_pb) ]
pb_tree_ext = "\\.con\\.tre$"
file_names_pb = list.files( path = trees_paths_pb, pattern = pb_tree_ext, full.names = TRUE )

# read pb trees
trees_pb = foreach( tree_file = file_names_pb ) %dopar%
           parse_tree_pb( tree_file, taxonomy_reference )

analyses_new = lapply(
	c( trees_iq, trees_pb ),
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

analyses_new$inference = c( rep( "ML", length( trees_iq ) ), rep( "Bayesian", length( trees_pb ) ))

# Summarize result
analyses_new$result = "Unresolved"
analyses_new$result[ ( analyses_new$inference == "Bayesian" ) & (analyses_new$support_porifera_sister >= posterior_prob_threshold) ] = "Porifera-sister"
analyses_new$result[ ( analyses_new$inference == "Bayesian" ) & (analyses_new$support_ctenophora_sister >= posterior_prob_threshold) ] = "Ctenophora-sister"
analyses_new$result[ ( analyses_new$inference == "ML" ) &       (analyses_new$support_porifera_sister >= bootstrap_threshold) ] = "Porifera-sister"
analyses_new$result[ ( analyses_new$inference == "ML" ) &       (analyses_new$support_ctenophora_sister >= bootstrap_threshold) ] = "Ctenophora-sister"
analyses_new$result = factor( analyses_new$result )

# Parse model components
analyses_new$model_summary = analyses_new$model
analyses_new$model_summary = factor( analyses_new$model_summary, levels=c("WAG", "GTR20", "poisson_C60", "WAG+C60+F+G", "LG+C60+F+G", "CAT+F81"))

# Partition comparison across matrices

n_total_partitions =
	partitions_all %>%
	group_by(matrix) %>%
	summarize("n_total_partitions"=n_distinct(partition))
n_compandBUSCO =
	partition_map_global %>%
		group_by(matrix) %>%
		summarize(
			"n_components"=n_distinct(component_number),
			"n_distinct_BUSCO"=n_distinct(BUSCO_ID)
		)
n_components_with_BUSCO =
	partition_map_global %>%
		filter(BUSCO_ID != "") %>%
		group_by( matrix, component_number ) %>%
		summarize(n()) %>%
		group_by(matrix) %>%
		tally(name="n_components_with_BUSCO")
n_ribo =
	partition_map_global %>%
		group_by( matrix, component_number ) %>%
		tally(ribo_found) %>%
		group_by(matrix) %>%
		tally()
discarded_parts =
	read_tsv("../reconciliation/blast/graphs/discarded_nodes.tsv") %>%
	group_by(matrix) %>%
	summarize("n_partitions_discarded"=n_distinct(partition_name))

partition_network_summary =
	n_total_partitions%>%
	# left_join(n_compandBUSCO,          by="matrix") %>%
	left_join(n_components_with_BUSCO, by="matrix") %>%
	left_join(n_ribo,                  by="matrix") %>%
	left_join(discarded_parts,         by="matrix") %>%
	mutate(n_partitions_discarded = replace_na(n_partitions_discarded, 0))
### Summarize phylogenetic signal by genes
Philippe2009  <- read_tsv("../trees_new/AU_test/Philippe2009_C60_genewise_lnL.txt")
Philippe2009$ABS <- abs(Philippe2009$`tr1_log-likelihood` - Philippe2009$`tr2_log-likelihood`)
Philippe2009$result[ (Philippe2009$tree_supported == "tr1" ) & (Philippe2009$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Philippe2009$result[ (Philippe2009$tree_supported == "tr1" ) & (Philippe2009$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Philippe2009$result[ (Philippe2009$tree_supported == "tr2" ) & (Philippe2009$ABS >= 1) ] = "Porifera-sister-strong-signal"
Philippe2009$result[ (Philippe2009$tree_supported == "tr2" ) & (Philippe2009$ABS <1) ] = "Porifera-sister-weak-signal"

Philippe2009 %>%
	count(result) %>%
	mutate(perc = n/ nrow(Philippe2009)) -> Philippe2009_count
Philippe2009_count$study <- 'Philippe2009_Opisthokonta_C60'


Philippe2009_only_choanozoa  <- read_tsv("../trees_new/AU_test/Philippe2009_only_choanozoa_C60_genewise_lnL.txt")
Philippe2009_only_choanozoa$ABS <- abs(Philippe2009_only_choanozoa$`tr1_log-likelihood` - Philippe2009_only_choanozoa$`tr2_log-likelihood`)
Philippe2009_only_choanozoa$result[ (Philippe2009_only_choanozoa$tree_supported == "tr1" ) & (Philippe2009_only_choanozoa$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Philippe2009_only_choanozoa$result[ (Philippe2009_only_choanozoa$tree_supported == "tr1" ) & (Philippe2009_only_choanozoa$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Philippe2009_only_choanozoa$result[ (Philippe2009_only_choanozoa$tree_supported == "tr2" ) & (Philippe2009_only_choanozoa$ABS >= 1) ] = "Porifera-sister-strong-signal"
Philippe2009_only_choanozoa$result[ (Philippe2009_only_choanozoa$tree_supported == "tr2" ) & (Philippe2009_only_choanozoa$ABS <1) ] = "Porifera-sister-weak-signal"

Philippe2009_only_choanozoa %>%
	count(result) %>%
	mutate(perc = n/ nrow(Philippe2009_only_choanozoa)) -> Philippe2009_only_choanozoa_count
Philippe2009_only_choanozoa_count$study <- 'Philippe2009_choanozoa_C60'


Ryan2013  <- read_tsv("../trees_new/AU_test/Ryan2013_est_C60_genewise_lnL.txt")
Ryan2013$ABS <- abs(Ryan2013$`tr1_log-likelihood` - Ryan2013$`tr2_log-likelihood`)
Ryan2013$result[ (Ryan2013$tree_supported == "tr1" ) & (Ryan2013$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Ryan2013$result[ (Ryan2013$tree_supported == "tr1" ) & (Ryan2013$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Ryan2013$result[ (Ryan2013$tree_supported == "tr2" ) & (Ryan2013$ABS >= 1) ] = "Porifera-sister-strong-signal"
Ryan2013$result[ (Ryan2013$tree_supported == "tr2" ) & (Ryan2013$ABS <1) ] = "Porifera-sister-weak-signal"

Ryan2013 %>%
	count(result) %>%
	mutate(perc = n/ nrow(Ryan2013)) -> Ryan2013_count
Ryan2013_count$study <- 'Ryan2013_Opisthokonta_C60'

Ryan2013_only_choanozoa  <- read_tsv("../trees_new/AU_test/Ryan2013_est_only_choanozoa_C60_genewise_lnL.txt")
Ryan2013_only_choanozoa$ABS <- abs(Ryan2013_only_choanozoa$`tr1_log-likelihood` - Ryan2013_only_choanozoa$`tr2_log-likelihood`)
Ryan2013_only_choanozoa$result[ (Ryan2013_only_choanozoa$tree_supported == "tr1" ) & (Ryan2013_only_choanozoa$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Ryan2013_only_choanozoa$result[ (Ryan2013_only_choanozoa$tree_supported == "tr1" ) & (Ryan2013_only_choanozoa$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Ryan2013_only_choanozoa$result[ (Ryan2013_only_choanozoa$tree_supported == "tr2" ) & (Ryan2013_only_choanozoa$ABS >= 1) ] = "Porifera-sister-strong-signal"
Ryan2013_only_choanozoa$result[ (Ryan2013_only_choanozoa$tree_supported == "tr2" ) & (Ryan2013_only_choanozoa$ABS <1) ] = "Porifera-sister-weak-signal"

Ryan2013_only_choanozoa %>%
	count(result) %>%
	mutate(perc = n/ nrow(Ryan2013_only_choanozoa)) -> Ryan2013_only_choanozoa_count
Ryan2013_only_choanozoa_count$study <- 'Ryan2013_choanozoa_C60'

Whelan2017_full  <- read_tsv("../trees_new/AU_test/Whelan2017_full_C60_genewise_lnL.txt")
Whelan2017_full$ABS <- abs(Whelan2017_full$`tr1_log-likelihood` - Whelan2017_full$`tr2_log-likelihood`)
Whelan2017_full$result[ (Whelan2017_full$tree_supported == "tr1" ) & (Whelan2017_full$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Whelan2017_full$result[ (Whelan2017_full$tree_supported == "tr1" ) & (Whelan2017_full$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Whelan2017_full$result[ (Whelan2017_full$tree_supported == "tr2" ) & (Whelan2017_full$ABS >= 1) ] = "Porifera-sister-strong-signal"
Whelan2017_full$result[ (Whelan2017_full$tree_supported == "tr2" ) & (Whelan2017_full$ABS <1) ] = "Porifera-sister-weak-signal"

Whelan2017_full %>%
	count(result) %>%
	mutate(perc = n/ nrow(Whelan2017_full)) -> Whelan2017_full_count
Whelan2017_full_count$study <- 'Whelan2017_holozoa_C60'

Whelan2017_full_only_choanozoa  <- read_tsv("../trees_new/AU_test/Whelan2017_full_only_choanozoa_C60_genewise_lnL.txt")
Whelan2017_full_only_choanozoa$ABS <- abs(Whelan2017_full_only_choanozoa$`tr1_log-likelihood` - Whelan2017_full_only_choanozoa$`tr2_log-likelihood`)
Whelan2017_full_only_choanozoa$result[ (Whelan2017_full_only_choanozoa$tree_supported == "tr1" ) & (Whelan2017_full_only_choanozoa$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Whelan2017_full_only_choanozoa$result[ (Whelan2017_full_only_choanozoa$tree_supported == "tr1" ) & (Whelan2017_full_only_choanozoa$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Whelan2017_full_only_choanozoa$result[ (Whelan2017_full_only_choanozoa$tree_supported == "tr2" ) & (Whelan2017_full_only_choanozoa$ABS >= 1) ] = "Porifera-sister-strong-signal"
Whelan2017_full_only_choanozoa$result[ (Whelan2017_full_only_choanozoa$tree_supported == "tr2" ) & (Whelan2017_full_only_choanozoa$ABS <1) ] = "Porifera-sister-weak-signal"

Whelan2017_full_only_choanozoa %>%
	count(result) %>%
	mutate(perc = n/ nrow(Whelan2017_full_only_choanozoa)) -> Whelan2017_full_only_choanozoa_count
Whelan2017_full_only_choanozoa_count$study <- 'Whelan2017_choanozoa_C60'



Philippe2009_WAG  <- read_tsv("../trees_new/AU_test/Philippe2009_WAG_genewise_lnL.txt")
Philippe2009_WAG$ABS <- abs(Philippe2009_WAG$`tr1_log-likelihood` - Philippe2009_WAG$`tr2_log-likelihood`)
Philippe2009_WAG$result[ (Philippe2009_WAG$tree_supported == "tr1" ) & (Philippe2009_WAG$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Philippe2009_WAG$result[ (Philippe2009_WAG$tree_supported == "tr1" ) & (Philippe2009_WAG$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Philippe2009_WAG$result[ (Philippe2009_WAG$tree_supported == "tr2" ) & (Philippe2009_WAG$ABS >= 1) ] = "Porifera-sister-strong-signal"
Philippe2009_WAG$result[ (Philippe2009_WAG$tree_supported == "tr2" ) & (Philippe2009_WAG$ABS <1) ] = "Porifera-sister-weak-signal"

Philippe2009_WAG %>%
	count(result) %>%
	mutate(perc = n/ nrow(Philippe2009_WAG)) -> Philippe2009_WAG_count
Philippe2009_WAG_count$study <- 'Philippe2009_Opisthokonta_WAG'


Philippe2009_WAG_only_choanozoa  <- read_tsv("../trees_new/AU_test/Philippe2009_only_choanozoa_WAG_genewise_lnL.txt")
Philippe2009_WAG_only_choanozoa$ABS <- abs(Philippe2009_WAG_only_choanozoa$`tr1_log-likelihood` - Philippe2009_WAG_only_choanozoa$`tr2_log-likelihood`)
Philippe2009_WAG_only_choanozoa$result[ (Philippe2009_WAG_only_choanozoa$tree_supported == "tr1" ) & (Philippe2009_WAG_only_choanozoa$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Philippe2009_WAG_only_choanozoa$result[ (Philippe2009_WAG_only_choanozoa$tree_supported == "tr1" ) & (Philippe2009_WAG_only_choanozoa$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Philippe2009_WAG_only_choanozoa$result[ (Philippe2009_WAG_only_choanozoa$tree_supported == "tr2" ) & (Philippe2009_WAG_only_choanozoa$ABS >= 1) ] = "Porifera-sister-strong-signal"
Philippe2009_WAG_only_choanozoa$result[ (Philippe2009_WAG_only_choanozoa$tree_supported == "tr2" ) & (Philippe2009_WAG_only_choanozoa$ABS <1) ] = "Porifera-sister-weak-signal"

Philippe2009_WAG_only_choanozoa %>%
	count(result) %>%
	mutate(perc = n/ nrow(Philippe2009_WAG_only_choanozoa)) -> Philippe2009_WAG_only_choanozoa_count
Philippe2009_WAG_only_choanozoa_count$study <- 'Philippe2009_choanozoa_WAG'


Ryan2013_WAG  <- read_tsv("../trees_new/AU_test/Ryan2013_est_WAG_genewise_lnL.txt")
Ryan2013_WAG$ABS <- abs(Ryan2013_WAG$`tr1_log-likelihood` - Ryan2013_WAG$`tr2_log-likelihood`)
Ryan2013_WAG$result[ (Ryan2013_WAG$tree_supported == "tr1" ) & (Ryan2013_WAG$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Ryan2013_WAG$result[ (Ryan2013_WAG$tree_supported == "tr1" ) & (Ryan2013_WAG$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Ryan2013_WAG$result[ (Ryan2013_WAG$tree_supported == "tr2" ) & (Ryan2013_WAG$ABS >= 1) ] = "Porifera-sister-strong-signal"
Ryan2013_WAG$result[ (Ryan2013_WAG$tree_supported == "tr2" ) & (Ryan2013_WAG$ABS <1) ] = "Porifera-sister-weak-signal"

Ryan2013_WAG %>%
	count(result) %>%
	mutate(perc = n/ nrow(Ryan2013_WAG)) -> Ryan2013_WAG_count
Ryan2013_WAG_count$study <- 'Ryan2013_Opisthokonta_WAG'

Ryan2013_WAG_only_choanozoa  <- read_tsv("../trees_new/AU_test/Ryan2013_est_only_choanozoa_WAG_genewise_lnL.txt")
Ryan2013_WAG_only_choanozoa$ABS <- abs(Ryan2013_WAG_only_choanozoa$`tr1_log-likelihood` - Ryan2013_WAG_only_choanozoa$`tr2_log-likelihood`)
Ryan2013_WAG_only_choanozoa$result[ (Ryan2013_WAG_only_choanozoa$tree_supported == "tr1" ) & (Ryan2013_WAG_only_choanozoa$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Ryan2013_WAG_only_choanozoa$result[ (Ryan2013_WAG_only_choanozoa$tree_supported == "tr1" ) & (Ryan2013_WAG_only_choanozoa$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Ryan2013_WAG_only_choanozoa$result[ (Ryan2013_WAG_only_choanozoa$tree_supported == "tr2" ) & (Ryan2013_WAG_only_choanozoa$ABS >= 1) ] = "Porifera-sister-strong-signal"
Ryan2013_WAG_only_choanozoa$result[ (Ryan2013_WAG_only_choanozoa$tree_supported == "tr2" ) & (Ryan2013_WAG_only_choanozoa$ABS <1) ] = "Porifera-sister-weak-signal"

Ryan2013_WAG_only_choanozoa %>%
	count(result) %>%
	mutate(perc = n/ nrow(Ryan2013_WAG_only_choanozoa)) -> Ryan2013_WAG_only_choanozoa_count
Ryan2013_WAG_only_choanozoa_count$study <- 'Ryan2013_choanozoa_WAG'

Whelan2017_full_WAG  <- read_tsv("../trees_new/AU_test/Whelan2017_full_WAG_genewise_lnL.txt")
Whelan2017_full_WAG$ABS <- abs(Whelan2017_full_WAG$`tr1_log-likelihood` - Whelan2017_full_WAG$`tr2_log-likelihood`)
Whelan2017_full_WAG$result[ (Whelan2017_full_WAG$tree_supported == "tr1" ) & (Whelan2017_full_WAG$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Whelan2017_full_WAG$result[ (Whelan2017_full_WAG$tree_supported == "tr1" ) & (Whelan2017_full_WAG$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Whelan2017_full_WAG$result[ (Whelan2017_full_WAG$tree_supported == "tr2" ) & (Whelan2017_full_WAG$ABS >= 1) ] = "Porifera-sister-strong-signal"
Whelan2017_full_WAG$result[ (Whelan2017_full_WAG$tree_supported == "tr2" ) & (Whelan2017_full_WAG$ABS <1) ] = "Porifera-sister-weak-signal"

Whelan2017_full_WAG %>%
	count(result) %>%
	mutate(perc = n/ nrow(Whelan2017_full_WAG)) -> Whelan2017_full_WAG_count
Whelan2017_full_WAG_count$study <- 'Whelan2017_holozoa_WAG'

Whelan2017_full_WAG_only_choanozoa  <- read_tsv("../trees_new/AU_test/Whelan2017_full_only_choanozoa_WAG_genewise_lnL.txt")
Whelan2017_full_WAG_only_choanozoa$ABS <- abs(Whelan2017_full_WAG_only_choanozoa$`tr1_log-likelihood` - Whelan2017_full_WAG_only_choanozoa$`tr2_log-likelihood`)
Whelan2017_full_WAG_only_choanozoa$result[ (Whelan2017_full_WAG_only_choanozoa$tree_supported == "tr1" ) & (Whelan2017_full_WAG_only_choanozoa$ABS >= 1) ] = "Ctenophore-sister-strong-signal"
Whelan2017_full_WAG_only_choanozoa$result[ (Whelan2017_full_WAG_only_choanozoa$tree_supported == "tr1" ) & (Whelan2017_full_WAG_only_choanozoa$ABS < 1) ] = "Ctenophore-sister-weak-signal"
Whelan2017_full_WAG_only_choanozoa$result[ (Whelan2017_full_WAG_only_choanozoa$tree_supported == "tr2" ) & (Whelan2017_full_WAG_only_choanozoa$ABS >= 1) ] = "Porifera-sister-strong-signal"
Whelan2017_full_WAG_only_choanozoa$result[ (Whelan2017_full_WAG_only_choanozoa$tree_supported == "tr2" ) & (Whelan2017_full_WAG_only_choanozoa$ABS <1) ] = "Porifera-sister-weak-signal"

Whelan2017_full_WAG_only_choanozoa %>%
	count(result) %>%
	mutate(perc = n/ nrow(Whelan2017_full_WAG_only_choanozoa)) -> Whelan2017_full_WAG_only_choanozoa_count
Whelan2017_full_WAG_only_choanozoa_count$study <- 'Whelan2017_choanozoa_WAG'

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
