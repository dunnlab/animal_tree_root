library( Biostrings )
library( "R.utils")
library( ape )
library( hutan )
library( tidyverse )


#######################################
## Classes


#' An S4 class to store multiple sequence alignments with partition data
#' 
#' @slot partitions  The partitions
setClass(
	Class = "PartitionedMultipleAlignment",
	representation = representation(
		partitions = "data.frame",
		manuscript_name = "character",
		matrix_name = "character",
		taxon_map = "vector"
	),
	contains = "MultipleAlignment"
)

#' Construct a PartitionedMultipleAlignment from an alignment file and a partition file
#'
#' @param alignment_file Path to the gene alignment
#' @param partition_file Path to the partition file in nexus format
#' @param taxon_map_global A dataframe where the first column corresponds to sequence names in the alignment and the second column to clade names for the sequences
#' @param partition_map_global A dataframe with columns manuscript, matrix, partition, global_partition_id that maps the locally defined partition to a global partition
#' @return A PartitionedMultipleAlignment object
#' @export
PartitionedMultipleAlignment = function( alignment_file, partition_file=NULL, taxon_map_global=NULL, partition_map_global=NULL, manuscript_name="", matrix_name="", alignment_format="phylip" ) {
	object = readAAMultipleAlignment( filepath = alignment_file, format=alignment_format )
	class(object) = "PartitionedMultipleAlignment"
	
	object@manuscript_name = manuscript_name
	object@matrix_name = matrix_name
		
	conn = file( partition_file, open="r")
	lines = readLines( conn )
	lines = str_replace( lines, "charset", "CHARSET" )
	lines = lines[ grepl("CHARSET", lines) ]
	close( conn )
	
	partitions = 
		lapply( 
			lines, 
			function( line ){
				fields = str_match(line, "CHARSET (.+?)\\s*=\\s*(\\d+)\\s*-\\s*(\\d+)")
				D = data.frame( 
					partition=fields[2], 
					start=as.integer(fields[3]), 
					stop=as.integer(fields[4]),
					stringsAsFactors = FALSE
					)
				return( D )
			}
		) %>%
		bind_rows()
	
	if( ! is.null(partition_map_global) ){
		partition_map_global %<>%
			filter( matrix == matrix_name )
		
		partitions %<>% 
			left_join( partition_map_global ) %>%
			select( -matrix )
		
		# Add names to the partitions that are not assigned to components. Take these from the original 
		# partition names, but prepend matrix name to be sure they are unique across matrices
		
		partitions$component_number[ is.na(partitions$component_number) ] = paste( matrix_name, partitions$partition[ is.na(partitions$component_number) ], "NA", sep="_")
		
	}
	
	object@partitions = partitions
	
	if( ! is.null(taxon_map_global) ){
		sequence_name = rownames( object )
		object@taxon_map = taxon_map_global[ match( sequence_name, taxon_map_global[[1]] ), 2 ][[1]]
	}
	
	object
}

#####################################################
## GGPlot extension of multiple sequence alignments

# Builds on https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html


#' Calculates rectangles that represent clades in an alignment
#'
#' @param data A PartitionedMultipleAlignment object
#' @return A dataframe
#' @export
clade_rects = function( data ){
	n_genes = nrow(data@partitions)
	clade_table = table( data@taxon_map ) %>% as.vector()
	clade_sums = as.vector(clade_table)
	clade_cumsums = cumsum(clade_sums)
	
	clades = factor( levels( data@taxon_map ), levels( data@taxon_map ))
	
	
	# Create a dataframe where each row is a clade and the columns describe how to draw the rectangle for the clade
	D = data.frame( clade=clades )
	
	D$xmin = rep( 0, length(clades) )
	D$xmax = rep( n_genes, length(clades) )
	D$ymin = c(0, clade_cumsums[1:length(clade_cumsums)-1])
	D$ymax = clade_cumsums
	D$manuscript = data@manuscript_name
	D$matrix = data@matrix_name
	D$manuscript_matrix = paste(D$manuscript, D$matrix, sep="_")
	
	D %<>% mutate( manuscript_matrix = as.factor(manuscript_matrix) )
	
	D
}


StatClades = 
	ggproto(
		"StatClades", 
		Stat,
		compute_group = 
			function(data, scales) {
				clade_rects( data )
			},
		required_aes = c( "x", "y" )
)

stat_clades = function(mapping = NULL, data = NULL, geom = "rect",
											 position = "identity", na.rm = FALSE, show.legend = NA, 
											 inherit.aes = TRUE, ...) {
	layer(
		stat = StatClades, data = data, mapping = mapping, geom = geom, 
		position = position, show.legend = show.legend, inherit.aes = inherit.aes,
		params = list(na.rm = na.rm, ...)
	)
}



compute_matrix_overlap = function( msa1, msa2 ){
	
	n_gene_overlap = sum( na.omit(msa1@partitions$gene) %in% na.omit(msa2@partitions$gene) )

	n_species_overlap = sum( rownames(msa1) %in% rownames(msa2) )
	
	tibble(
		manuscript_name_1 = msa1@manuscript_name,
		matrix_name_1 = msa1@matrix_name,
		n_species_1 = nrow( msa1 ),
		n_partitions_1 = nrow( msa1@partitions ),
		manuscript_name_2 = msa2@manuscript_name,
		matrix_name_2 = msa2@matrix_name,
		n_species_2 = nrow( msa2 ),
		n_partitions_2 = nrow( msa2@partitions ),
		n_species_overlap = n_species_overlap,
		n_gene_overlap = n_gene_overlap
		
	)
}

overlap_rects = function( msa1, msa2 ){
	overlap = compute_matrix_overlap( msa1, msa2 )
	
	# Create a dataframe where each row is one of the matrices in the pairwise comparison and 
	# the columns describe how to draw the rectangle for matrix
	
	# First row is msa1 rectangle, second row is msa2 rectangle
	D = data.frame( MSA=c("Matrix_1", "Matrix_2") )
	
	D$xmin = c( 0, overlap$n_partitions_1 - overlap$n_gene_overlap )
	D$xmax = c( overlap$n_partitions_1, overlap$n_partitions_1 + overlap$n_partitions_2 - overlap$n_gene_overlap )
	D$ymin = c(0, overlap$n_species_1 - overlap$n_species_overlap )
	D$ymax = c( overlap$n_species_1, overlap$n_species_1 + overlap$n_species_2 - overlap$n_species_overlap )
	D$manuscript_1 = msa1@manuscript_name
	D$manuscript_2 = msa2@manuscript_name
	D$matrix_1 = paste(msa1@manuscript_name, msa1@matrix_name, sep="_")
	D$matrix_2 = paste(msa2@manuscript_name, msa2@matrix_name, sep="_")
	
	D
	
	
}

#####################################################
## Tree Parsing
#' Parse iqtree or pb trees and return an overloaded ape::phylo with additional items:
#'  (additional named items NA if not otherwise defined)
#' 	
#' matrix: matrix name
#' model: model name
#' modelfinder: TRUE if best model found
#' clades: clades of the tips of the tree
#' sampling = list of outgroups sampled
#' ctenophora_sister: support for ctenophora as sister
#' porifera_sister: support for porifera as sister
#' bpdiff: bpcomp's maxdiff
#' CtPlBiCn: set of taxa in the group Ctenophora+Placozoa+Bilateria+Cnidaria, the clade that exists under Proifera-sister
#' PoPlBiCn: set of taxa in the group Porifera+Placozoa+Bilateria+Cnidaria, the clade that exists under Ctenophora-sister
#'
#' @param tree_file a string with path to a tree file
#' @param taxonomy_reference data.frame or tibble generated from taxon_table.tsv
#' @return ape::phylo with additional items
#' @export

parse_tree = 	function( tree_file,  taxonomy_reference ){
	
	tree = read.tree( file = tree_file )

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
	
	if( any( is.na( outgroup ) ) ){
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
		
	# initialize additional named items
	tree$matrix = NA
	tree$model = NA
	tree$modelfinder = NA
	tree$clades = clades
	tree$sampling = sampling
	tree$ctenophora_sister = NA
	tree$porifera_sister = NA
	tree$bpdiff = NA
	
	if( any( is.na( clades ) ) ){
		missing_clades = tree$tip.label[ is.na( clades ) ]
		stop( str_c( c( "Some tips have missing clade names:", missing_clades ), collapse=" " ) )
	}
	
	if( ! any( clades == "Ctenophora" ) ){
		stop( "There are no ctenophores in the tree, can't test rooting topology" )
	} 
	
	if( ! any( clades == "Porifera" ) ){
		stop( "There are no sponges in the tree, can't test rooting topology" )
	} 
	
	# Get the set of taxa in the group Ctenophora+Placozoa+Bilateria+Cnidaria, the clade that exists under Proifera-sister
	tree$CtPlBiCn = tree$tip.label[ tree$clades %in% c("Ctenophora", "Placozoa", "Bilateria", "Cnidaria") ]
	
	# Get the set of taxa in the group Porifera+Placozoa+Bilateria+Cnidaria, the clade that exists under Ctenophora-sister
	tree$PoPlBiCn = tree$tip.label[ tree$clades %in% c("Porifera", "Placozoa", "Bilateria", "Cnidaria") ]
	
	return( tree )
}

read_pb_treelist = function( trees_file, burnin, subsample_every ){
	tree_strings = read_lines( trees_file, skip = burnin)
	tree_strings = tree_strings[c(TRUE, rep(FALSE,subsample_every-1))]
	return(tree_strings)
}

parse_tree_pb = function( tree_file, taxonomy_reference ){
	# set up filenames
	tree_file_path = normalizePath( tree_file )
	print( tree_file_path )
	tree_ext = "\\.con\\.tre$"
	bpdiff_file_path = sub( tree_ext, ".bpdiff", tree_file_path )
	pb_treelist = sprintf(
		sub( str_c("bpcomp",tree_ext), "Chain%d.treelist", tree_file_path ),
		1:2
	)
	
	if(!file.exists(bpdiff_file_path)){
		warning( str_c( "Missing bpdiff file: ", bpdiff_file_path ) )
		return( NA )
	}
	
	mdiff_line = strsplit( 
		grep( "^maxdiff", readLines( bpdiff_file_path ), value = TRUE), 
		":") %>% unlist()
	mdiff = as.numeric(mdiff_line[2])
	
	tree = parse_tree( tree_file_path, taxonomy_reference )
	tree$bpdiff = mdiff
	
	# 
	burnins = round( as.vector(sapply( pb_treelist, countLines ) * 0.2 ))
	keep_every = 100
	treelist_text = c( unlist( mapply( read_pb_treelist, pb_treelist, burnins, keep_every ) ) )
	sample_trees = read.tree( text = treelist_text )
	
	Ctenophora_sister = lapply(
		sample_trees,
		function( phy ){
			is_monophyletic( phy, tree$PoPlBiCn )
		}
	) %>% 
		unlist() %>%
		mean()
	
	Porifera_sister = lapply(
		sample_trees,
		function( phy ){
			is_monophyletic( phy, tree$CtPlBiCn )
		}
	) %>% 
		unlist() %>%
		mean()
	
	tree$ctenophora_sister = Ctenophora_sister
	tree$porifera_sister = Porifera_sister
	
	return(tree)
}

parse_tree_iqtree = function( tree_file, taxonomy_reference ){
	tree_file_path = normalizePath( tree_file )
	print( tree_file_path )
	tree_wd = dirname( tree_file_path )
	filename_parts = strsplit( basename( tree_file ), "\\." ) %>% unlist()
	matrix_name = filename_parts[1]
	model_name = filename_parts[2]
	log_file_path = file.path( tree_wd, str_c( matrix_name, model_name, "log", sep = "." ) )
	bootstrap_file_path = file.path( tree_wd, str_c( matrix_name, model_name, "ufboot", sep = "." ) )
	
	if(!file.exists(log_file_path)){
		warning( str_c( "Missing log file: ", log_file_path ) )
		return( NA )
	}
	
	if(!file.exists(bootstrap_file_path)){
		warning( str_c( "Missing bootstrap file: ", bootstrap_file_path ) )
		return( NA )
	}
	
	tree = parse_tree( tree_file_path, taxonomy_reference )
	
	# Read the bootstrap trees
	bootstrap_trees = read.tree(file = bootstrap_file_path)
	
	Ctenophora_sister = lapply(
		bootstrap_trees,
		function( phy ){
			is_monophyletic( phy, tree$PoPlBiCn )
		}
	) %>% 
		unlist() %>%
		mean()
	
	Porifera_sister = lapply(
		bootstrap_trees,
		function( phy ){
			is_monophyletic( phy, tree$CtPlBiCn )
		}
	) %>% 
		unlist() %>%
		mean()				
	
	# Stuff a few other things into the tree object
	tree$matrix = matrix_name
	tree$model = model_name
	tree$modelfinder = FALSE
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
	return(tree)
}


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

generate_constraint_trees = function(seq_matrix){
	
	clades = unique(seq_matrix@taxon_map)
	
	# Write constraint trees for hypothesis testing
	
	# Ctenophora+Placozoa+Bilateria+Cnidaria, the clade that exists under Proifera-sister
	CtPlBiCn = rownames(seq_matrix)[ seq_matrix@taxon_map %in% c("Ctenophora", "Placozoa", "Bilateria", "Cnidaria") ]
	CtPlBiCn_not = rownames(seq_matrix)[ ! rownames(seq_matrix) %in% CtPlBiCn ]
	porifera_sister_constraint_tree = generate_constaint_tree( CtPlBiCn, CtPlBiCn_not )
	write.tree( 
		porifera_sister_constraint_tree, 
		file = paste( constraint_tree_path, seq_matrix@matrix_name, ".porifera_sister_constraint.tree", sep="") 
	)
	
	# Porifera+Placozoa+Bilateria+Cnidaria, the clade that exists under Ctenophora-sister
	PoPlBiCn = rownames(seq_matrix)[ seq_matrix@taxon_map %in% c("Porifera", "Placozoa", "Bilateria", "Cnidaria") ]
	PoPlBiCn_not = rownames(seq_matrix)[ ! rownames(seq_matrix) %in% PoPlBiCn ]
	ctenophora_sister_constraint_tree = generate_constaint_tree( PoPlBiCn, PoPlBiCn_not )
	write.tree( 
		ctenophora_sister_constraint_tree, 
		file = paste( constraint_tree_path, seq_matrix@matrix_name, ".ctenophora_sister_constraint.tree", sep="") 
	)	
	
	# Ctenophora+Cnidaria, the clade that exists under Coelenterata
	CtCn = rownames(seq_matrix)[ seq_matrix@taxon_map %in% c("Ctenophora", "Cnidaria") ]
	CtCn_not = rownames(seq_matrix)[ ! rownames(seq_matrix) %in% CtCn ]
	coelenterata_constraint_tree = generate_constaint_tree( CtCn, CtCn_not )
	write.tree( 
		coelenterata_constraint_tree, 
		file = paste( constraint_tree_path, seq_matrix@matrix_name, ".coelenterata_constraint.tree", sep="") 
	)		
	
	
}

#' Get the occupancy of a matrix by species and martition
#'
#' @param seq_matrix A PartitionedMultipleAlignment object
#' @return A numeric matrix with rows species, columns partitions, and value the number of sampled amino acids
#' @export
get_matrix_occupancy = function( seq_matrix ){
	n_species = nrow( seq_matrix )
	n_partitions = nrow( seq_matrix@partitions )
	
	seq_matrix_raw = as.matrix(seq_matrix)
	
	occupancy_aa = matrix( NA, nrow=n_species, ncol=n_partitions )
	for( i in 1:n_partitions){
		start = seq_matrix@partitions$start[i]
		stop  = seq_matrix@partitions$stop[i]
		sub_seq_matrix = seq_matrix_raw[1:n_species, start:stop]
		not_gap_matrix = sub_seq_matrix != "-"
		occupied = rowSums( not_gap_matrix )
		occupancy_aa[,i] = occupied
	}

	colnames( occupancy_aa ) = seq_matrix@partitions$partition
	
	rownames( occupancy_aa ) = rownames(seq_matrix)
	
	return( occupancy_aa )
}