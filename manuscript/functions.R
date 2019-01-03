library(Biostrings)
library(tidyverse)


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
	lines = lines[ grepl("CHARSET", lines) ]
	close( conn )
	
	partitions = 
		lapply( 
			lines, 
			function( line ){
				fields = str_match(line, "CHARSET (\\w+)\\s*=\\s*(\\d+)\\s*-\\s*(\\d+)")
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
		
		# Add names to the partitions that are not assigned to genes. Take these from the original 
		# partition names, but prepend matrix name to be sure they are unique across matrices
		
		partitions$gene[ is.na(partitions$gene) ] = paste( matrix_name, partitions$partition[ is.na(partitions$gene) ], collapse = '_' )
		
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
