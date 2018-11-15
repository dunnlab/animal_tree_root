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
PartitionedMultipleAlignment = function( alignment_file, partition_file=NULL, taxon_map_global=NULL, partition_map_global=NULL, manuscript_name=NULL, matrix_name=NULL, alignment_format="phylip" ) {
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
			filter( manuscript == manuscript_name ) %>%
			filter( matrix == matrix_name )
		
		
			
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





