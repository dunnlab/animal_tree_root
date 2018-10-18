library(Biostrings)
library(tidyverse)

#' An S4 class to store multiple sequence alignments with partition data
#' 
#' @slot partitions  The partitions
setClass(
	Class = "PartitionedMultipleAlignment",
	representation = representation(
		partitions = "data.frame",
		taxon_map = "vector",
		partition_map = "vector"
	),
	contains = "MultipleAlignment"
)

#' Construct a PartitionedMultipleAlignment from an alignment file and a partition file
#'
#' @param alignment_file Path to the gene alignment
#' @param partition_file Path to the partition file in nexus format
#' @param taxon_map A dataframe where the first column corresponds to sequence names in the alignment and the second column to clade names for the sequences
#' @param partition_map A dataframe where the first column corresponds to partition names and the second column to names for groups of partitions
#' @return A PartitionedMultipleAlignment object
#' @export
PartitionedMultipleAlignment = function( alignment_file, partition_file=NULL, taxon_map=NULL, partition_map=NULL, alignment_format="phylip" ) {
	object = readAAMultipleAlignment( filepath = alignment_file, format=alignment_format )
	class(object) = "PartitionedMultipleAlignment"
	
	
		
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
	
	object@partitions = partitions
	
	if( ! is.null(taxon_map) ){
		sequence_name = rownames( object )
		object@taxon_map = taxon_map[ match( sequence_name, taxon_map[[1]] ), 2 ][[1]]
	}
	
	object
}



