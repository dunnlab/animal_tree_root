library(Biostrings)
library(tidyverse)

#' An S4 class to store multiple sequence alignments with partition data
#' 
#' @slot partitions  The partitions
setClass(
	Class = "PartitionedMultipleAlignment",
	representation = representation(
		partitions = "data.frame"
	),
	contains = "MultipleAlignment"
)

#' Construct a PartitionedMultipleAlignment from an alignment file and a partition file
#'
#' @param data_list A list containing the expression data
#' @return An Expression object
#' @export
PartitionedMultipleAlignment = function( alignment_file, partition_file=NULL, alignment_format="phylip" ) {
	object = readAAMultipleAlignment( filepath = alignment_file, format=alignment_format )
	class(object) = "PartitionedMultipleAlignment"
	
	if( ! is.null(partition_file) ){
		
		conn = file( partition_file, open="r")
		lines = readLines( conn )
		lines = lines[ grepl("CHARSET", lines) ]
		
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
		
	}
	
	object
}



a_file = "~/Desktop/animal_root/blast/phylips/Chang2015_Chang2015.phy"
p_file = "~/Desktop/animal_root/blast/phylips/Chang2015_Chang2015.nex"
msa = PartitionedMultipleAlignment( a_file, p_file )
rownames( msa )
