library(ape)


setOldClass("phylo")

#######################################
## Classes


#' An S4 class to store a single sample from a phylobayes chain file
#' 
#' @slot allocation  The allocation vector
setClass(
  "PhylobayesSample",
  slots = c(
    
    path = "character",
    generation = "numeric",
    
    # These are shown in the order of a phylobayes MPI record. One item per line
    # unless otherwise noted
    
    tree = "phylo",
    x2 = "numeric", # value 1 appears constant between iterations
    x3 = "numeric", # could be NSPR kappa or kmax, not in any other file
    x4 = "numeric", # appears to be alpha parameter, found in trace file
    x5 = "numeric", # ncat?
    x6 = "numeric", # 20 numbers, possibly some starting values for each amino acid (last tab is blank)
    # Blank line
    frequencies = "matrix", # number of columns is number of states (eg 20 for amino acids), number of rows is the number of categories (one line per category in the file)
    allocation = "numeric" # integer that indicates the category of each site
  )
)

#' Construct a PhylobayesSample from the text for a single generation from a phylobayes chain
#'
#' @param sample_text text for a single sample from a phylobayes chain file
#' @param path path to the original file
#' @param generation the generation of the sample in the original file
#' @return A PhylobayesSample object, or NULL if the record is incomplete
#' @export
PhylobayesSample = function( 
  text, 
  path = NA,
  generation = NA
  ) {
  
  
  if( text[7] != "" ){
    # The seventh line should be blank
    return( NULL )
  }
  
  if ( ! grepl( "\\(", text[1] ) ) {
    # Record does not start with a tree
    return( NULL )
  }
  
  
  object = new( 
    "PhylobayesSample",
    path = path, 
    generation = generation,
    tree = ape::read.tree( text=text[1] ),
    x2 = as.numeric(text[2]),
    x3 = as.numeric(text[3]),
    x4 = as.numeric(text[4]),
    x5 = as.numeric(text[5]),
    x6 = text[6] %>% str_split("\t") %>% unlist() %>% as.numeric(),
    frequencies = text[8:(length(text)-1)] %>% str_split("\t") %>% lapply( as.numeric ) %>% do.call( rbind, . ), 
    allocation = text[length(text)] %>% str_split("\t") %>% unlist() %>% as.numeric()
  )
  
  if( object@x5 != nrow(object@frequencies) ){
     # Unexpected number of category frequency profiles
    return( NULL )
  }
  
  object
}

#######################################
## Functions

#' Parse a phylobayes mpi chain file
#'
#' @param fine_name name of file to parse
#' @return A PhylobayesSample object
#' @export

parse_phylobayes_chain = function( file_name ){
  
  generation = 0
  samples = list()
  
  # This will be very inefficient for large files, but will work for a first pass
  lines = read_lines( file_name )
  n = 0
  sample_text = character()
  
  # Each record starts with a tree line, which starts with a (
  # Since many files are quite large, it may be apprpriate to tail them,
  # in which case the first record may be incomplete. Set a flag to note when the
  # first complete record has started, and don't parse before then.
  start = FALSE
  
  for (line in lines) {
    n = n+1
    if( n == 1 ) {
      if ( grepl( "\\(", line ) ) {
        generation = generation + 1
        sample_text = append( sample_text, line )
      } else {
        stop( "Chain file does not start with a tree" )
      }
    } else {
      if ( grepl( "\\(", line ) ){
        
        samples = append( 
          samples,
          PhylobayesSample( sample_text, file_name, generation )
          )
        
        # Get ready to parse new record, starting with this line
        generation = generation + 1
        sample_text = character()
      }
      sample_text = append( sample_text, line )
    }
  }
  
  # Add the final record
  samples = append( 
    samples,
    PhylobayesSample( sample_text, file_name, generation )
  )
  
  samples
}

