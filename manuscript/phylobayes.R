library(ape)


setOldClass("phylo")

#######################################
## Classes


#' An S4 class to store a single sample from a phylobayes chain file
#' In an unconstrained CAT run, the number of distinct categories, and 
#' therefore rows in the frequencies matrix, is 5000, but only a subset 
#' of them are allocated in the allocation vector
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
    branchalpha = "numeric", # value 1 appears constant between iterations
    branchbeta = "numeric", # could be NSPR kappa or kmax, not in any other file
    alpha = "numeric", # appears to be alpha parameter, found in trace file
    ncat = "numeric", # ncat?
    rr = "numeric", # present in CAT runs, but not NCAT
    dirweight = "numeric", # 20 numbers, possibly some starting values for each amino acid (last tab is blank)
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
  
  # Chain files from nCAT and CAT runs have slightly different formats. Specifically, 
  # nCAT files have a tree line, then 4 lines with one scalar each, then a line with a vector.
  # CAT files are the same, but have 5 lines with one scalar each
  text_6 = NA
  text_7 = NA
  freq_start = NA
  
  if( text[7] == "" ){
    text_7 =text[6]
    freq_start = 8
  } else if ( text[8] == ""  ) {
    text_6 = text[6]
    text_7 = text[7]
    freq_start = 9
  } else {
    # The seventh or eight line should be blank
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
    branchalpha = as.numeric(text[2]),
    branchbeta = as.numeric(text[3]),
    alpha = as.numeric(text[4]),
    alpha = as.numeric(text[5]),
    ncat = as.numeric(text_6),
    dirweight = text_7 %>% str_split("\t") %>% unlist() %>% as.numeric(),
    frequencies = text[freq_start:(length(text)-1)] %>% str_split("\t") %>% lapply( as.numeric ) %>% do.call( rbind, . ), 
    allocation = text[length(text)] %>% str_split("\t") %>% unlist() %>% as.numeric()
  )
  
  rownames( object@frequencies ) = 0:( nrow( object@frequencies ) - 1 )
  
  if( object@ncat != nrow(object@frequencies) ){
     # Unexpected number of category frequency profiles
    return( NULL )
  }
  
  # The order of the amino acid vector in PhyloBayes is specified here:
  # https://github.com/bayesiancook/pbmpi/blob/ba03ab5146140c9beab673735e826a10ac0ea6ab/sources/BiologicalSequences.h#L41
  # The trailing gap character is then removed here:
  # https://github.com/bayesiancook/pbmpi/blob/ba03ab5146140c9beab673735e826a10ac0ea6ab/sources/StateSpace.cpp#L131
  
  colnames( object@frequencies ) = 
    c( "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )

  
  object
}

#######################################
## Functions

#' Parse a phylobayes mpi chain file. Incomplete records at the start or end of the file 
#' are ignored
#'
#' @param fine_name name of file to parse
#' @return A PhylobayesSample object
#' @export

parse_phylobayes_chain = function( file_name ){
  
  # Since many files are quite large, it may be apprpriate to tail or head them,
  # in which case the first or last record may be incomplete. Those will be NULL records,
  # which append does not add to a list
  
  
  generation = 0
  samples = list()
  
  # This will be very inefficient for large files, but will work for a first pass
  lines = read_lines( file_name )
  n = 0
  sample_text = character()
  
  # Each record starts with a tree line, which starts with a (

  
  for (line in lines) {
    n = n+1
    if ( grepl( "\\(", line ) && n != 1 ){
      
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
  
  # Add the final record
  samples = append( 
    samples,
    PhylobayesSample( sample_text, file_name, generation )
  )
  
  samples
}

