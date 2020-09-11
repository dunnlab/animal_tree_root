library(ape)
library(tidyverse)


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
    ncat = as.numeric(text[5]),
    rr = as.numeric(text_6),
    dirweight = text_7 %>% str_split("\t") %>% unlist() %>% as.numeric(),
    frequencies = text[freq_start:(length(text)-1)] %>% str_split("\t", 20) %>% lapply( as.numeric ) %>% do.call( rbind, . ), 
    allocation = text[length(text)] %>% str_trim() %>% str_split("\t") %>% unlist() %>% as.numeric()
  )
  
  rownames( object@frequencies ) = 0:( nrow( object@frequencies ) - 1 )
  
  # The order of the amino acid vector in PhyloBayes is specified here:
  # https://github.com/bayesiancook/pbmpi/blob/ba03ab5146140c9beab673735e826a10ac0ea6ab/sources/BiologicalSequences.h#L41
  # The trailing gap character is then removed here:
  # https://github.com/bayesiancook/pbmpi/blob/ba03ab5146140c9beab673735e826a10ac0ea6ab/sources/StateSpace.cpp#L131
  
  colnames( object@frequencies ) = 
    c( "A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y" )
  
  if( object@ncat != nrow(object@frequencies) ){
     # Unexpected number of category frequency profiles
    return( NULL )
  }
  
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


parse_phylobayes_last_sample = function( file_name ){
  chain = parse_phylobayes_chain( file_name )
  chain[[length(chain)]]
}


#' Summarise each category
#'
#' @param pb A PhylobayesSample object
#' @return A tibble, with one row per equilibrium frequency category
#' @export

summarise_sample = function( pb ){
  
  allocation_counts = table(pb@allocation)
  index = allocation_counts %>% names() %>% as.numeric()
  
  if( ! all( index %in% rownames(pb@frequencies) ) ){
    stop("Allocations are outside available category range")
  }
  
  # Consider only the subset of frequencies that are allocated
  frequencies = pb@frequencies[ match(index , rownames(pb@frequencies) ), ]
  
  # Building on https://www.statmethods.net/advstats/mds.html
  d = dist( frequencies )
  fit = cmdscale(d,eig=TRUE, k=2)
  
  categories = 
    tibble(
      index = index,
      count = c(allocation_counts),
      x = fit$points[,1],
      y = fit$points[,2]
    ) 
  
  if( ! all( categories$index == rownames(frequencies) ) ){
    stop("category indexes and frequency vector names do not match")
  }
  
  frequencies_tib = 
    as_tibble( frequencies )
  
  names( frequencies_tib ) = paste0( "aa.", names( frequencies_tib ) )
  
  categories %<>% 
    bind_cols( frequencies_tib ) %>%
    arrange(desc(count)) %>% 
    mutate( cumsum = cumsum(count)) %<>% 
    mutate(rank = 1:nrow(categories))
  
}

# Return the allocation midpoint, ie the percent of sites that are allocated to the
# 50% of the least frequent categories
allocation_midpoint = function ( pb_summary ) {
  x = pb_summary %>% pull( count ) %>% sort(decreasing = FALSE)
  x = x/ sum(x)
  y = cumsum( x )
  return ( round( median( y ) *100, 2 ) )
}
