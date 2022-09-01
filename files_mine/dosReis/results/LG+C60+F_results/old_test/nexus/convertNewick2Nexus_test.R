### conda activate r_ape
## Description: This script converts newick to nexus format

#install.packages("ape")
library(ape)

args = commandArgs(trailingOnly=TRUE)
if ( length(args) == 2 ) {
	message("Changing newick tree: ",args[1])
	message(" to nexus format:     ",args[2])
} else {
	stop("You need to give two files: 1. newick tree and 2. out nexus file name.", call.=FALSE)
}

args <- commandArgs(TRUE)
tree <- read.tree(args[1])
outfile <- paste0(make.names(args[2]))
write.nexus(tree, file=outfile)
