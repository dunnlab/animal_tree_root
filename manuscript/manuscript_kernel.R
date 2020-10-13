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

source( "functions.R" )
source( "phylobayes.R" )

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
    Borowiec2015 Total1080
    Chang2015 Chang2015
    Dunn2008  Dunn2008
    Hejnol2009 Hejnol2009
    Moroz2014 ED3a
    Nosenko2013 nonribosomal_9187_smatrix
    Nosenko2013 ribosomal_11057_smatrix
    Philippe2009 Philippe2009
    Ryan2013 est.opisthokonta
    Ryan2013 genome.opisthokonta
    Simion2017 supermatrix_97sp_401632pos_1719genes
    Whelan2015 Metazoa_Choano_RCFV_strict",
    header = TRUE,
    stringsAsFactors = FALSE)

bootstrap_threshold = 90
posterior_prob_threshold = 95

# Load data

papers = read_tsv( "../data_processed/tables/previously_published_manuscripts.tsv" )
datasets = read_tsv( "../data_processed/tables/previously_published_matrices.tsv" )
analyses_published =
  read_tsv( "../data_processed/tables/previously_published_analyses.tsv") %>%
  mutate( clade = factor( clade, levels = c( "Choanimalia", "Holozoa", "Opisthokonta" ) ) )

taxonomy_reference = read_tsv("../reconciliation/taxonomy_info/taxon_table.tsv")

# raw columns:
# component_number, matrix, partition_name, edges, nodes_in_component, component_density, BUSCO_ID, BUSCO_description,
# SwissProt_accession, SwissProt_description, GO_annotations, ribo_found
partition_map_global =
  read_tsv("../reconciliation/blast/graphs/partition_components_split_annotated.tsv") %>%
  dplyr::rename(partition = partition_name) %>%
  mutate( component_number = as.character(component_number) )

analyses_published$result = "Unresolved"
analyses_published$result[ 
  ( analyses_published$inference == "Bayesian" ) & 
  (analyses_published$support_porifera_sister >= posterior_prob_threshold) ] = "Porifera-sister"
analyses_published$result[ 
  ( analyses_published$inference == "Bayesian" ) & 
  (analyses_published$support_ctenophora_sister >= posterior_prob_threshold) ] = "Ctenophora-sister"
analyses_published$result[ 
  ( analyses_published$inference == "ML" ) & 
  (analyses_published$support_porifera_sister >= bootstrap_threshold) ] = "Porifera-sister"
analyses_published$result[ 
  ( analyses_published$inference == "ML" ) & 
  (analyses_published$support_ctenophora_sister >= bootstrap_threshold) ] = "Ctenophora-sister"
analyses_published$result = factor( analyses_published$result )

analyses_published$model_combined = 
  factor( 
    analyses_published$model_combined, 
    levels = c( "WAG", "LG", "GTR", "data partitioning", "Recoding + GTR", 
                "Recoding + GTR + CAT", "Poisson + CAT", "GTR + CAT" ) )

# Matrix taxon composition

clades = c( "Fungi", "Ichthyosporea", "Filasterea", "Choanoflagellida", "Ctenophora", 
            "Porifera", "Placozoa", "Bilateria", "Cnidaria" )

taxa =
  taxonomy_reference %>%
  distinct( relabelled_name, clade_assignment, ncbi_tax_id ) %>%
  dplyr::rename( taxon = relabelled_name, clade = clade_assignment ) %>%
  mutate( clade = factor( clade, levels = clades )  )

matrix_path = "../data_processed/matrices"
phylip_file_names = list.files(path = matrix_path, pattern = ".+\\.phy$", full.names = TRUE)

sequence_matrices = foreach( phylip_file = phylip_file_names) %dopar%
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
  summarise( n_sequences = n() )


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
  summarise( BUSCO = names(which(table(Description) == max(table(Description)))[1]) )

partition_map_global %<>% left_join( partition_to_busco_map, by = c("matrix", "partition") )


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
    by = c("matrix")
  )

matrix_summary$manuscript = str_split( matrix_summary$matrix, "_", simplify = TRUE )[, 1]

cluster_summary =
  partition_map_global %>%
  group_by( component_number ) %>%
  summarise(
    n_partitions = n(),
    n_matrices = length(unique(matrix)),
    n_with_busco = sum(! is.na(BUSCO)),
    n_unique_busco = length(unique(na.omit(BUSCO)))
  )


# Matrix overlap
matrix_overlap =
  lapply(sequence_matrices, function(x) lapply(sequence_matrices, function(y) compute_matrix_overlap(x, y))) %>%
  unlist(recursive = FALSE) %>%
  bind_rows()

# Remove reciprocal comparisons and comparisons to self
n = nrow(matrix_overlap)
mask = lower.tri(matrix(nrow = sqrt(n), ncol = sqrt(n)))
dim(mask) = NULL
matrix_overlap = matrix_overlap[mask, ]


# New analyses of published matrices

# read iqtrees
trees_path_iqtree = "../trees_new/iqtree"
iqtree_ext = "\\.treefile$"
file_names_iqtree = list.files( path = trees_path_iqtree, pattern = iqtree_ext, full.names = TRUE )

trees_iq = foreach( tree_file = file_names_iqtree ) %dopar%
  parse_tree_iqtree( tree_file, taxonomy_reference )

# read phylobayes 
trees_path_sensitive = "../trees_new/sensitive"
pb_tree_ext = "\\.con\\.tre$"

# read pb trees
trees_path_pb = "../trees_new/phylobayes"
file_names_pb = list.files( path = trees_path_pb, pattern = pb_tree_ext, full.names = TRUE )

trees_pb = foreach( tree_file = file_names_pb ) %dopar%
           parse_tree_pb( tree_file, taxonomy_reference )

# sensitivity analyses
file_names_sensitive = list.files( path = trees_path_sensitive, pattern = pb_tree_ext, full.names = TRUE )
trees_sensitive = foreach( tree_file = file_names_sensitive ) %dopar%
  parse_tree_pb( tree_file, taxonomy_reference )

# sensitivity tibble
analyses_sensitive = lapply(
  trees_sensitive,
  function( tree ){
    data.frame(
      matrix = tree$matrix,
      model = tree$model,
      clade = tree$sampling,
      support_ctenophora_sister = tree$ctenophora_sister * 100,
      support_porifera_sister = tree$porifera_sister * 100,
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

analyses_sensitive$inference = rep( "Bayesian", length( trees_sensitive ) )
analyses_sensitive$result = "Unresolved"

analyses_sensitive$result[ (analyses_sensitive$support_porifera_sister >= posterior_prob_threshold) ] = 
  "Porifera-sister"

analyses_sensitive$result[ (analyses_sensitive$support_ctenophora_sister >= posterior_prob_threshold) ] = 
  "Ctenophora-sister"

analyses_sensitive$result = factor( analyses_sensitive$result )
analyses_sensitive$model_summary = factor( 
  analyses_sensitive$model,
  levels = c("Poisson+CAT60", "Poisson+CAT70", "Poisson+CAT80", "Poisson+CAT90", "Poisson+CAT100", 
             "Poisson+CAT110", "Poisson+CAT120", "Poisson+CAT150", "Poisson+CAT180", "Poisson+CAT270", 
             "Poisson+CAT340", "Poisson+CAT360", "Poisson+CAT380", "Poisson+CAT400", "Poisson+CAT420", 
             "Poisson+CAT440", "Poisson+CAT460", "Poisson+CAT480")
)

# new trees tibble
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
analyses_new$result[ ( analyses_new$inference == "Bayesian" ) & 
                       (analyses_new$support_porifera_sister >= posterior_prob_threshold) ] = "Porifera-sister"
analyses_new$result[ ( analyses_new$inference == "Bayesian" ) & 
                       (analyses_new$support_ctenophora_sister >= posterior_prob_threshold) ] = "Ctenophora-sister"
analyses_new$result[ ( analyses_new$inference == "ML" ) & 
                       (analyses_new$support_porifera_sister >= bootstrap_threshold) ] = "Porifera-sister"
analyses_new$result[ ( analyses_new$inference == "ML" ) & 
                       (analyses_new$support_ctenophora_sister >= bootstrap_threshold) ] = "Ctenophora-sister"
analyses_new$result = factor( analyses_new$result )

# Parse model components
analyses_new$model_summary = analyses_new$model
analyses_new$model_summary = factor( analyses_new$model_summary, levels = c("WAG", "GTR20", "Poisson+C60", "WAG+C60", 
                                                                          "LG+C60", "CAT+F81"))

# Partition comparison across matrices

n_total_partitions =
  partitions_all %>%
  group_by(matrix) %>%
  summarize("n_total_partitions" = n_distinct(partition))
n_components_with_BUSCO =
  partition_map_global %>%
    filter(BUSCO_ID != "") %>%
    group_by( matrix, component_number ) %>%
    summarize(n()) %>%
    group_by(matrix) %>%
    tally(name = "n_components_with_BUSCO")
n_ribo =
  partition_map_global %>%
    group_by( matrix, component_number ) %>%
    tally(ribo_found) %>%
    group_by(matrix) %>%
    tally()
discarded_parts =
  read_tsv("../reconciliation/blast/graphs/discarded_nodes.tsv") %>%
  group_by(matrix) %>%
  summarize("n_partitions_discarded" = n_distinct(partition_name))

partition_network_summary =
  n_total_partitions %>%
  left_join(n_components_with_BUSCO, by = "matrix") %>%
  left_join(n_ribo,                  by = "matrix") %>%
  left_join(discarded_parts,         by = "matrix") %>%
  mutate(n_partitions_discarded = replace_na(n_partitions_discarded, 0))

### Summarize phylogenetic signal by genes
au_tests = parse_au_gene_tests()

### Summarize categories from pbmpi

# Parse the last sample from chain 1 of each analysis
phil_cat_c1 =
  parse_phylobayes_last_sample("../trees_new/frequency/subsampled_Philippe2009_only_choanozoa.phy_Poisson_CAT_Chain_1.chain")

phil_cat60_c1 =
  parse_phylobayes_last_sample("../trees_new/frequency/subsampled_Philippe2009_only_choanozoa.phy_Poisson_nCAT60_Chain_1.chain")

whel_cat_c1 =
  parse_phylobayes_last_sample("../trees_new/frequency/subsampled_Whelan2017_strict.phy_Poisson_CAT_Chain1.chain")

whel_cat60_c1 =
  parse_phylobayes_last_sample("../trees_new/frequency/subsampled_Whelan2017_strict.phy_Poisson_CAT60_Chain_1.chain")

# Create a single tibble with summaries of all analyses
pb_summaries = list(
  summarise_sample(phil_cat_c1) %>%
    mutate ( chain=1, matrix="Philippe2009", model="Poisson+CAT"), 
  summarise_sample(phil_cat60_c1) %>%
    mutate ( chain=1, matrix="Philippe2009", model="Poisson+nCAT60"), 
  summarise_sample(whel_cat_c1) %>%
    mutate ( chain=1, matrix="Whelan2017_strict", model="Poisson+CAT"), 
  summarise_sample(whel_cat60_c1) %>%
    mutate ( chain=1, matrix="Whelan2017_strict", model="Poisson+nCAT60")
) %>%
  bind_rows() %>%
  mutate( analysis = paste(matrix, model, sep=" ")) %>%
  mutate( model=factor(model, levels=c("Poisson+nCAT60", "Poisson+CAT")) )

pb_frequencies = pb_summaries %>% select( starts_with("aa.") ) %>% data.matrix()

# Identify the midpoint of each set of allocations
mid_phil_cat =
  pb_summaries %>%
  filter( matrix=="Philippe2009" ) %>%
  filter( model=="Poisson+CAT" ) %>%
  allocation_midpoint()

mid_phil_ncat60 =
  pb_summaries %>%
  filter( matrix=="Philippe2009" ) %>%
  filter( model=="Poisson+nCAT60" ) %>%
  allocation_midpoint()

mid_whel_cat =
  pb_summaries %>%
  filter( matrix=="Whelan2017_strict" ) %>%
  filter( model=="Poisson+CAT" ) %>%
  allocation_midpoint()

mid_whel_ncat60 =
  pb_summaries %>%
  filter( matrix=="Whelan2017_strict" ) %>%
  filter( model=="Poisson+nCAT60" ) %>%
  allocation_midpoint()

n_categories_phil_cat =
  pb_summaries %>%
  filter( matrix=="Philippe2009" ) %>%
  filter( model=="Poisson+CAT" ) %>%
  nrow()

n_categories_whel_cat =
  pb_summaries %>%
  filter( matrix=="Whelan2017_strict" ) %>%
  filter( model=="Poisson+CAT" ) %>%
  nrow()


# Perform global MDS analysys
fit = cmdscale( dist( pb_frequencies ) ,eig=TRUE, k=2)
pb_summaries %<>% mutate( x_global = fit$points[,1], y_global = fit$points[,2] )

## Record information about the session
session_info_kernel = sessionInfo()
system_time_kernel = Sys.time()

commit_kernel =
  system("git log | head -n 1", intern = TRUE) %>%
  str_replace("commit ", "")

time_stop = Sys.time()
time_run = time_stop - time_start

## Write the results to prepare them for manuscript.rmd

save.image("manuscript.RData")

## Write the results from RData to prepare them for supplementary tables
write_csv(partition_map_global, "./Supplementary_tables/Supplementary_Table_1.csv", na = "NA", quote_escape = "double")
write_csv(analyses_published, "./Supplementary_tables/Supplementary_Table_2.csv", na = "NA", quote_escape = "double")
write_csv(analyses_new, "./Supplementary_tables/Supplementary_Table_3.csv", na = "NA", quote_escape = "double")
write_csv(analyses_sensitive, "./Supplementary_tables/Supplementary_Table_6.csv", na = "NA", quote_escape = "double")
write_csv(au_tests, "./Supplementary_tables/Supplementary_Table_7.csv", na = "NA", quote_escape = "double")

## Move all read functions to kernel so all vairables are stored in Rdata

cat_categories = read_tsv("../data_processed/tables/cat_categories.tsv")
table_study_summary = read_tsv("../data_processed/tables/study_summary.tsv")
taxa_map_whelan=read_tsv("../reconciliation/taxonomy_info/taxon_table.tsv") %>% filter(original_matrix=="../considered_data/Whelan2017/strict.phy")
taxa_map_philippe=read_tsv("../reconciliation/taxonomy_info/taxon_table.tsv") %>% filter(original_matrix=="../considered_data/Philippe2009/Philippe2009.nex")
tree1= read.tree("../trees_new/iqtree/Whelan2017_strict.model_test.treefile")
tree2=read.tree("../trees_new/sensitive/Whelan2017_strict.phy_Poisson_CAT60.con.tre")
tree3=read.tree("../trees_new/sensitive/Whelan2017_strict.phy_Poisson_CAT90.con.tre")
tree4=read.tree("../trees_new/phylobayes/Whelan2017_strict.bpcomp.con.tre")
tree5=read.tree("../trees_new/phylobayes/Whelan2017_strict.phy_GTR_CAT.con.tre")
tree6= read.tree("../trees_new/iqtree/Philippe2009_only_choanozoa.WAG+C60.treefile")
tree7=read.tree("../trees_new/sensitive/Philippe2009_only_choanozoa.phy_Poisson_CAT60.con.tre")
tree8=read.tree("../trees_new/sensitive/Philippe2009_only_choanozoa.phy_Poisson_CAT150.con.tre")
tree9=read.tree("../trees_new/phylobayes/Philippe2009_only_choanozoa.bpcomp.con.tre")
tree10=read.tree("../trees_new/phylobayes/Philippe2009_only_choanozoa.phy_GTR_CAT.con.tre")
ribo = read_tsv("../data_processed/tables/ribosomal_gene.tsv")
busco = read_tsv("../data_processed/tables/busco_gene.tsv")
cross_validation =
  read_tsv("../data_processed/tables/cross_validation.tsv") %>%
  gather("model", "score", `nCAT60`, `Poisson-CAT`)



