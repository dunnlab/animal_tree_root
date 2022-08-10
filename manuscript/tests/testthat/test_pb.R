library( ape )
library( tidyverse )
library( magrittr )
library( stringr )


source( "../../phylobayes.R" )

context("testing framework")

test_that("Testing works", {
  expect_equal( 1, 1 )
})


# example.chain was created by manually editing 
# Whelan2017_strict.phy_Poisson_CAT60_Chain_1_last500.chain
# to obtain the last seven generations of the chain. This is 
# an nCAT=60 chain.
#
# example_single.chain was created by taking a single record
# from example.chain
#
# example_cat_single.chain was created by taking a single 
# record from Whelan2017_strict.phy_Chain1_last50000.chain.
# It is an unconstrained CAT run



context("single nCAT phylobayes chain record")

lines = read_lines( "example_single.chain" )

pb = PhylobayesSample( lines, "example_single.chain", 1 )


test_that("Tree has the correct number of tips", {
  expect_equal( length(pb@tree$tip.label), 76 )
})

test_that("Tree has the correct number of internal nodes", {
  expect_equal( pb@tree$Nnode, 74 )
})


test_that("Missing tree returns NULL", {
  expect_null( PhylobayesSample( lines[-1], "example_single.chain", 1 )  )
})

test_that("Missing blank line returns NULL", {
  expect_null( PhylobayesSample( lines[-7], "example_single.chain", 1 )  )
})

test_that("Missing last line returns NULL", {
  expect_null( PhylobayesSample( lines[-length(lines)], "example_single.chain", 1 )  )
})




context("phylobayes nCAT chain files")

chain = parse_phylobayes_chain("example.chain")

test_that("Chain has correct number of samples", {
  expect_equal( length(chain), 7 )
})

test_that("Allocation vector lengths are correct", {
  x = lapply(chain, function(x){length(x@allocation)}) %>% unlist()
  expect_true( all(x == 49388) )
})

test_that("Number of tips in trees are correct", {
  x = lapply(chain, function(x){length( x@tree$tip.label )}) %>% unlist()
  expect_true( all(x == 76) )
})


context("single CAT phylobayes chain record")

lines = read_lines( "example_cat_single.chain" )

pb = PhylobayesSample( lines, "example_single.chain", 1 )

test_that("Tree has the correct number of tips", {
  expect_equal( length(pb@tree$tip.label), 76 )
})
