library( ape )
library( tidyverse )
library( magrittr )
library( stringr )


source( "../../phylobayes.R" )

context("testing framework")

# The example file was gcreated by manually editing 
# Whelan2017_strict.phy_Poisson_CAT60_Chain_1_last500.chain
# to obtain the last three generations of the chain

example_name = "example.chain"

test_that("Testing works", {
  expect_equal( 1, 1 )
})


context("single phylobayes chain record")

lines = read_lines( "example_single.chain" )

pb = PhylobayesSample( lines, "example_single.chain", 1 )

test_that("Tree has the correct number of tips", {
  expect_equal( pb@tree$Nnode, 74 )
})

context("phylobayes chain files")

chain = parse_phylobayes_chain("example.chain")

test_that("Tree has the correct number of tips", {
  expect_equal( length(chain), 7 )
})

