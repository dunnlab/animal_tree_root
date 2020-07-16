library( ape )
library( tidyverse )
library( magrittr )
library( stringr )


source( "../../phylobayes.R" )

context("testing framework")

test_that("Testing works", {
  expect_equal( 1, 1 )
})

context("phylobayes chain files")