library(tidyverse)
library(testthat)
source('/gpfs/afm/moxon/thomas2/APAtrap/sub_snakemake/APAtrap/average_bedgraph.R')

merged_bedgraph = tibble(chromosome=c(1,1,1,2,2,2,2),
			start=c(50,80,110,30,60,90,150),
			stop=c(60,90,120,40,70,100,160),
			cov1=c(1,2,5,6,3,4,6),
			cov2=c(1,10,3,6,4,5,6)
)

avg_bedgraph = AvgBedgraph(merged_bedgraph)
print(avg_bedgraph)

context('test that average values are correct')

test_that("avg is a true mean average of 2 coverage columns", {
  expect_equal(avg_bedgraph$avg[1], 1)
  expect_equal(avg_bedgraph$avg[2], 6)
  expect_equal(avg_bedgraph$avg[3], 4)
  expect_equal(avg_bedgraph$avg[4], 6)
  expect_equal(avg_bedgraph$avg[5], 3.5)
  expect_equal(avg_bedgraph$avg[6], 4.5)
  expect_equal(avg_bedgraph$avg[7], 6)
})

merged_bedgraph2 = tibble(chromosome=c(1,1,1,2,2,2,2),
                        start=c(50,80,110,30,60,90,150),
                        stop=c(60,90,120,40,70,100,160),
                        cov1=c(1,2,5,6,3,4,6),
                        cov2=c(1,10,3,6,4,5,6),
			cov3=c(1,10,3,6,4,5,6)
)

avg_bedgraph2 = AvgBedgraph(merged_bedgraph2)

test_that("avg is a true mean average of 3 coverage columns", {
  expect_equal(avg_bedgraph2$avg[1], 1)
  expect_equal(avg_bedgraph2$avg[2], 7.33, tolerance=0.1)
  expect_equal(avg_bedgraph2$avg[3], 3.67, tolerance=0.1)
  expect_equal(avg_bedgraph2$avg[4], 6)
  expect_equal(avg_bedgraph2$avg[5], 3.67, tolerance=0.1)
  expect_equal(avg_bedgraph2$avg[6], 4.67, tolerance=0.1)
  expect_equal(avg_bedgraph2$avg[7], 6)
})

merged_bedgraph3 = tibble(chromosome=c(1,1,1,2,2,2,2),
                        start=c(50,80,110,30,60,90,150),
                        stop=c(60,90,120,40,70,100,160),
                        cov1=c(1,2,5,6,3,4,6),
)

avg_bedgraph3 = AvgBedgraph(merged_bedgraph3)

test_that("avg is a true mean average of 1 coverage columns", {
  expect_equal(avg_bedgraph3$avg[1], 1)
  expect_equal(avg_bedgraph3$avg[2], 2)
  expect_equal(avg_bedgraph3$avg[3], 5)
  expect_equal(avg_bedgraph3$avg[4], 6)
  expect_equal(avg_bedgraph3$avg[5], 3)
  expect_equal(avg_bedgraph3$avg[6], 4)
  expect_equal(avg_bedgraph3$avg[7], 6)
})
