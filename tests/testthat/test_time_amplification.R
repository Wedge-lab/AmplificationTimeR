# Test time_amplification.R
library(testthat)
library(covr)

# test input
test_that("AmplificationTimeR runs and produces a data.frame when all input is correct", {
  expect_equal(class(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status)),"data.frame")
})

# test copy number input
test_that("AmplificationTimeR fails and produces an error message when copy number input is incorrect class.", {
  expect_error(time_amplification(cn_data = as.matrix(test_data_cn),
                                    multiplicity_data = test_data_mult,
                                    sample_id = test_data_id,
                                    amplification_chrom = test_data_chrom,
                                    amplification_start = test_data_start,
                                    amplification_stop = test_data_stop,
                                    is_WGD = test_data_status))
})

test_that("AmplificationTimeR fails and produces an error message when copy number input has incorrect columns.", {
  expect_error(time_amplification(cn_data = test_data_cn_incorrect,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
test_that("AmplificationTimeR fails and produces an error message when copy number input is empty.", {
  expect_error(time_amplification(cn_data = test_data_cn_empty,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
# test multiplicity input
test_that("AmplificationTimeR fails and produces an error message when multiplicity input is incorrect class.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = as.matrix(test_data_mult),
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})

test_that("AmplificationTimeR fails and produces an error message when multiplicity input has incorrect columns.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult_incorrect,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
test_that("AmplificationTimeR fails and produces an error message when multiplicity input is empty.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult_empty,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
# test chromosome input
test_that("AmplificationTimeR fails and produces an error message when chromosome input is incorrect.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = 30,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
# test start input
test_that("AmplificationTimeR fails and produces an error message when start input is incorrect or not found in file.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = 1000000000,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
  
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = "one",
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
# test stop input - add in a test informed by cytobands to make sure people don't specify something outrageous?
test_that("AmplificationTimeR fails and produces an error message when stop input is incorrect.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = "one",
                                  is_WGD = test_data_status))
})
# test WGD input
test_that("AmplificationTimeR fails and produces an error message when is_WGD input is incorrect.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = "yes"))
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = 1))
  
})
# test name input
test_that("AmplificationTimeR fails and produces an error message when sample_id input is incorrect.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = 1,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
  
  
})

# test that code deals with nMaj2_A AND nMin2_A in input

test_that("AmplificationTimeR runs and produces a data.frame when nMaj2_A and nMin2_A are included", {
  expect_equal(class(time_amplification(cn_data = test_data_cn_n2,
                                        multiplicity_data = test_data_mult,
                                        sample_id = test_data_id,
                                        amplification_chrom = test_data_chrom,
                                        amplification_start = test_data_start,
                                        amplification_stop = test_data_stop,
                                        is_WGD = test_data_status)),"data.frame")
})

# test that code fails if region supplied is not amplified

test_that("AmplificationTimeR fails when region is not amplified in diploid sample.", {
  expect_error(time_amplification(cn_data = test_data_cn_not_amp,
                                        multiplicity_data = test_data_mult,
                                        sample_id = test_data_id,
                                        amplification_chrom = test_data_chrom,
                                        amplification_start = test_data_start,
                                        amplification_stop = test_data_stop,
                                        is_WGD = test_data_status))
})

# test multiplicity has mutations in region

test_that("AmplificationTimeR fails when there are no mutations in the amplified region", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult_no_muts_in_region,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})

# test that code fails if region supplied is not amplified with n2A present in diploid sample
test_that("AmplificationTimeR fails when region is not amplified in diploid sample with n2A.", {
  expect_error(time_amplification(cn_data = test_data_cn_not_amp,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})

# test that code fails if region supplied is not amplified - WGD
test_that("AmplificationTimeR fails when region is not amplified in WGD sample.", {
  expect_error(time_amplification(cn_data = test_data_cn_not_amp_wgd,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE))
})

# test that code fails if region supplied is not amplified with n2A present in WGD sample
test_that("AmplificationTimeR fails when region is not amplified in WGD sample with n2A.", {
  expect_error(time_amplification(cn_data = test_data_cn_not_amp_wgd_n2A,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE))
})

# test that n2A_sum will be chosen if it is higher than n1A_sum
test_that("AmplificationTimeR runs and produces a data.frame when nMaj2_A and nMin2_A are included and n2 has the higher CN value", {
  expect_equal(time_amplification(cn_data = test_data_cn_n2_higher,
                                        multiplicity_data = test_data_mult,
                                        sample_id = test_data_id,
                                        amplification_chrom = test_data_chrom,
                                        amplification_start = test_data_start,
                                        amplification_stop = test_data_stop,
                                        is_WGD = test_data_status)[,"highest_copy_number"],"2+1")
})

# test output
test_that("AmplificationTimeR runs and produces a data.frame with correct columns when all input is correct", {
                                           expect_equal(colnames(time_amplification(cn_data = test_data_cn,
                                                                                    multiplicity_data = test_data_mult,
                                                                                    sample_id = test_data_id,
                                                                                    amplification_chrom = test_data_chrom,
                                                                                    amplification_start = test_data_start,
                                                                                    amplification_stop = test_data_stop,
                                                                                    is_WGD = test_data_status)),c("sample","region","highest_copy_number","event_order",
                                                                                                                  "t_1","t_1_lower_ci","t_1_upper_ci",
                                                                                                                  "t_2","t_2_lower_ci","t_2_upper_ci",
                                                                                                                  "t_3","t_3_lower_ci","t_3_upper_ci",
                                                                                                                  "t_4","t_4_lower_ci","t_4_upper_ci",
                                                                                                                  "t_5","t_5_lower_ci","t_5_upper_ci",
                                                                                                                  "t_6","t_6_lower_ci","t_6_upper_ci",
                                                                                                                  "t_7","t_7_lower_ci","t_7_upper_ci",
                                                                                                                  "t_8","t_8_lower_ci","t_8_upper_ci",
                                                                                                                  "t_9","t_9_lower_ci","t_9_upper_ci",
                                                                                                                  "t_10","t_10_lower_ci","t_10_upper_ci"))
                                         })

# *** Issue with this test is that bootstrap isn't reproducible because of seed - fix in later versions ***
test_that("AmplificationTimeR produces expected output when all input is correct", {
  expect_equal(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status), test_output)
})

########
# create tests for different copy numbers and test output
###
# 2+1
test_that("AmplificationTimeR runs and produces 2+2 highest output", {
  expect_equal(time_amplification(cn_data = test_data_cn_2_2,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status)[,"highest_copy_number"],"2+2")
})

# 3+1
test_that("AmplificationTimeR runs and produces 3+1 highest output", {
  expect_equal(time_amplification(cn_data = test_data_cn_3_1,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status)[,"highest_copy_number"],"3+1")
})


# 3+2
test_that("AmplificationTimeR runs and produces 3+2 highest output", {
  expect_equal(time_amplification(cn_data = test_data_cn_3_2,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status)[,"highest_copy_number"],"3+2")
})


# 4+1
test_that("AmplificationTimeR runs and produces 4+1 highest output", {
  expect_equal(time_amplification(cn_data = test_data_cn_4_1,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status)[,"highest_copy_number"],"4+1")
})

# 4+2
test_that("AmplificationTimeR runs and produces 4+2 highest output and right order", {
  expect_equal(time_amplification(cn_data = test_data_cn_4_2,
                                  multiplicity_data = test_data_mult_4_2_WGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number",],"4+2")
  expect_equal(time_amplification(cn_data = test_data_cn_4_2,
                                  multiplicity_data = test_data_mult_4_2_WGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order",],"WGG") # fix once added
  expect_equal(time_amplification(cn_data = test_data_cn_4_2,
                                  multiplicity_data = test_data_mult_4_2_GW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number",],"4+2")
  expect_equal(time_amplification(cn_data = test_data_cn_4_2,
                                  multiplicity_data = test_data_mult_4_2_GW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order",],"GW") # fix once added
})


# 4+3
test_that("AmplificationTimeR runs and produces 4+3 highest output", {
  expect_equal(time_amplification(cn_data = test_data_cn_4_3,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status)[,"highest_copy_number"],"4+3") # fix once added
})

# 4+4
test_that("AmplificationTimeR runs and produces 4+4 highest output", {
  expect_equal(time_amplification(cn_data = test_data_cn_4_4,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status)[,"highest_copy_number"],"4+4")
})

# 5+1
test_that("AmplificationTimeR runs and produces 5+1 highest output and order", {
  expect_equal(time_amplification(cn_data = test_data_cn_5_1,
                                  multiplicity_data = test_data_mult_5_1_WGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"5+1")
  expect_equal(time_amplification(cn_data = test_data_cn_5_1,
                                  multiplicity_data = test_data_mult_5_1_WGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"WGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_5_1,
                                  multiplicity_data = test_data_mult_5_1_GWG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"5+1")
  expect_equal(time_amplification(cn_data = test_data_cn_5_1,
                                  multiplicity_data = test_data_mult_5_1_GWG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GWG")
})

# 5+2
test_that("AmplificationTimeR runs and produces 5+2 highest output and order", {
  expect_equal(time_amplification(cn_data = test_data_cn_5_2,
                                  multiplicity_data = test_data_mult_5_2_WGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"5+2")
  expect_equal(time_amplification(cn_data = test_data_cn_5_2,
                                  multiplicity_data = test_data_mult_5_2_WGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"WGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_5_2,
                                  multiplicity_data = test_data_mult_5_2_GWG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"5+2")
  expect_equal(time_amplification(cn_data = test_data_cn_5_2,
                                  multiplicity_data = test_data_mult_5_2_GWG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GWG")
})

# 6+1
test_that("AmplificationTimeR runs and produces 6+1 highest output and order", {
  expect_equal(time_amplification(cn_data = test_data_cn_6_1,
                                  multiplicity_data = test_data_mult_6_1_WGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"6+1")
  expect_equal(time_amplification(cn_data = test_data_cn_6_1,
                                  multiplicity_data = test_data_mult_6_1_WGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"WGGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_6_1,
                                  multiplicity_data = test_data_mult_6_1_GWGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"6+1")
  expect_equal(time_amplification(cn_data = test_data_cn_6_1,
                                  multiplicity_data = test_data_mult_6_1_GWGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GWGG")
  expect_equal(time_amplification(cn_data = test_data_cn_6_1,
                                  multiplicity_data = test_data_mult_6_1_GGW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"6+1")
  expect_equal(time_amplification(cn_data = test_data_cn_6_1,
                                  multiplicity_data = test_data_mult_6_1_GGW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GGW")
})

# 6+2
test_that("AmplificationTimeR runs and produces 6+2 highest output and order", {
  expect_equal(time_amplification(cn_data = test_data_cn_6_2,
                                  multiplicity_data = test_data_mult_6_2_WGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"6+2")
  expect_equal(time_amplification(cn_data = test_data_cn_6_2,
                                  multiplicity_data = test_data_mult_6_2_WGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"WGGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_6_2,
                                  multiplicity_data = test_data_mult_6_2_GWGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"6+2")
  expect_equal(time_amplification(cn_data = test_data_cn_6_2,
                                  multiplicity_data = test_data_mult_6_2_GWGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GWGG")
  expect_equal(time_amplification(cn_data = test_data_cn_6_2,
                                  multiplicity_data = test_data_mult_6_2_GGW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"6+2")
  expect_equal(time_amplification(cn_data = test_data_cn_6_2,
                                  multiplicity_data = test_data_mult_6_2_GGW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GGW")
})

# 9+2
test_that("AmplificationTimeR runs and produces 9+2 highest output and order", {
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_WGGGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"9+2")
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_WGGGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"WGGGGGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_GWGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"9+2")
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_GWGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GWGGGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_GGWGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"9+2")
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_GGWGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GGWGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_GGGWG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"9+2")
  expect_equal(time_amplification(cn_data = test_data_cn_9_2,
                                  multiplicity_data = test_data_mult_9_2_GGGWG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GGGWG")
})

# 10+2
test_that("AmplificationTimeR runs and produces 10+2 highest output and order", {
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_WGGGGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"10+2")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_WGGGGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"WGGGGGGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GWGGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"10+2")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GWGGGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GWGGGGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GGWGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"10+2")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GGWGGGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GGWGGGG")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GGGWGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"10+2")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GGGWGG,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GGGWGG")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GGGGW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"highest_copy_number"],"10+2")
  expect_equal(time_amplification(cn_data = test_data_cn_10_2,
                                  multiplicity_data = test_data_mult_10_2_GGGGW,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = TRUE)[,"event_order"],"GGGGW")
})
