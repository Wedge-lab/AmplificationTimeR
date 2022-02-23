# Test time_amplification.R
library(testthat)
library(covr)

# set up test data that should work (2+1 not WGD)
test_data_cn <- as.data.frame(matrix(nrow = 1, ncol = 5, data = c(1,1,20,2,1)))
colnames(test_data_cn) <- c("chr","startpos","endpos","nMaj1_A","nMin1_A")
test_data_cn$startpos <- as.numeric(test_data_cn$startpos)
test_data_cn$endpos <- as.numeric(test_data_cn$endpos)
test_data_cn$nMaj1_A <- as.numeric(test_data_cn$nMaj1_A)
test_data_cn$nMin1_A <- as.numeric(test_data_cn$nMin1_A)
test_data_mult <- as.data.frame(matrix(nrow = 20, ncol = 3, data = c(rep("1",20),1:20,rep(2,5),rep(1,15))))
colnames(test_data_mult) <- c("chr","end","no.chrs.bearing.mut")
test_data_mult$end <- as.numeric(test_data_mult$end)
test_data_mult$no.chrs.bearing.mut <- as.numeric(test_data_mult$no.chrs.bearing.mut)

test_data_id <- "test_data"
test_data_chrom <- "1"
test_data_start <- 1
test_data_stop <- 20
test_data_status <- FALSE

test_output <- as.data.frame(matrix(nrow = 1, ncol = 14, data = c("test_data","1:1-20","2+1","G","0.5882353","NA","NA","NA","NA","NA","NA","NA","NA","NA")))
colnames(test_output) <- c("sample","region","highest_copy_number","event_order","t_1","t_2","t_3","t_4","t_5","t_6","t_7","t_8","t_9","t_10")
test_output$t_1 <- as.numeric(test_output$t_1)
test_output[,6:14]<- as.logical(test_output[,6:14])

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
test_data_cn_incorrect <- test_data_cn
colnames(test_data_cn_incorrect) <- c("chr","startpos","endpos","nMaj1_A","x")
test_that("AmplificationTimeR fails and produces an error message when copy number input has incorrect columns.", {
  expect_error(time_amplification(cn_data = test_data_cn_incorrect,
                                  multiplicity_data = test_data_mult,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
test_data_cn_empty <- test_data_cn[-1,]
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
test_data_mult_incorrect <- test_data_mult
colnames(test_data_mult_incorrect) <- c("chr","end","x")
test_that("AmplificationTimeR fails and produces an error message when multiplicity input has incorrect columns.", {
  expect_error(time_amplification(cn_data = test_data_cn,
                                  multiplicity_data = test_data_mult_incorrect,
                                  sample_id = test_data_id,
                                  amplification_chrom = test_data_chrom,
                                  amplification_start = test_data_start,
                                  amplification_stop = test_data_stop,
                                  is_WGD = test_data_status))
})
test_data_mult_empty <- test_data_mult[-c(1:20),]
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
test_data_cn_n2 <- test_data_cn
test_data_cn_n2$nMaj2_A <- 2
test_data_cn_n2$nMin2_A <- 0
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
test_data_cn_not_amp <- test_data_cn
test_data_cn_not_amp$nMaj1_A <- 1
test_data_cn_not_amp$nMin1_A <- 1
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
test_data_mult_no_muts_in_region <- test_data_mult
test_data_mult_no_muts_in_region$end <- test_data_mult_no_muts_in_region$end+40
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
test_data_cn_not_amp_n2A <- test_data_cn_not_amp
test_data_cn_not_amp$nMaj2_A <- 2
test_data_cn_not_amp$nMin2_A <- 0
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
test_data_cn_not_amp_wgd <- test_data_cn
test_data_cn_not_amp_wgd$nMaj1_A <- 2
test_data_cn_not_amp_wgd$nMin1_A <- 2
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
test_data_cn_not_amp_wgd_n2A <- test_data_cn_not_amp_wgd
test_data_cn_not_amp_wgd_n2A$nMaj2_A <- 2
test_data_cn_not_amp_wgd_n2A$nMin2_A <- 0
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
test_data_cn_n2_higher <- test_data_cn
test_data_cn_n2_higher$nMaj1_A <- 1
test_data_cn_n2_higher$nMin1_A <- 1
test_data_cn_n2_higher$nMaj2_A <- 2
test_data_cn_n2_higher$nMin2_A <- 1
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
test_that("AmplificationTimeR runs and produces a data.frame with columns 'sample','region','highest_copy_number','event_order',
                                         't_1','t_2','t_3','t_4','t_5',
                                         't_6','t_7','t_8','t_9','t_10' when all input is correct", {
                                           expect_equal(colnames(time_amplification(cn_data = test_data_cn,
                                                                                    multiplicity_data = test_data_mult,
                                                                                    sample_id = test_data_id,
                                                                                    amplification_chrom = test_data_chrom,
                                                                                    amplification_start = test_data_start,
                                                                                    amplification_stop = test_data_stop,
                                                                                    is_WGD = test_data_status)),c("sample","region","highest_copy_number","event_order",
                                                                                                                  "t_1","t_2","t_3","t_4","t_5",
                                                                                                                  "t_6","t_7","t_8","t_9","t_10"))
                                         })

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
# create test data sets for different copy numbers and test output
###
# 2+2
test_data_cn_2_2 <- test_data_cn
test_data_cn_2_2$nMin1_A <- 2
test_data_cn_2_2$nMaj1_A <- 2
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
test_data_cn_3_1 <- test_data_cn
test_data_cn_3_1$nMaj1_A <- 3
test_data_cn_3_1$nMin1_A <- 1
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
test_data_cn_3_2 <- test_data_cn
test_data_cn_3_2$nMaj1_A <- 3
test_data_cn_3_2$nMin1_A <- 2
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
test_data_cn_4_1 <- test_data_cn
test_data_cn_4_1$nMaj1_A <- 4
test_data_cn_4_1$nMin1_A <- 1
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
test_data_cn_4_2 <- test_data_cn
test_data_cn_4_2$nMaj1_A <- 4
test_data_cn_4_2$nMin1_A <- 2
test_data_mult_4_2_WGG <- test_data_mult
test_data_mult_4_2_WGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_4_2_GW <- test_data_mult
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
test_data_cn_4_3 <- test_data_cn
test_data_cn_4_3$nMaj1_A <- 4
test_data_cn_4_3$nMin1_A <- 3
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
test_data_cn_4_4 <- test_data_cn
test_data_cn_4_4$nMaj1_A <- 4
test_data_cn_4_4$nMin1_A <- 4
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
test_data_cn_5_1 <- test_data_cn
test_data_cn_5_1$nMaj1_A <- 5
test_data_cn_5_1$nMin1_A <- 1
test_data_mult_5_1_WGGG <- test_data_mult
test_data_mult_5_1_WGGG$no.chrs.bearing.mut[1:2] <- 4
test_data_mult_5_1_GWG <- test_data_mult
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
test_data_cn_5_2 <- test_data_cn
test_data_cn_5_2$nMaj1_A <- 5
test_data_cn_5_2$nMin1_A <- 2
test_data_mult_5_2_WGGG <- test_data_mult
test_data_mult_5_2_WGGG$no.chrs.bearing.mut[1:2] <- 4
test_data_mult_5_2_GWG <- test_data_mult
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
test_data_cn_6_1 <- test_data_cn
test_data_cn_6_1$nMaj1_A <- 6
test_data_cn_6_1$nMin1_A <- 1
test_data_mult_6_1_WGGGG <- test_data_mult
test_data_mult_6_1_WGGGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_6_1_WGGGG$no.chrs.bearing.mut[3:4] <- 5
test_data_mult_6_1_GWGG <- test_data_mult
test_data_mult_6_1_GWGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_6_1_GGW <- test_data_mult

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
test_data_cn_6_2 <- test_data_cn
test_data_cn_6_2$nMaj1_A <- 6
test_data_cn_6_2$nMin1_A <- 2
test_data_mult_6_2_WGGGG <- test_data_mult
test_data_mult_6_2_WGGGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_6_2_WGGGG$no.chrs.bearing.mut[3:4] <- 5
test_data_mult_6_2_GWGG <- test_data_mult
test_data_mult_6_2_GWGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_6_2_GGW <- test_data_mult
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

# 10+2
test_data_cn_10_2 <- test_data_cn
test_data_cn_10_2$nMaj1_A <- 10
test_data_cn_10_2$nMin1_A <- 2
test_data_mult_10_2_WGGGGGGGG <- test_data_mult
test_data_mult_10_2_WGGGGGGGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_10_2_WGGGGGGGG$no.chrs.bearing.mut[3:4] <- 5
test_data_mult_10_2_WGGGGGGGG$no.chrs.bearing.mut[5:6] <- 7
test_data_mult_10_2_WGGGGGGGG$no.chrs.bearing.mut[7:8] <- 9
test_data_mult_10_2_GWGGGGGG <- test_data_mult
test_data_mult_10_2_GWGGGGGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_10_2_GWGGGGGG$no.chrs.bearing.mut[3:4] <- 5
test_data_mult_10_2_GWGGGGGG$no.chrs.bearing.mut[5:6] <- 7
test_data_mult_10_2_GGWGGGG <- test_data_mult
test_data_mult_10_2_GGWGGGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_10_2_GGWGGGG$no.chrs.bearing.mut[3:4] <- 5
test_data_mult_10_2_GGGWGG <- test_data_mult
test_data_mult_10_2_GGGWGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_10_2_GGGGW <- test_data_mult

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
