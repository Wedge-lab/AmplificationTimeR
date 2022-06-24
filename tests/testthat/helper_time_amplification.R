# creating test data for testing AmplificationTimeR

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

test_output <- as.data.frame(matrix(nrow = 1, ncol = 34, data = c("test_data","1:1-20","2+1","G","0.63958596","0.4899845","0.73335304","NA","NA","NA","NA","NA","NA","NA",
                                                                  "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA")))
colnames(test_output) <- c("sample","region","highest_copy_number","event_order","t_1","t_1_lower_ci","t_1_upper_ci",
                           "t_2","t_2_lower_ci","t_2_upper_ci","t_3","t_3_lower_ci","t_3_upper_ci",
                           "t_4","t_4_lower_ci","t_4_upper_ci","t_5","t_5_lower_ci","t_5_upper_ci",
                           "t_6","t_6_lower_ci","t_6_upper_ci","t_7","t_7_lower_ci","t_7_upper_ci",
                           "t_8","t_8_lower_ci","t_8_upper_ci","t_9","t_9_lower_ci","t_9_upper_ci",
                           "t_10","t_10_lower_ci","t_10_upper_ci")
test_output$t_1 <- as.numeric(test_output$t_1)
test_output$t_1_lower_ci <- as.numeric(test_output$t_1_lower_ci)
test_output$t_1_upper_ci <- as.numeric(test_output$t_1_upper_ci)
test_output[,8:34]<- as.logical(test_output[,8:34])

test_data_cn_incorrect <- test_data_cn
colnames(test_data_cn_incorrect) <- c("chr","startpos","endpos","nMaj1_A","x")

test_data_cn_empty <- test_data_cn[-1,]

test_data_mult_incorrect <- test_data_mult
colnames(test_data_mult_incorrect) <- c("chr","end","x")

test_data_mult_empty <- test_data_mult[-c(1:20),]

test_data_cn_n2 <- test_data_cn
test_data_cn_n2$nMaj2_A <- 2
test_data_cn_n2$nMin2_A <- 0

# test that code fails if region supplied is not amplified
test_data_cn_not_amp <- test_data_cn
test_data_cn_not_amp$nMaj1_A <- 1
test_data_cn_not_amp$nMin1_A <- 1

# test multiplicity has mutations in region
test_data_mult_no_muts_in_region <- test_data_mult
test_data_mult_no_muts_in_region$end <- test_data_mult_no_muts_in_region$end+40

# test that code fails if region supplied is not amplified with n2A present in diploid sample
test_data_cn_not_amp_n2A <- test_data_cn_not_amp
test_data_cn_not_amp$nMaj2_A <- 2
test_data_cn_not_amp$nMin2_A <- 0

# test that code fails if region supplied is not amplified - WGD
test_data_cn_not_amp_wgd <- test_data_cn
test_data_cn_not_amp_wgd$nMaj1_A <- 2
test_data_cn_not_amp_wgd$nMin1_A <- 2

# test that code fails if region supplied is not amplified with n2A present in WGD sample
test_data_cn_not_amp_wgd_n2A <- test_data_cn_not_amp_wgd
test_data_cn_not_amp_wgd_n2A$nMaj2_A <- 2
test_data_cn_not_amp_wgd_n2A$nMin2_A <- 0


# test that n2A_sum will be chosen if it is higher than n1A_sum
test_data_cn_n2_higher <- test_data_cn
test_data_cn_n2_higher$nMaj1_A <- 1
test_data_cn_n2_higher$nMin1_A <- 1
test_data_cn_n2_higher$nMaj2_A <- 2
test_data_cn_n2_higher$nMin2_A <- 1

########
# create test data sets for different copy numbers and test output
###
# 2+2
test_data_cn_2_2 <- test_data_cn
test_data_cn_2_2$nMin1_A <- 2
test_data_cn_2_2$nMaj1_A <- 2

# 3+1
test_data_cn_3_1 <- test_data_cn
test_data_cn_3_1$nMaj1_A <- 3
test_data_cn_3_1$nMin1_A <- 1

# 3+2
test_data_cn_3_2 <- test_data_cn
test_data_cn_3_2$nMaj1_A <- 3
test_data_cn_3_2$nMin1_A <- 2

# 4+1
test_data_cn_4_1 <- test_data_cn
test_data_cn_4_1$nMaj1_A <- 4
test_data_cn_4_1$nMin1_A <- 1

# 4+2
test_data_cn_4_2 <- test_data_cn
test_data_cn_4_2$nMaj1_A <- 4
test_data_cn_4_2$nMin1_A <- 2
test_data_mult_4_2_WGG <- test_data_mult
test_data_mult_4_2_WGG$no.chrs.bearing.mut[1:2] <- 3
test_data_mult_4_2_GW <- test_data_mult

# 4+3
test_data_cn_4_3 <- test_data_cn
test_data_cn_4_3$nMaj1_A <- 4
test_data_cn_4_3$nMin1_A <- 3

# 4+4
test_data_cn_4_4 <- test_data_cn
test_data_cn_4_4$nMaj1_A <- 4
test_data_cn_4_4$nMin1_A <- 4

# 5+1
test_data_cn_5_1 <- test_data_cn
test_data_cn_5_1$nMaj1_A <- 5
test_data_cn_5_1$nMin1_A <- 1
test_data_mult_5_1_WGGG <- test_data_mult
test_data_mult_5_1_WGGG$no.chrs.bearing.mut[1:2] <- 4
test_data_mult_5_1_GWG <- test_data_mult

# 5+2
test_data_cn_5_2 <- test_data_cn
test_data_cn_5_2$nMaj1_A <- 5
test_data_cn_5_2$nMin1_A <- 2
test_data_mult_5_2_WGGG <- test_data_mult
test_data_mult_5_2_WGGG$no.chrs.bearing.mut[1:2] <- 4
test_data_mult_5_2_GWG <- test_data_mult

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

