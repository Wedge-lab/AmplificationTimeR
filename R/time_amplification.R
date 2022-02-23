# Function for timing amplifications

#' Amplification timing
#'
#' Time individual amplification events for a locus that has been gained multiple times.
#' @param cn_data Data frame containing copy number information.
#' @param multiplicity_data Data frame containing multiplicity information.
#' @param sample_id ID name or label for sample. 
#' @param amplification_chrom Chromosome of the amplified region.
#' @param amplification_start Start position of the amplified region.
#' @param amplification_stop End position of the amplified region.  
#' @param is_WGD logical: TRUE indicates that the sample has been whole genome duplicated.
#' @return A data frame containing approximate timing of each amplification, and the most likely order of events. 
#' @export

time_amplification <- function(cn_data, multiplicity_data, 
                               sample_id,
                               amplification_chrom, 
                               amplification_start, 
                               amplification_stop, 
                               is_WGD){
  
  ########
  # Check input
  
  # Input is right type
  if(class(cn_data)[1] != "data.frame"){
    stop("'cn_data' must be an object of type 'data.frame'")
  }
  if(class(multiplicity_data)[1] != "data.frame"){
    stop("'multiplicity_data' must be an object of class 'data.frame'")
  }
  if(!is.character(sample_id)){
    stop("'sample_id' must be an object of type 'character'")
  }
  if((amplification_chrom %in% c("1","2","3","4","5","6","7","8","9","10",
                                 "11","12","13","14","15","16","17","18","19","20",
                                 "21","22","X","Y")) == FALSE){
    stop("Invalid chromosome supplied")
  }
  if(!is.numeric(amplification_start)){
    stop("'amplification_start' must be an object of type 'numeric'")
  }
  if(!is.numeric(amplification_stop)){
    stop("'amplification_stop' must be an object of type 'numeric'")
  }
  if(!is.logical(is_WGD)){
    stop("'is_WGD' must be either 'TRUE' or 'FALSE'")
  }
  
  # Input has right columns
  if(!all(c("chr","startpos","endpos","nMaj1_A","nMin1_A") %in% colnames(cn_data))){
    stop("Incorrect column names in 'cn_data'")
  }
  if(!all(c("chr","end","no.chrs.bearing.mut") %in% colnames(multiplicity_data))){
    stop("Incorrect column names in 'multiplicity_data'")
  }
  
  # Input has non-zero number of rows
  if(nrow(cn_data) == 0){
    stop("'cn_file' has no content")
  }
  if(nrow(multiplicity_data) == 0){
    stop("'multiplicity_data' has no content")
  }
  
  ##############################################################################
  # subset copy number file for amplified region
  tmp_cn <- subset(cn_data, chr == amplification_chrom &
                     startpos <= amplification_stop &
                     endpos >= amplification_start)
  
  if(nrow(tmp_cn) == 0){
    stop("Cannot subset 'cn_data' for this region.  Double-check the region specified, and that copy number has been called in this region. ")
  }
  
  tmp_cn$length <- (tmp_cn$end - tmp_cn$start)+1
  
  tmp_cn$n1A <- paste(tmp_cn$nMaj1_A, "+", tmp_cn$nMin1_A, sep = "")
  tmp_cn$n1A_sum <- tmp_cn$nMaj1_A + tmp_cn$nMin1_A
  
  if(c("nMaj2_A") %in% colnames(tmp_cn)){
    tmp_cn$n2A <- paste(tmp_cn$nMaj2_A, "+", tmp_cn$nMin2_A, sep = "")
    tmp_cn$n2A_sum <- tmp_cn$nMaj2_A + tmp_cn$nMin2_A
  }
  
  # subset multiplicity for amplified region
  tmp_mult <- subset(multiplicity_data, chr == amplification_chrom & 
                       end >= amplification_start & 
                       end <= amplification_stop )
  if(nrow(tmp_mult) == 0){
    stop("Cannot subset 'multiplicity_data' for this region.  Cannot run analysis if there are no mutations in this region.")
  }
  ###
  # Assuming region is spanned by 1 copy number segment
  # adjust later to account for partially gained regions
  ###
  
  # Check if region is amplified
  # Cannot run code if region is not amplified
  is_amplified <- NA
  
  if(c("n2A_sum") %in% colnames(tmp_cn)){
    if(is_WGD == FALSE & tmp_cn$n1A_sum <= 2 & tmp_cn$n2A_sum <= 2){
      stop("Region not amplified in sample")
      
    }else if(is_WGD == TRUE & tmp_cn$n1A_sum <= 4 & tmp_cn$n2A_sum <= 4){
      stop("Region not amplified in sample")
      
    }else{
      is_amplified <- TRUE
    }
  }else{
    if(is_WGD == FALSE & tmp_cn$n1A_sum <= 2){
      stop("Region not amplified in sample")
      
    }else if(is_WGD == TRUE & tmp_cn$n1A_sum <= 4){
      stop("Region not amplified in sample")
      
    }else{
      is_amplified <- TRUE
    }
  }
  
  # get the highest copy number value - assume lower value has same timing
  if(c("n2A_sum") %in% colnames(tmp_cn)){
    max_amplification <- names(which.max(tmp_cn[,c("n1A_sum","n2A_sum")]))
    if(max_amplification == c("n1A_sum")){
      max_amplification_split <- tmp_cn$n1A
    }else if(max_amplification == c("n2A_sum")){
      max_amplification_split <- tmp_cn$n2A
    }
  }else{
    max_amplification_split <- tmp_cn$n1A
  }
  
  if(is_amplified == TRUE){
    tmp_mult$no.chrs.bearing.mut.ceiling <- ceiling(tmp_mult$no.chrs.bearing.mut)
    ############################################################################
    ### calculate the multiplicity of all of the mutations
    ## First create a table to house your multiplicity values
    
    if(nrow(tmp_mult > 0)){
      n1 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 1))
      n2 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 2))
      n3 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 3))
      n4 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 4))
      n5 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 5))
      n6 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 6))
      n7 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 7))
      n8 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 8))
      n9 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 9))
      n10 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 10))
      n11 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 11))
      n12 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 12))
      n13 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 13))
      n14 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 14))
      n15 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 15))
      n16 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 16))
      n17 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 17))
      n18 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 18))
      n19 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 19))
      n20 <- nrow(subset(tmp_mult, no.chrs.bearing.mut.ceiling == 20))
    }
    
    ############################################################################
    # make your output data frame
    amplification_results <- as.data.frame(matrix(nrow=1, ncol = 14))
    colnames(amplification_results) <- c("sample","region","highest_copy_number","event_order",
                                         "t_1","t_2","t_3","t_4","t_5",
                                         "t_6","t_7","t_8","t_9","t_10")
    ############################################################################
    ### start calculating timing
    ### 
    amplification_results$sample <- sample_id
    amplification_results$region <- paste(amplification_chrom,":",amplification_start,"-",amplification_stop, sep = "")
    amplification_results$highest_copy_number <- max_amplification_split
    
    if(max_amplification_split == c("2+1") & is_WGD == FALSE){ ######## start with non-WGD
      
      t_1 <- (2*n2)/(n1+2)
      amplification_results$event_order <- "G"
      amplification_results$t_1 <- t_1
      
    }else if(max_amplification_split == c("2+2") & is_WGD == FALSE){
      
      amplification_results$event_order <- "Cannot be timed"

    }else if(max_amplification_split == c("3+1") & is_WGD == FALSE){
      
      t_1 <- ((4*n3)/(n1 + 2*n2 + 3*n3))
      t_2 <- ((4*(n2+n3))/(n1 + 2*n2 + 3*n3))
      amplification_results$event_order <- "GG"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      
    }else if(max_amplification_split == c("3+2") & is_WGD == FALSE){
      
      amplification_results$event_order <- "Cannot be timed"
      
    }else if(max_amplification_split == c("4+1") & is_WGD == FALSE){
      
      t_1 <- ((5*(n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
      t_2 <- ((5*(n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
      t_3 <- ((5*(n2 + n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
      amplification_results$event_order <- "GGG"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      amplification_results$t_3 <- t_3
      
    }else if(max_amplification_split == c("4+2") & is_WGD == TRUE){ ####### WGD onwards, with further options
      # WGG
      if(n3 > 0){
        t_1 <- ((6*(n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
        t_2 <- ((6*(n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
        t_3 <- ((6*(n2 + n3))/(n1 + 2*n2 + 3*n3 + 4*n4))
        amplification_results$event_order <- "WGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
      }else if(n3 == 0){# GW
        t_1 <- ((6*(n4))/(n1 + 2*n2 + 4*n4))
        t_2 <- ((2*(n2 + 2*n4))/(n1 + 2*n2 + 4*n4))
        amplification_results$event_order <- "GW"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
      }
      
    }else if(max_amplification_split == c("4+3") & is_WGD == TRUE){
      
    }else if(max_amplification_split == c("4+4") & is_WGD == TRUE){
      
    }else if(max_amplification_split == c("5+1") & is_WGD == TRUE){
      # WGGG
      if(n4 > 0){
        t_1 <- ((6*(n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        t_2 <- ((6*(n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        t_3 <- ((6*(n3 + n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        t_4 <- ((6*(n2 + n3 + n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        amplification_results$event_order <- "WGGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
      }else if(n4 == 0){# GWG
        t_1 <- ((6*(n5))/(n1 + 2*n2 + 3*n3 + 5*n5))
        t_2 <- ((6*(n3 + n5))/(n1 + 2*n2 + 3*n3 + 5*n5))
        t_3 <- ((6*(n2 + n5))/(n1 + 2*n2 + 3*n3 + 5*n5))
        amplification_results$event_order <- "GWG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
      }
      
    }else if(max_amplification_split == c("5+2") & is_WGD == TRUE){
      # WGGG
      if(n4 > 0){
        t_1 <- ((7*(n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        t_2 <- ((7*(n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        t_3 <- ((7*(n3 + n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        t_4 <- ((7*(n2 + n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
        amplification_results$event_order <- "WGGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
      }else if(n4 == 0){# GWG
        t_1 <- ((7*(n5))/(n1 + 2*n2 + 3*n3 + 5*n5))
        t_2 <- ((7*(n3 + n5))/(n1 + 2*n2 + 3*n3 + 5*n5))
        t_3 <- (((n2 - n3))/(n1 + 2*n2 + 3*n3 + 5*n5))
        amplification_results$event_order <- "GWG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
      }
      
    }else if(max_amplification_split == c("6+1") & is_WGD == TRUE){
      # WGGGG
      if((n3 > 0) & (n5 > 0)){
        t_1 <- ((7*(n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_2 <- ((7*(n5 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_3 <- ((7*(n4 + n5 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_4 <- ((7*(n3 + n4 + n5 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_5 <- ((7*(n2 + n3 + n4 + n5 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        amplification_results$event_order <- "WGGGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
        amplification_results$t_5 <- t_5
      }else if((n3 > 0) & (n5 == 0)){# GWGG
        t_1 <- ((7*(n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
        t_2 <- ((7*(n4 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
        t_3 <- ((7*(n3 + n4 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
        t_4 <- ((7*(n2 + n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
        amplification_results$event_order <- "GWGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
      }else if((n3 == 0) & (n5 == 0)){# GGW
        t_1 <- ((7*(n6))/(n1 + 2*n2 + 6*n4 + 6*n6))
        t_2 <- ((7*(n4 + n6))/(n1 + 2*n2 + 6*n4 + 6*n6))
        t_3 <- ((7*(n2 + 3*n4 + 3*n6))/(n1 + 2*n2 + 6*n4 + 6*n6))
        amplification_results$event_order <- "GGW"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
      }
      
    }else if(max_amplification_split == c("6+2") & is_WGD == TRUE){
      # WGGGG
      if((n3 > 0) & (n5 > 0)){
        t_1 <- ((8*(n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_2 <- ((8*(n5 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_3 <- ((8*(n4 + n5 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_4 <- ((8*(n3 + n4 + n5 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        t_5 <- ((8*(n2 + n3 + n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6))
        amplification_results$event_order <- "WGGGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
        amplification_results$t_5 <- t_5
      }else if((n3 > 0) & (n5 == 0)){# GWGG
        t_1 <- ((8*(n6))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
        t_2 <- ((8*(n4 + n6))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
        t_3 <- ((8*(n3 + n4 + n6))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
        t_4 <- ((8*(n2 + n3 - n4))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
        amplification_results$event_order <- "GWGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
      }else if((n3 == 0) & (n5 == 0)){# GGW
        t_1 <- ((8*(n6))/(n1 + 2*n2 + 4*n4 + 6*n6))
        t_2 <- ((8*(n4 + n6))/(n1 + 2*n2 + 4*n4 + 6*n6))
        t_3 <- ((2*(n2 + 2*n4 + 3*n6))/(n1 + 2*n2 + 4*n4 + 6*n6))
        amplification_results$event_order <- "GGW"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
      }
    }else if(max_amplification_split == c("10+2") & is_WGD == TRUE){
      # WGGGGGGGG  
      if((n3 > 0) & (n5 > 0) & (n7 > 0) & (n9 > 0)){
        t_1 <- ((12*(n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_2 <- ((12*(n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_3 <- ((12*(n8 + n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_4 <- ((12*(n7 + n8 + n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_5 <- ((12*(n6 + n7 + n8 + n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_6 <- ((12*(n5 + n6 + n7 + n8 + n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_7 <- ((12*(n4 + n5 + n6 + n7 + n8 + n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_8 <- ((12*(n3 + n4 + n5 + n6 + n7 + n8 + n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        t_9 <- ((12*(n2 + n3 + n4 + n5 + n6 + n7 + n8 + n9 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 9*n9 + 10*n10))
        amplification_results$event_order <- "WGGGGGGGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
        amplification_results$t_5 <- t_5
        amplification_results$t_6 <- t_6
        amplification_results$t_7 <- t_7
        amplification_results$t_8 <- t_8
        amplification_results$t_9 <- t_9
      }else if((n3 > 0) & (n5 > 0) & (n7 > 0) & (n9 == 0)){# GWGGGGGG 
        t_1 <- ((12*(n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        t_2 <- ((12*(n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        t_3 <- ((12*(n7 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        t_4 <- ((12*(n6 + n7 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        t_5 <- ((12*(n5 + n6 + n7 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        t_6 <- ((12*(n4 + n5 + n6 + n7 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        t_7 <- ((12*(n3 + n4 + n5 + n6 + n7 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        t_8 <- ((12*(n2 + n3 + n4 + n5 + n6 + n7 - n8))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 7*n7 + 8*n8 + 10*n10))
        amplification_results$event_order <- "GWGGGGGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
        amplification_results$t_5 <- t_5
        amplification_results$t_6 <- t_6
        amplification_results$t_7 <- t_7
        amplification_results$t_8 <- t_8
      }else if((n3 > 0) & (n5 > 0) & (n7 == 0) & (n9 == 0)){# GGWGGGG
        t_1 <- ((12*(n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 10*n8 + 10*n10))
        t_2 <- ((12*(n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 10*n8 + 10*n10))
        t_3 <- ((12*(n6 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 10*n8 + 10*n10))
        t_4 <- ((12*(n5 + n6 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 10*n8 + 10*n10))
        t_5 <- ((12*(n4 + n5 + n6 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 10*n8 + 10*n10))
        t_6 <- ((12*(n3 + n4 + n5 + n6 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 10*n8 + 10*n10))
        t_7 <- ((12*(n2 + n3 + n4 + n5 - 2*n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5 + 6*n6 + 10*n8 + 10*n10))
        amplification_results$event_order <- "GGWGGGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
        amplification_results$t_5 <- t_5
        amplification_results$t_6 <- t_6
        amplification_results$t_7 <- t_7
      }else if((n3 > 0) & (n5 == 0) & (n7 == 0) & (n9 == 0)){# GGGWGG
        t_1 <- ((12*(n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 10*n6 + 10*n8 + 10*n10))
        t_2 <- ((12*(n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 10*n6 + 10*n8 + 10*n10))
        t_3 <- ((12*(n6 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 10*n6 + 10*n8 + 10*n10))
        t_4 <- ((12*(n4 + n6 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 10*n6 + 10*n8 + 10*n10))
        t_5 <- ((12*(n3 + n4 + n6 + n8 + n10))/(n1 + 2*n2 + 3*n3 + 4*n4 + 10*n6 + 10*n8 + 10*n10))
        t_6 <- ((12*(n2 + n3 - 3*n4))/(n1 + 2*n2 + 3*n3 + 4*n4 + 10*n6 + 10*n8 + 10*n10))
        amplification_results$event_order <- "GGGWGG"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
        amplification_results$t_5 <- t_5
        amplification_results$t_6 <- t_6
      }else if((n3 == 0) & (n5 == 0) & (n7 == 0) & (n9 == 0)){# GGGGW
        t_1 <- ((12*(n10))/(n1 + 2*n2 + 10*n4 + 10*n6 + 10*n8 + 10*n10))
        t_2 <- ((12*(n8 + n10))/(n1 + 2*n2 + 10*n4 + 10*n6 + 10*n8 + 10*n10))
        t_3 <- ((12*(n6 + n8 + n10))/(n1 + 2*n2 + 10*n4 + 10*n6 + 10*n8 + 10*n10))
        t_4 <- ((12*(n4 + n6 + n8 + n10))/(n1 + 2*n2 + 10*n4 + 10*n6 + 10*n8 + 10*n10))
        t_5 <- ((2*(n2 + 5*n4 + 5*n6 + 5*n8 + 5*n10))/(n1 + 2*n2 + 10*n4 + 10*n6 + 10*n8 + 10*n10))
        amplification_results$event_order <- "GGGGW"
        amplification_results$t_1 <- t_1
        amplification_results$t_2 <- t_2
        amplification_results$t_3 <- t_3
        amplification_results$t_4 <- t_4
        amplification_results$t_5 <- t_5
      }else{
        amplification_results$event_order <- "Something went wrong"
      }
    }
  }
  return(amplification_results)
}

