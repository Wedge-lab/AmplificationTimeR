# Function for timing amplifications

#' Calculate maths for timing
#' 
#' Hidden function for calculating the maths used in AmplificationTimeR - not intended for use by end users.  Called by time_amplification
#' @param mult_data List of observed multiplicity values
#' @param max_amp Maximum copy number state, split into major and minor copy number
#' @param is_WGD Whole genome duplication status
#' @param ordering_event Order of events determined based on copy number and multiplicity of mutations
#' @return A data frame containing approximate timing of each amplification, and the most likely order of events. 
#' @keywords internal
#' @export

time_amplification_maths <- function(mult_data, max_amp, is_WGD, ordering_event){
  max_amplification_split <- max_amp
  order_event <- ordering_event
  ############################################################################
  # make your output data frame
  amplification_results <- as.data.frame(matrix(nrow=1, ncol = 12))
  colnames(amplification_results) <- c("highest_copy_number","event_order",
                                       "t_1","t_2","t_3","t_4","t_5",
                                       "t_6","t_7","t_8","t_9","t_10")
  amplification_results$highest_copy_number <- max_amplification_split
  ############################################################################
  ### calculate the multiplicity of all of the mutations
  ## First create a table to house your multiplicity values
  
  if(length(mult_data > 0)){
    n1 <- sum(mult_data == 1) 
    n2 <- sum(mult_data == 2) 
    n3 <- sum(mult_data == 3) 
    n4 <- sum(mult_data == 4) 
    n5 <- sum(mult_data == 5) 
    n6 <- sum(mult_data == 6) 
    n7 <- sum(mult_data == 7) 
    n8 <- sum(mult_data == 8) 
    n9 <- sum(mult_data == 9) 
    n10 <- sum(mult_data == 10) 
    n11 <- sum(mult_data == 11) 
    n12 <- sum(mult_data == 12) 
    n13 <- sum(mult_data == 13) 
    n14 <- sum(mult_data == 14) 
    n15 <- sum(mult_data == 15) 
    n16 <- sum(mult_data == 16) 
    n17 <- sum(mult_data == 17) 
    n18 <- sum(mult_data == 18) 
    n19 <- sum(mult_data == 19) 
    n20 <- sum(mult_data == 20) 
  }
  
  
  ############################################################################
  ### start calculating timing
  ### 
  
  if(max_amplification_split == c("2+1") & is_WGD == FALSE & order_event == "G"){ ######## start with non-WGD
    
    t_1 <- (2*n2)/(n1+2)
    amplification_results$event_order <- "G"
    amplification_results$t_1 <- t_1
    
  }else if(max_amplification_split == c("2+2") & is_WGD == FALSE & order_event == "Cannot be timed"){
    
    amplification_results$event_order <- "Cannot be timed"
    
  }else if(max_amplification_split == c("3+1") & is_WGD == FALSE & order_event == "GG"){
    
    t_1 <- ((4*n3)/(n1 + 2*n2 + 3*n3))
    t_2 <- ((4*(n2+n3))/(n1 + 2*n2 + 3*n3))
    amplification_results$event_order <- "GG"
    amplification_results$t_1 <- t_1
    amplification_results$t_2 <- t_2
    
  }else if(max_amplification_split == c("3+2") & is_WGD == FALSE & order_event == "Cannot be timed"){
    
    amplification_results$event_order <- "Cannot be timed"
    
  }else if(max_amplification_split == c("4+1") & is_WGD == FALSE & order_event == "GGG"){
    
    t_1 <- ((5*(n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
    t_2 <- ((5*(n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
    t_3 <- ((5*(n2 + n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
    amplification_results$event_order <- "GGG"
    amplification_results$t_1 <- t_1
    amplification_results$t_2 <- t_2
    amplification_results$t_3 <- t_3
    
  }else if(max_amplification_split == c("4+2") & is_WGD == TRUE){ ####### WGD onwards, with further options
    # WGG
    if(n3 > 0 & order_event == "WGG"){
      t_1 <- ((6*(n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
      t_2 <- ((6*(n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4))
      t_3 <- ((6*(n2 + n3))/(n1 + 2*n2 + 3*n3 + 4*n4))
      amplification_results$event_order <- "WGG"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      amplification_results$t_3 <- t_3
    }else if(n3 == 0 & order_event == "GW"){# GW
      t_1 <- ((6*(n4))/(n1 + 2*n2 + 4*n4))
      t_2 <- ((2*(n2 + 2*n4))/(n1 + 2*n2 + 4*n4))
      amplification_results$event_order <- "GW"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
    }
    
  }else if(max_amplification_split == c("4+3") & is_WGD == TRUE & order_event == "Cannot be timed"){
    
  }else if(max_amplification_split == c("4+4") & is_WGD == TRUE & order_event == "Cannot be timed"){
    
  }else if(max_amplification_split == c("5+1") & is_WGD == TRUE){
    # WGGG
    if(n4 > 0 & order_event == "WGGG"){
      t_1 <- ((6*(n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      t_2 <- ((6*(n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      t_3 <- ((6*(n3 + n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      t_4 <- ((6*(n2 + n3 + n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      amplification_results$event_order <- "WGGG"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      amplification_results$t_3 <- t_3
      amplification_results$t_4 <- t_4
    }else if(n4 == 0 & order_event == "GWG"){# GWG
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
    if(n4 > 0 & order_event == "WGGG"){
      t_1 <- ((7*(n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      t_2 <- ((7*(n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      t_3 <- ((7*(n3 + n4 + n5))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      t_4 <- ((7*(n2 + n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4 + 5*n5))
      amplification_results$event_order <- "WGGG"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      amplification_results$t_3 <- t_3
      amplification_results$t_4 <- t_4
    }else if(n4 == 0 & order_event == "GWG"){# GWG
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
    if((n3 > 0) & (n5 > 0) & order_event == "WGGGG"){
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
    }else if((n3 > 0) & (n5 == 0) & order_event == "GWGG"){# GWGG
      t_1 <- ((7*(n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
      t_2 <- ((7*(n4 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
      t_3 <- ((7*(n3 + n4 + n6))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
      t_4 <- ((7*(n2 + n3 + n4))/(n1 + 2*n2 + 3*n3 + 4*n4 + 6*n6))
      amplification_results$event_order <- "GWGG"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      amplification_results$t_3 <- t_3
      amplification_results$t_4 <- t_4
    }else if((n3 == 0) & (n5 == 0) & order_event == "GGW"){# GGW
      t_1 <- ((7*(n6))/(n1 + 2*n2 + 4*n4 + 6*n6))
      t_2 <- ((7*(n4 + n6))/(n1 + 2*n2 + 4*n4 + 6*n6))
      t_3 <- ((7*(n2 + 2*n4 + 3*n6))/(n1 + 2*n2 + 4*n4 + 6*n6))
      amplification_results$event_order <- "GGW"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      amplification_results$t_3 <- t_3
    }
    
  }else if(max_amplification_split == c("6+2") & is_WGD == TRUE){
    # WGGGG
    if((n3 > 0) & (n5 > 0) & order_event == "WGGGG"){
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
    }else if((n3 > 0) & (n5 == 0) & order_event == "GWGG"){# GWGG
      t_1 <- ((8*(n6))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
      t_2 <- ((8*(n4 + n6))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
      t_3 <- ((8*(n3 + n4 + n6))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
      t_4 <- ((8*(n2 + n3 - n4))/(n1 + 2*n2 - n3 + 4*n4 + 6*n6))
      amplification_results$event_order <- "GWGG"
      amplification_results$t_1 <- t_1
      amplification_results$t_2 <- t_2
      amplification_results$t_3 <- t_3
      amplification_results$t_4 <- t_4
    }else if((n3 == 0) & (n5 == 0) & order_event == "GGW"){# GGW
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
    if((n3 > 0) & (n5 > 0) & (n7 > 0) & (n9 > 0) & order_event == "WGGGGGGGG"){
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
    }else if((n3 > 0) & (n5 > 0) & (n7 > 0) & (n9 == 0) & order_event == "GWGGGGGG"){# GWGGGGGG 
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
    }else if((n3 > 0) & (n5 > 0) & (n7 == 0) & (n9 == 0) & order_event == "GGWGGGG"){# GGWGGGG
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
    }else if((n3 > 0) & (n5 == 0) & (n7 == 0) & (n9 == 0) & order_event == "GGGWGG"){# GGGWGG
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
    }else if((n3 == 0) & (n5 == 0) & (n7 == 0) & (n9 == 0) & order_event == "GGGGW"){# GGGGW
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
  return(amplification_results)
}

#' Amplification timing
#'
#' Wrapper function that times individual amplification events for a locus that has been gained multiple times.
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
  tmp_mult$no.chrs.bearing.mut.ceiling <- ceiling(tmp_mult$no.chrs.bearing.mut)
  
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
  
  
  
  # Get all of the multiplicities to feed into time_amplification_maths
  tmp_values <- tmp_mult$no.chrs.bearing.mut.ceiling
  
  if(is_amplified == TRUE){
    amplification_results_ci <- as.data.frame(matrix(nrow=1, ncol = 34))
    colnames(amplification_results_ci) <- c("sample","region","highest_copy_number","event_order",
                                            "t_1","t_1_lower_ci","t_1_upper_ci",
                                            "t_2","t_2_lower_ci","t_2_upper_ci",
                                            "t_3","t_3_lower_ci","t_3_upper_ci",
                                            "t_4","t_4_lower_ci","t_4_upper_ci",
                                            "t_5","t_5_lower_ci","t_5_upper_ci",
                                            "t_6","t_6_lower_ci","t_6_upper_ci",
                                            "t_7","t_7_lower_ci","t_7_upper_ci",
                                            "t_8","t_8_lower_ci","t_8_upper_ci",
                                            "t_9","t_9_lower_ci","t_9_upper_ci",
                                            "t_10","t_10_lower_ci","t_10_upper_ci")
    
    amplification_results_ci$sample <- sample_id
    amplification_results_ci$region <- paste(amplification_chrom,":",amplification_start,"-",amplification_stop, sep = "")
    amplification_results_ci$highest_copy_number <- max_amplification_split
    
    ##############################################################################
    # Get event order
    if(length(tmp_values > 0)){
      n1 <- sum(tmp_values == 1) 
      n2 <- sum(tmp_values == 2) 
      n3 <- sum(tmp_values == 3) 
      n4 <- sum(tmp_values == 4) 
      n5 <- sum(tmp_values == 5) 
      n6 <- sum(tmp_values == 6) 
      n7 <- sum(tmp_values == 7) 
      n8 <- sum(tmp_values == 8) 
      n9 <- sum(tmp_values == 9) 
      n10 <- sum(tmp_values == 10) 
      n11 <- sum(tmp_values == 11) 
      n12 <- sum(tmp_values == 12) 
      n13 <- sum(tmp_values == 13) 
      n14 <- sum(tmp_values == 14) 
      n15 <- sum(tmp_values == 15) 
      n16 <- sum(tmp_values == 16) 
      n17 <- sum(tmp_values == 17) 
      n18 <- sum(tmp_values == 18) 
      n19 <- sum(tmp_values == 19) 
      n20 <- sum(tmp_values == 20) 
    }
    event_ordering <- NA
    if(max_amplification_split == c("2+1") & is_WGD == FALSE){ ######## start with non-WGD
      
      event_ordering <- "G"
      
    }else if(max_amplification_split == c("2+2") & is_WGD == FALSE){
      
      event_ordering <- "Cannot be timed"
      
    }else if(max_amplification_split == c("3+1") & is_WGD == FALSE){
      
      event_ordering <- "GG"
      
    }else if(max_amplification_split == c("3+2") & is_WGD == FALSE){
      
      event_ordering <- "Cannot be timed"
      
    }else if(max_amplification_split == c("4+1") & is_WGD == FALSE){
      
      event_ordering <- "GGG"
      
    }else if(max_amplification_split == c("4+2") & is_WGD == TRUE){ ####### WGD onwards, with further options
      # WGG
      if(n3 > 0){
        event_ordering <- "WGG"
      }else if(n3 == 0){# GW
        event_ordering <- "GW"
      }
      
    }else if(max_amplification_split == c("4+3") & is_WGD == TRUE){
      
    }else if(max_amplification_split == c("4+4") & is_WGD == TRUE){
      
    }else if(max_amplification_split == c("5+1") & is_WGD == TRUE){
      # WGGG
      if(n4 > 0){
        event_ordering <- "WGGG"
      }else if(n4 == 0){# GWG
        event_ordering <- "GWG"
      }
      
    }else if(max_amplification_split == c("5+2") & is_WGD == TRUE){
      # WGGG
      if(n4 > 0){
        event_ordering <- "WGGG"
      }else if(n4 == 0){# GWG
        event_ordering <- "GWG"
      }
      
    }else if(max_amplification_split == c("6+1") & is_WGD == TRUE){
      # WGGGG
      if((n3 > 0) & (n5 > 0)){
        event_ordering <- "WGGGG"
      }else if((n3 > 0) & (n5 == 0)){# GWGG
        event_ordering <- "GWGG"
      }else if((n3 == 0) & (n5 == 0)){# GGW
        event_ordering <- "GGW"
      }
      
    }else if(max_amplification_split == c("6+2") & is_WGD == TRUE){
      # WGGGG
      if((n3 > 0) & (n5 > 0)){
        event_ordering <- "WGGGG"
      }else if((n3 > 0) & (n5 == 0)){# GWGG
        event_ordering <- "GWGG"
      }else if((n3 == 0) & (n5 == 0)){# GGW
        event_ordering <- "GGW"
      }
    }else if(max_amplification_split == c("10+2") & is_WGD == TRUE){
      # WGGGGGGGG  
      if((n3 > 0) & (n5 > 0) & (n7 > 0) & (n9 > 0)){
        event_ordering <- "WGGGGGGGG"
      }else if((n3 > 0) & (n5 > 0) & (n7 > 0) & (n9 == 0)){# GWGGGGGG 
        event_ordering <- "GWGGGGGG"
      }else if((n3 > 0) & (n5 > 0) & (n7 == 0) & (n9 == 0)){# GGWGGGG
        event_ordering <- "GGWGGGG"
      }else if((n3 > 0) & (n5 == 0) & (n7 == 0) & (n9 == 0)){# GGGWGG
        event_ordering <- "GGGWGG"
      }else if((n3 == 0) & (n5 == 0) & (n7 == 0) & (n9 == 0)){# GGGGW
        event_ordering <- "GGGGW"
      }else{
        event_ordering <- "Something went wrong"
      }
    }
    ##############################################################################
    
    single_time <- time_amplification_maths(mult_data = tmp_values, 
                                            max_amp = max_amplification_split, 
                                            is_WGD = is_WGD, 
                                            ordering_event = event_ordering)
    amplification_results_ci$event_order <- single_time$event_order                                                                                         
    
    
    bootstrap_amplification <- as.data.frame(matrix(nrow = 500, ncol = 12))
    colnames(bootstrap_amplification) <- c("highest_copy_number","event_order",
                                           "t_1","t_2","t_3","t_4","t_5",
                                           "t_6","t_7","t_8","t_9","t_10")
    
    for(b in 1:nrow(bootstrap_amplification)){
      set.seed(b) # makes bootstrap results reproducible, so multiple runs will produce same values
      temp_multsample <- sample(tmp_values, length(tmp_values), replace = TRUE)
      bootstrap_amplification[b,1:12] <- time_amplification_maths(mult_data = temp_multsample, 
                                                                  max_amp = max_amplification_split,
                                                                  is_WGD = is_WGD,
                                                                  ordering_event = event_ordering)
    }
    
    
    
    if(!is.na(single_time$t_1)){
      amplification_results_ci$t_1 <- mean(bootstrap_amplification$t_1, na.rm = TRUE)
      t_1_lower_ci <- confint(lm(t_1 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_1_lower_ci_adj <- ((length(tmp_values)*t_1_lower_ci)/(5 + length(tmp_values)))
      t_1_upper_ci <- confint(lm(t_1 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_1_upper_ci_adj <- ((5 + (length(tmp_values)*t_1_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_1_lower_ci <- t_1_lower_ci_adj
      amplification_results_ci$t_1_upper_ci <- t_1_upper_ci_adj
    }
    if(!is.na(single_time$t_2)){
      amplification_results_ci$t_2 <- mean(bootstrap_amplification$t_2, na.rm = TRUE)
      t_2_lower_ci <- confint(lm(t_2 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_2_lower_ci_adj <- ((length(tmp_values)*t_2_lower_ci)/(5 + length(tmp_values)))
      t_2_upper_ci <- confint(lm(t_2 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_2_upper_ci_adj <- ((5 + (length(tmp_values)*t_2_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_2_lower_ci <- t_2_lower_ci_adj
      amplification_results_ci$t_2_upper_ci <- t_2_upper_ci_adj
    }
    if(!is.na(single_time$t_3)){
      amplification_results_ci$t_3 <- mean(bootstrap_amplification$t_3, na.rm = TRUE)
      t_3_lower_ci <- confint(lm(t_3 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_3_lower_ci_adj <- ((length(tmp_values)*t_3_lower_ci)/(5 + length(tmp_values)))
      t_3_upper_ci <- confint(lm(t_3 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_3_upper_ci_adj <- ((5 + (length(tmp_values)*t_3_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_3_lower_ci <- t_3_lower_ci_adj
      amplification_results_ci$t_3_upper_ci <- t_3_upper_ci_adj
    }
    if(!is.na(single_time$t_4)){
      amplification_results_ci$t_4 <- mean(bootstrap_amplification$t_4, na.rm = TRUE)
      t_4_lower_ci <- confint(lm(t_4 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_4_lower_ci_adj <- ((length(tmp_values)*t_4_lower_ci)/(5 + length(tmp_values)))
      t_4_upper_ci <- confint(lm(t_4 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_4_upper_ci_adj <- ((5 + (length(tmp_values)*t_4_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_4_lower_ci <- t_4_lower_ci_adj
      amplification_results_ci$t_4_upper_ci <- t_4_upper_ci_adj
    }
    if(!is.na(single_time$t_5)){
      amplification_results_ci$t_5 <- mean(bootstrap_amplification$t_5, na.rm = TRUE)
      t_5_lower_ci <- confint(lm(t_5 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_5_lower_ci_adj <- ((length(tmp_values)*t_5_lower_ci)/(5 + length(tmp_values)))
      t_5_upper_ci <- confint(lm(t_5 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_5_upper_ci_adj <- ((5 + (length(tmp_values)*t_5_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_5_lower_ci <- t_5_lower_ci_adj
      amplification_results_ci$t_5_upper_ci <- t_5_upper_ci_adj
    }
    if(!is.na(single_time$t_6)){
      amplification_results_ci$t_6 <- mean(bootstrap_amplification$t_6, na.rm = TRUE)
      t_6_lower_ci <- confint(lm(t_6 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_6_lower_ci_adj <- ((length(tmp_values)*t_6_lower_ci)/(5 + length(tmp_values)))
      t_6_upper_ci <- confint(lm(t_6 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_6_upper_ci_adj <- ((5 + (length(tmp_values)*t_6_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_6_lower_ci <- t_6_lower_ci_adj
      amplification_results_ci$t_6_upper_ci <- t_6_upper_ci_adj
    }
    if(!is.na(single_time$t_7)){
      amplification_results_ci$t_7 <- mean(bootstrap_amplification$t_7, na.rm = TRUE)
      t_7_lower_ci <- confint(lm(t_7 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_7_lower_ci_adj <- ((length(tmp_values)*t_7_lower_ci)/(5 + length(tmp_values)))
      t_7_upper_ci <- confint(lm(t_7 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_7_upper_ci_adj <- ((5 + (length(tmp_values)*t_7_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_7_lower_ci <- t_7_lower_ci_adj
      amplification_results_ci$t_7_upper_ci <- t_7_upper_ci_adj
    }
    if(!is.na(single_time$t_8)){
      amplification_results_ci$t_8 <- mean(bootstrap_amplification$t_8, na.rm = TRUE)
      t_8_lower_ci <- confint(lm(t_8 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_8_lower_ci_adj <- ((length(tmp_values)*t_8_lower_ci)/(5 + length(tmp_values)))
      t_8_upper_ci <- confint(lm(t_8 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_8_upper_ci_adj <- ((5 + (length(tmp_values)*t_8_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_8_lower_ci <- t_8_lower_ci_adj
      amplification_results_ci$t_8_upper_ci <- t_8_upper_ci_adj
    }
    if(!is.na(single_time$t_9)){
      amplification_results_ci$t_9 <- mean(bootstrap_amplification$t_9, na.rm = TRUE)
      t_9_lower_ci <- confint(lm(t_9 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_9_lower_ci_adj <- ((length(tmp_values)*t_9_lower_ci)/(5 + length(tmp_values)))
      t_9_upper_ci <- confint(lm(t_9 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_9_upper_ci_adj <- ((5 + (length(tmp_values)*t_9_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_9_lower_ci <- t_9_lower_ci_adj
      amplification_results_ci$t_9_upper_ci <- t_9_upper_ci_adj
    }
    if(!is.na(single_time$t_10)){
      amplification_results_ci$t_10 <- mean(bootstrap_amplification$t_10, na.rm = TRUE)
      t_10_lower_ci <- confint(lm(t_10 ~ 1, bootstrap_amplification), level = 0.95)[1]
      t_10_lower_ci_adj <- ((length(tmp_values)*t_10_lower_ci)/(5 + length(tmp_values)))
      t_10_upper_ci <- confint(lm(t_10 ~ 1, bootstrap_amplification), level = 0.95)[2]
      t_10_upper_ci_adj <- ((5 + (length(tmp_values)*t_10_upper_ci))/(5 + length(tmp_values)))
      amplification_results_ci$t_10_lower_ci <- t_10_lower_ci_adj
      amplification_results_ci$t_10_upper_ci <- t_10_upper_ci_adj
    }
    
  }
  
  return(amplification_results_ci)
}

