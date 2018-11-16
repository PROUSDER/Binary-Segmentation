partialSum <- function(sour, len){
  sum(sour[1:len])
}

Z_k <- function(sour, k_1, k_2){
  #this function creates and returns all the test score values in a given interval.
  
  test.stat <- NULL
  
  for (k in (k_1+1):(k_2-1)){
    
    #for each k in the interval we create a test score, 
    #|z_k|
    
    coeff <- (k -k_1)/(k_2 - k_1)
    max_difference <- (partialSum(sour, k_2)-partialSum(sour, k_1))
    const <- coeff*max_difference  
    
    top_fraction <- const - (partialSum(sour, k)-partialSum(sour, k_1))
    bottom_fraction <- sqrt((k - k_1)*(1-(coeff)))
    test.stat <- c(test.stat, abs((top_fraction)/(bottom_fraction)))
    
  }
  return(test.stat)
}  

Binary.segmentation.operation <- function(sample, passed_segments = NULL){
  
  
  finalised_segments <- NULL
  boundary <- length(sample)^(3/8)
  
  if(is.null(passed_segments)){ #the first operation on a dataset will require initalisation
    working_segments <- c(1, length(sample))
  } else {
      working_segments <- passed_segments
  }
  finalised_segments <- working_segments
  
  for (i in 1:(length(working_segments)-1)){ # n points gives n-1 segments
    if(max(Z_k(sample, working_segments[i],working_segments[i+1]), na.rm = TRUE)>= boundary){
      #we take the maximum as a estimated change point
     if(!(which.max(Z_k(sample, working_segments[i],working_segments[i+1])) %in% passed_segments)){
       #we test if the changepoinit is already accounted for 
       finalised_segments <- c(finalised_segments, which.max(Z_k(sample, working_segments[i],working_segments[i+1])))}
    } 
  }
  finalised_segments <- sort(finalised_segments)
  return(finalised_segments)
}

Binary.Segmentation <- function(sample){
  held_segments <- (1)
  passed_segments <- NULL
  
  passed_segments <- Binary.segmentation.operation(sample)
  while(!setequal(held_segments, passed_segments)){
    # binary segmentation runs until the operation doesn't produce a new change-point
    held_segments <- passed_segments
    passed_segments <- Binary.segmentation.operation(sample, held_segments)

  }
  
  return(passed_segments)
}