


Z_k <- function(data, k_1, k_2){
  #this function creates and returns all the test score values in a given interval.
  
  test.stat <- NULL
  cx <- cumsum(data)
  for (k in (k_1+1):(k_2-1)){
    
    #for each k in the interval we create a test score, 
    #|z_k|
    
    coeff <- (k - k_1)/(k_2 - k_1)
    max_difference <- (cx[k_2]-cx[k_1])
    const <- coeff*max_difference  
    
    top_fraction <- const - (cx[k]-cx[k_1])
    bottom_fraction <- sqrt((k - k_1)*(1-(coeff)))
    test.stat <- c(test.stat, abs((top_fraction)/(bottom_fraction)))
    
  }
  return(test.stat)
}  

construct.model <- function(sampl, ofs){
  
  est.signal <- NULL
  means <- NULL
  
  for(i in 1:(length(ofs)-1)){
    a <- ofs[i]
    b <- ofs[i+1]
    print(a)
    print(b)
    est.signal <- c(est.signal, rep(mean(sampl[a:b]), (b-a)))
    print(est.signal)
    means <- c(means, mean(sampl[a:b]))
    }
  est.signal <- c(est.signal[1],est.signal)
  return(list(est.signal=est.signal, est.means=means))
}

Binary.segmentation.operation <- function(sample, passed_segments = NULL){
  
  changepoint_scores <- NULL
  estimated_changepoints <- NULL
  finalised_segments <- NULL
  boundary <- length(sample)^(3/8)
  
  if(is.null(passed_segments)){ #the first operation on a dataset will require initalisation
    working_segments <- c(1, length(sample))
    finalised_segments <- c(1,length(sample))
  } else {
      working_segments <- sort(passed_segments)
      finalised_segments <- passed_segments
  }
  
  
  
  #we need to order the estimated change-points so that the largest Z_k is first
  for (i in 1:(length(working_segments)-1)){ # n points gives n-1 segments
    if(max(Z_k(sample, working_segments[i],working_segments[i+1]), na.rm = TRUE)>= boundary){
      #we take the maximum as a estimated change point
     if(!((which.max(Z_k(sample, working_segments[i],working_segments[i+1]))+working_segments[i]) %in% passed_segments)){
       #we test if the changepoinit is already accounted for 
       estimated_changepoints <- c(estimated_changepoints, (which.max(Z_k(sample, working_segments[i],working_segments[i+1]))+working_segments[i]))
       changepoint_scores <- c(changepoint_scores, max(Z_k(sample, working_segments[i],working_segments[i+1]), na.rm = TRUE))
     }
    } 
  }
  print(estimated_changepoints)
  if(!(is.null(estimated_changepoints))){
    #print(estimated_changepoints)
    #print(changepoint_scores)
    bound <- cbind(estimated_changepoints, changepoint_scores)
    #print(bound)
    estimated_changepoints <- bound[order(-bound[,2]), 1]
    #print(estimated_changepoints)
  }
  finalised_segments <- c(finalised_segments, estimated_changepoints)
  return(unname(finalised_segments))
}


Binary.Segmentation <- function(sampl){
  held_segments <- (1)
  passed_segments <- NULL
  est.signal <- NULL
  
  passed_segments <- Binary.segmentation.operation(sampl)
  while(!setequal(held_segments, passed_segments)){
    # binary segmentation runs until the operation doesn't produce a new change-point
    held_segments <- passed_segments
    passed_segments <- Binary.segmentation.operation(sampl, held_segments)

  }
  print("in order of appearance the estimated changepoints are ")
  ordered_finalised_segments <- sort(passed_segments)
  #ordered_finalised_segments <- optimise.model(sampl,passed_segments)
  print(ordered_finalised_segments)
  output <- construct.model(sampl, ordered_finalised_segments)
  ofs <- list(est.changepoints = ordered_finalised_segments)
  output <- c(output, ofs)
  
  return(output)
}
