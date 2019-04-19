


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

Bootstrap.of.error <- function(data, percentile = 99, simulation.number = 1000){
  #we aim to find what effect we can expect the error to have on the CUSUM (Z_k) we are using, 
  #assuming the error is standard normal. Then we take the 99%-ile of the simulated CUSUMs
  #and use this to determine if a the maximum CUSUM score of a binary-segmentation operation is significant.
  
  
  max.scores.attained <- NULL
  size <- length(data)
  for(i in 1:simulation.number){
    max.scores.attained <- c(max.scores.attained, max(Z_k(rnorm(size), 1, size)))
  }
  
  high.percentile <- quantile(max.scores.attained, percentile/100)
  return(high.percentile)
}


construct.model <- function(sampl, ofs){
  
  est.signal <- NULL
  means <- NULL
  
  for(i in 1:(length(ofs)-1)){
    a <- ofs[i]
    b <- ofs[i+1]

    est.signal <- c(est.signal, rep(mean(sampl[a:b]), (b-a)))
    
    means <- c(means, mean(sampl[a:b]))
    }
  est.signal <- c(est.signal[1],est.signal)
  return(list(est.signal=est.signal, est.means=means))
}

Binary.segmentation.operation <- function(sampl, passed_segments = NULL, boundary){

  changepoint_scores <- NULL
  estimated_changepoints <- NULL
  finalised_segments <- NULL

  
  if(is.null(passed_segments)){ #the first operation on a dataset will require initalisation
    working_segments <- c(1, length(sampl))
    finalised_segments <- c(1,length(sampl))
  } else {
      working_segments <- sort(passed_segments)
      finalised_segments <- passed_segments
  }
  

  
  #we need to order the estimated change-points so that the largest Z_k is first
  for (i in 1:(length(working_segments)-1)){ # n changepoints gives n-1 segments
    #we take the maximum as a estimated change point
    # the zero removes the error message, but has no affect on the output of the function
    if(max(c(Z_k(sampl, working_segments[i],working_segments[i+1]), 0), na.rm = TRUE)>= boundary){
      ch <- which.max(Z_k(sampl, working_segments[i],working_segments[i+1]))+working_segments[i]
      #we test if the changepoinit is already accounted for
     if(!( ch %in% passed_segments)){
        
       estimated_changepoints <- c(estimated_changepoints, ch)
       changepoint_scores <- c(changepoint_scores, max(Z_k(sampl, working_segments[i],working_segments[i+1]), na.rm = TRUE))
     }
    } 
  }
 
  if(!(is.null(estimated_changepoints))){
    
    bound <- cbind(estimated_changepoints, changepoint_scores)
    
    estimated_changepoints <- bound[order(-bound[,2]), 1]
   
  }
  finalised_segments <- c(finalised_segments, estimated_changepoints)
  return(unname(finalised_segments))
}


Binary.Segmentation <- function(sampl, boundary = NULL, Bootstrap.settings = c(99 ,1000)){
  #
  # the function will apply binary segmentation to the supplied data.
  # the function will use stored values for data with length less than 5000. 
  #
  # inputs   {sample data, significant value for the CUSUM}
  #         @ sampl, boundary
  #
  # returns  {estimated signal, estimated means, estimated changepoints}
  #         @ $est.signal, $est.means, $est.changepoints
  #
  #
  
  held_segments <- (1)
  passed_segments <- NULL
  est.signal <- NULL
  
      # as var(X_t - X_t-1) = 2var(X_t) 
  sigma <- mad(diff(sampl)/sqrt(2)) 
  
  if(is.null(boundary)){
    if(ceiling(length(sampl)/50)>length(gr)){
      significant_value <- Bootstrap.of.error(sampl, Bootstrap.settings[1], Bootstrap.settings[2]) * sigma
      
      # the differencing in the median absolute deviation removes a linear trend if it exists from the data.
      # Further improvements of the binary segmentation operation
      # could see the addition of some methods to control for
      # other trend structures in data (seasonality, sinusoidal models, etc)
      
      break
      }else{significant_value <- gr[ceiling(length(sampl)/50)] * sigma}
    }else{
    significant_value <- boundary
  }
  #print(paste("using a significant value of",significant_value))
  passed_segments <- Binary.segmentation.operation(sampl, boundary = significant_value)
  while(!setequal(held_segments, passed_segments)){
    # binary segmentation runs until the operation doesn't produce a new change-point
    held_segments <- passed_segments
    passed_segments <- Binary.segmentation.operation(sampl, held_segments, significant_value)

  }
  #print("in order of occurence the estimated changepoints are ")
  ordered_finalised_segments <- sort(passed_segments)

  #print(ordered_finalised_segments)
  output <- construct.model(sampl, ordered_finalised_segments)
  ofs <- list(est.changepoints = ordered_finalised_segments)
  output <- c(output, ofs)
  
  return(output)
}

printOut <- function(script){
  nth<-paste0(script, c(rep(",", length(script))))
  print(nth)
}  
  
  
 gr <-c( 
   #values calculated from static_bootstrap_error
   #need to use these to create a expected coefficient for each length of data      
   1.59355109511403, 2.27434764081231, 1.80701533003434, 2.8033070905804,  2.47273118124658,
   2.2678094766454,  1.70293648265862, 1.75824435066782, 2.39994853939748, 2.05618556035496,
   1.73291957577073, 3.11727859002214, 2.10249651986484, 2.67586029410928, 2.76878761718571,
   2.75785692702678, 4.87509641573421, 2.22512584062143, 2.03842524505992, 2.37116906755947,
   1.5474823207671,  2.55967603965411, 2.25032612965191, 2.24384714649027, 2.67298196948062,
   2.75206181348804, 2.3844308779272,  1.86951557520666, 2.90789763990639, 2.39609079059703,
   2.36630467691164, 1.65735179688503, 1.92602291823894, 3.16666079285206, 2.55407626276911,
   2.71336027282355, 2.70211719086013, 2.74081541080886, 2.10870397995379, 2.58190811840651,
   2.91100111458016, 2.5324076524782,  2.01231573699737, 1.80634053283447, 2.73435942159519,
   2.80245913502402, 2.64791378874629, 3.00301313866373, 2.69935421670575, 1.94230421909558,
   2.16407636957608, 1.44827011693701, 2.72442517260475, 1.99176484428726, 2.49710108434088,
   2.35418251227281, 2.57755205493521, 3.25433632142714, 2.14215706008829, 1.94472238634253,
   2.88476115509466, 1.82735593220973, 2.2774540267339,  1.67420231173741, 3.44292470484641,
   2.40759330752955, 2.68237696493158, 2.17008872899264, 2.79188401240115, 2.13069343600699,
   2.66100024669349, 1.80807501956531, 2.50112856124061, 2.0369523705705,  2.21318887583437,
   2.41189321549167, 2.99824201791627, 3.79203106205445, 1.67307440964618, 2.51447446833737,
   1.77702995440098, 1.5505891998448,  2.58992961159756, 2.07853452213866, 1.7625054324503, 
   2.62990767012543, 2.23781938670282, 2.62106475697716, 2.21957314044713, 3.04026112958293,
   1.88810850395291, 2.1071713219079,  2.69014396012727, 2.30362612067092, 2.46792439170224,
   2.67455884182336, 1.65455090250905, 3.87022184984458, 2.92549773370628, 2.16536099742458
   )
