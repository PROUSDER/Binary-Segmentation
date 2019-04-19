
optimise.model <- function(data, modelpar){
  # modelpar is the location of changepoints such that the entries are sorted 
  # descendingly by round of binary segmetation and then descendingly by stat.score size.
  # Optimises model using BIC
  
  combined.BIC.score <- NULL
  size <- length(data)
  
  
  
  number.of.models <- length(modelpar)
  for(i in 2:number.of.models){
    
    means <- NULL
    changepoints_tested <- sort(modelpar[1:i])
    segments <- diff(changepoints_tested)
    for(j in 1:(length(changepoints_tested)-1)){
      a <-changepoints_tested[j]
      b <-changepoints_tested[j+1]
      #print(a)
      #print(b)
      means <- c(means, (mean(data[a:b])))
    }
    
    model.to.be.tested <- cbind(segments, means)
    
    #print(model.to.be.tested)
    est.signal <- rep(means, times = segments)
    est.signal <- c(est.signal[1], est.signal)
    #print(est.signal)
    RSS <- sum((data - est.signal)^2)
    #print(RSS)
    BIC.score <- (size/2) * log((1/size)*RSS) + i*log(size)
    combined.BIC.score <- c(combined.BIC.score, BIC.score)
  }
  #print(combined.BIC.score)
  best.model.size <- which.min(combined.BIC.score)+1
  #print(best.model.size)
  best.model <- modelpar[1:best.model.size]
 return(sort(best.model))
}

cx <- cumsum(data)


Z_k <- function(sour, k_1, k_2){
  #this function creates and returns all the test score values in a given interval.
  
  test.stat <- NULL
  
  cx <- cumsum(sour)
  
  for (k in (k_1+1):(k_2-1)){
    
    #for each k in the interval we create a test score, 
    #|z_k|
    
    coeff <- (k -k_1)/(k_2 - k_1)
    max_difference <- (cx[k_2]-cx[k_1])
    const <- coeff*max_difference  
    
    top_fraction <- const - (cx[k]-cx[k_1])
    bottom_fraction <- sqrt((k - k_1)*(1-(coeff)))
    test.stat <- c(test.stat, abs((top_fraction)/(bottom_fraction)))
    
  }
  return(test.stat)
}  

OBSE <- function(data){
  est.error <- NULL
  for(i in 1:1000){
  simulated_values <- round(runif(round(length(data)/20), min = 0, max = length(data)))
  simulated_values <- data[simulated_values]
  est.error <- c(est.error, max(Z_k(simulated_values, 1, round(length(data)/20))))
  } 
  len <- length(est.error)
  significant.quartile <- round(len*0.99)
  est.error[significant.quartile]
}

variance.est <- function(data){
  # using 'first lag MAD, 
  #using the first lag will remove the effect of the signal from the data'
  variance.estimate <- mad(diff(data))
  return(variance.estimate)
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

static_bootstrap_error <- function(){
  stored_errors <- NULL
  
  for(j in 1:100){
    stored_errors <- c(stored_errors, Bootstrap.of.error(rnorm(j * 50)))
    print(j)
  }
  return(stored_errors)
}