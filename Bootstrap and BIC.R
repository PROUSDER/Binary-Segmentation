## Residuals for binary-segmentation

error <- (data - signal)^2



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

Bootstrap.of.error <- function(data){
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

variance.est <- function(data){# using 'MAD'
  x.tilde <- median(data)
  centralised.x <- data - x.tilde
  variance.estimate <- median(abs(centralised.x))
  
  return(variance.estimate)
}