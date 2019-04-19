## Louis Forbes Wright

#using _-_- as a basis for the algorithm created here


# we start by creating some data which contains change points.

means <- c(0,1)
changepoints <- c(0,50,100)
  
data <- rep(means, diff(changepoints))
par(mfrow = c(2,1))

plot(0:(length(data)-1), data, type = 'l', xlab ="i", ylab = "a_i")
abline(h = mean(means)) ### NEED TO ADD CHANGEPOINTS TO GRAPH.

# now we need to create the function defined in the Introduction
#to show the what the test statistic looks like


cx <- cumsum(data)


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
  test.stat <- c(0, test.stat, 0)
  return(test.stat)
}  
output <- Z_k(data, 1,length(data)) # we create a vector of all the test scores

#estChangePoints <- NULL
#estChangePoints <- c(estChangePoints, max.col(test.stat))

plot(0:(length(output)-1), output, type = 'l', xlab = "x", ylab ="f(x)") # we plot the test scores against the index of that data
