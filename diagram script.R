## Louis Forbes Wright

#using _-_- as a basis for the algorithm created here


# we start by creating some data which contains change points.

means <- c(0,3,4,5,10,4,0)
changePoints <- c(0,100,200,300,400,500,600,700)

#means <- c(rnorm(n = 15, mean = 20, sd = 100), rep(0,135))
#changePoints <- c(rchisq(n = 15, df = 130), rep(0,135))

data <- NULL
 
for (i in 1:length(changePoints)-1){
  data <- c(data, rnorm(n = (changePoints[i+1] - changePoints[i]), mean = means[i], sd = 0))
  } 
par(mfrow = c(2,1))

plot(0:(length(data)-1), data, type = 'l')
abline(h = mean(means)) ### NEED TO ADD CHANGEPOINTS TO GRAPH.

# now we need to create the function defined in the Introduction
#to show the what the test statistic looks like

partialSum <- function(sour, len){
  sum(sour[1:len])
}

Z_k <- function(sour, k_1, k_2){
  #this function creates and returns all the test score values in a given interval.
  
  test.stat <- NULL
  
  for (k in k_1+1:k_2-1){
  
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
output <- Z_k(data, 0,length(data)) # we create a vector of all the test scores

#estChangePoints <- NULL
#estChangePoints <- c(estChangePoints, max.col(test.stat))

plot(0:(length(output)-1), output, type = 'l') # we plot the test scores against the index of that data
boundary <- length(output)^(3/8)               # the critical value for the size of the interval tested
abline(h = boundary)