## Louis Forbes Wright

#using _-_- as a basis for the algorithm created here


# we start by creating some data which contains a change point.

means <- c(0,3,4,5,10,4,0)
changePoints <- c(100,200,300,400,500,600,700)

j <- 1
data <- NULL
 
for (i in 1:length(means)){
  data <- c(data, rnorm(100 ,mean = means[i]))
  } 
par(mfrow = c(2,1))

plot(0:(length(data)-1), data, type = 'l')

# now we need to create the function defined in the Introduction
#to show the what the test statistic looks like

partialSum <- function(len){
  sum(data[1:len])
}

Z_k <- function(k, minimum, maximum){
  
  coeff <- (k -minimum)/(maximum - minimum)
  (coeff*(partialSum(maximum)-partialSum(minimum) - (partialSum(k)-partialSum(minimum))/(sqrt((k - minimum)*(1 - (coeff))))))
}
test.stat <- NULL
for (i in 0:length(data)){
  test.stat <-c(test.stat, Z_k(i, 0, length(data)))
}

plot(0:(length(test.stat)-1), test.stat, type = 'l')
