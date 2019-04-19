##
#' A function for creating the diagram for Lemma 3.2
#'
#' @param n the function creates n changepoints with the partition of [0,1]
#' 


diagram_script <- function(n, m = 1000){
  
  par(mfrow = c(2,1))
  
  x <- seq(0,1,1/m)
  scale <- sqrt(x*(1-x))
  y <- rep(0, length(x))
  
  changepoints <- trunc(sort(c(0, runif(n), 1))*m)/m
  scalar <- diff(changepoints)
  # we aim to create a step function such that over when we integrate the step function over x we get 0
  
  # to do this we first set step values such that the sum of the step values = 0
  step_value <-  c(rep(1/length(scalar),length(scalar)/2),rep(-1*(1/length(scalar)),length(scalar)/2))
  # an line to ensure that the number of step values is equal to number of distances between changepoints
  step_value <- if(length(step_value) != length(scalar)){
    c(step_value[1:(length(step_value)/2)] ,0 , step_value[(length(step_value)/2 + 1):length(step_value)])
  }else{step_value}
  #we then scale each step value so that when multiplied with the length of the step we obtain the original value
  # so that the sum is still 0
  step_value <- step_value/scalar
  lambda <- rep(step_value, diff(trunc(changepoints * m)))
  
  cx <- cumsum(lambda)
  
  f <- function(j){
    cx[j]/scale[j]
  }
  
  out <- f(2:(length(x)-1))
  
  plot(lambda, type = 'l', xaxp = c(0,0,1), xlab =  TeX("$$"), ylab = TeX("$$"))
  axis(side = 1, at = seq(0,10,2)*100, labels = seq(0,10,2) * 0.1)
  plot(out, type = 'l', xaxp = c(0,0,1), ylab = TeX("$f(x)$"), xlab = TeX("$x$"))
  axis(side = 1, at = seq(0,10,2)*100, labels = seq(0,10,2) * 0.1)
}
