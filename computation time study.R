## computation time study for the Binary Segmentation algorithm

# With thanks to:
# 
# The Tail Greedy unbalanced haar paper
# and
# https://stackoverflow.com/questions/6262203/measuring-function-execution-time-in-r

# We measure the compute time of the algorithm a process X 
# where each X_i is distributed with a standard normal distribution

study <- function(max.size){
  
  
  time.taken <- NULL
  process_lengths <- 10^(1:max.size)
  
  for( i in length(process_lengths)){
    
    print(i)
    test <- rnorm(process_lengths[i])
    
    start.time <- Sys.time()
    Binary.Segmentation(test)
    end.time <- Sys.time()
    
    time.taken <- c(time.taken, end.time - start.time)
  }
  
  return(time.taken)
}