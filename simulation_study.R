## simulation study
#
#
#inputs: [list of models, list of methods, list of method parameters ]
#outputs: tables for each model [method , <N N >N, Exp[N_hat - N], MSE, d(mu, mu*)]
#
# We will be using the methods: WBS, BS, PELT, SMUCE, and FDRSeg. We are aiming to use
# a balance of both iterated single change-point as well as global methods in the study
#
# To generate the data, we will be using the models: blocks, fms, mix, teeth10, stairs10, and extreme teeth,
# with a selection of parameters

# thanks goes to () for their work in the implementation and the associated papers ()
# as well as the 'mosum'

library('changepoint')
library('stepR')
library('FDRSeg')
library('wbs')
library('mosum')

#blocks model
blocks <- function(variance = 1,changepoints, means, size){
  if(all(diff(changepoints)>0)){
    corrected.changepoints <- c(0, changepoints, size)
  }else{
    print("changepoints must be linearly ordered")
  }
  
  error <- rnorm(n = size, sd=sqrt(variance))
  if(((length(means)-1) == length(changepoints))){
    signal <- rep(means, times = diff(corrected.changepoints))
  }
  model <- signal + error
  return(model)
}
determine <- function(test_value, true_value){
  foo <- test_value - true_value
  logical <- c(foo < 0, foo == 0, foo > 0)
  
  return(as.numeric(logical))
}

handle <- function(foo){
  output <- NULL
  sample.mean <- NULL
  sample.var <- NULL
  a <- c(-1,0,1)
  N <- length(foo[1,])
  
  for(i in 1:length(foo[,1])){
    sample.mean <- c(sample.mean, mean(rep(a, times =  foo[i,])))
    
    sample.var <-c(sample.var, var(rep(a, times =  foo[i,])))
  }
  return(cbind(sample.mean, sample.var))
}


store <- list(designation = "This is a list of the testing metrics for the simulation study")
modelList <- c("blocks", "fms", "teeth10", "stairs10", "mix")

MCsimulation <- function(models = NULL, N){
  
  modelList <- models
  
  store <- list(designation = "This is a list of the testing metrics for the simulation study")
  if(is.null(models)){modelList <- c("blocks", "fms", "teeth10", "stairs10", "mix")}
  zero.triple <- c(0, 0, 0)
  
  
  for(i in 1:length(modelList)){
    
    difference_in_changepoints <- rbind(zero.triple, zero.triple, zero.triple, zero.triple, zero.triple)
    difference_in_means <- c(0,0,0,0,0)
    changepoint.count <- NULL
    
    obs <- mosum::testData(model = modelList[i])
    qfs <- simulQuantile(alpha = 0.5, length(obs), type = "fdrseg")
    dfs <- simulQuantile(alpha = 0.5, length(obs), type = "smuce") 
    
    for(j in 1:N){
      
      if((j %in% (ceiling(1:N/2)*2)) | j <= 10){
        counter <- c(i,j)
        print(counter)
        
      }
     
      obs <- (mosum::testData(model = modelList[i]))
      
      sd <- sdrobnorm(obs)
      sigma <- mad(obs)
      
      # using alpha = 0.5 here, as suggested by the authors of the package
      
      obs.wbs <- wbs(obs)
      
      # the PELT method requires the data to variance to be standardised
      obs.PELT <- changepoint::cpt.mean(obs*(1/sd), method = "PELT")
      obs.FDRSeg <- fdrseg(obs, qfs, sd = sigma) 
      obs.myBinSeg <- Binary.Segmentation(obs)
      obs.SMUCE <- stepFit(obs, q = dfs)
      
 # The PELT and SMUCE methods use a different definition of changepoint to our other models
 # so we must change it
        #obs.PELT@cpts <- obs.PELT@cpts - 1
        obs.SMUCE$leftEnd <- obs.SMUCE$leftEnd - 1
 # We want to compare the model's estimations of the number of changepoints
 # against the true number of changepoints.
        
      true_changepoints <- which(abs(diff(mosum::testSignal(model = modelList[i])$mu_t))!=0)
      no.true_changepoints <- length(true_changepoints)
      
      A <- determine(obs.wbs$cpt$no.cpt.th, no.true_changepoints)
      B <- determine(length(obs.PELT@cpts), no.true_changepoints)
      C <- determine(length(obs.FDRSeg$value), no.true_changepoints)
      D <- determine(length(obs.myBinSeg$est.means), no.true_changepoints)
      E <- determine(length(obs.SMUCE$leftEnd), no.true_changepoints)
      
      difference_in_changepoints <- unname(difference_in_changepoints) + rbind(A, B, C, D, E)
      changepoint.count <- rbind(changepoint.count, cbind(obs.wbs$cpt$no.cpt.th, length(obs.PELT@cpts), length(obs.FDRSeg$value), length(obs.myBinSeg$est.means), length(obs.SMUCE$leftEnd)))
      
 # We also want to compare the estimated model against the true model.
 # some of the results require preprocessing
      
      obs.FDRSeg$left[1] <- 0
      obs.FDRSeg$left <- c(obs.FDRSeg$left, obs.FDRSeg$n)
      
      true_signal <- mosum::testSignal(model = modelList[i])$mu_t
      
      A <- sum(abs(means.between.cpt(obs, cpt = obs.wbs$cpt$cpt.th[[1]]) - true_signal))/length(obs)
      B <- sum(abs(rep(obs.PELT@param.est$mean, diff(c(0,obs.PELT@cpts))) - true_signal))/length(obs)
      
      h <- rep(obs.FDRSeg$value, times = diff(obs.FDRSeg$left))
      j <- abs(h - true_signal)
      C <- sum(j)/length(obs)
      
    
      D <- sum(abs(obs.myBinSeg$est.signal - true_signal))/length(obs)
      E <- sum(abs(rep(obs.SMUCE$value, diff(c(0, obs.SMUCE$rightEnd))) - true_signal))/length(obs)
    
      
      difference_in_means <- difference_in_means + rbind(A, B, C, D, E)
      
      
    }
    
    foo <- changepoint.count
    #print(foo)

    mean_difference <- rbind(mean(foo[,1]), mean(foo[,2]), mean(foo[,3]), mean(foo[,4]), mean(foo[,5]))
    var_difference <- rbind(var(foo[,1]), var(foo[,2]), var(foo[,3]), var(foo[,4]), var(foo[,5]))
    
    
    pass_in <- cbind(mean_difference, var_difference)
    
    MAE <- difference_in_means / N
    
    interim <- cbind(difference_in_changepoints, pass_in, MAE)

    store <- c(store, list(interim))
    print(interim)
  }
  names(store) <- c("designation", modelList)
  return(store)
}
