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

library('changepoint')
library('stepR')
library('FDRSeg')
library('wbs')


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

store <- list(designation = "This is a list of the testing metrics for the simulation study")
modelList <- list("blocks", "fms", "teeth10", "stairs10", "mix")

MCsimulation <- function(li1, modelList, li3, N){
  
  store <- NULL
  modelList <- c("blocks", "fms", "teeth10", "stairs10", "mix")
  zero.triple <- c(0, 0, 0)
  
  for(i in 1:length(ModelList)){
    
    difference_in_changepoints <- rbind(zero.triple, zero.triple, zero.triple, zero.triple, zero.triple)
    difference_in_means <- c(0,0,0,0,0)
    
    for(j in 1:N){
     
      
      obs <- (mosum::testData(model = modelList[i]))
      
      sd <- sdrobnorm(obs)
      sigma <- mad(obs)
      qfs <- simulQuantile(alpha = 0.5, length(obs), type = "fdrseg")
      dfs <- simulQuantile(alpha = 0.5, length(obs), type = "smuce") 
      
      # using alpha = 0.5 here not sure if that is appropiate
      
      obs.wbs <- wbs(obs)
      obs.PELT <- changepoint::cpt.mean(obs*(1/sigma), method = "PELT")
      obs.FDRSeg <- fdrseg(obs, qfs, sd = sd)
      obs.myBinSeg <- Binary.Segmentation(obs)
      obs.SMUCE <- stepFit(obs, q = dfs)
      
 # The PELT and SMUCE methods use a different definition of changepoint to our other models
 # so we must change it
        obs.PELT@cpts <- obs.PELT@cpts -rep(1, length(obs.PELT@cpts))
        obs.SMUCE$rightEnd <- obs.SMUCE$rightEnd - rep(1, length(obs.SMUCE$rightEnd))
 # We want to compare the model's estimations of the number of changepoints
 # against the true number of changepoints.
        
      true_changepoints <- which(abs(diff(mosum::testSignal(model = modelList[i])$mu_t)=!0))
      no.true_changepoints <- length(true_changepoints)
      
      A <- determine(obs.wbs$cpt$no.cpt.th, no.true_changepoints)
      B <- determine(length(obs.PELT@cpts), no.true_changepoints)
      C <- determine(length(obs.FDRSeg$value), no.true_changepoints)
      D <- determine(length(obs.myBinSeg$est.means), no.true_changepoints)
      E <- determine(length(obs.SMUCE$rightEnd), no.true_changepoints)
      
      difference_in_changepoints <- unname(difference_in_changepoints) + rbind(A, B, C, D, E)
      
 # We also want to compare the estimated model against the true model.
      
      true_signal <- mosum::testSignal(model = modelList[i])$mu_t
      
      A <- sum(abs(means.between.cpt(obs, cpt = obs.wbs$cpt$cpt.th[[1]]) - true_signal))/length(obs)
      B <- sum(abs(rep(obs.PELT@param.est$mean, obs.PELT@cpts) - true_signal))/length(obs)
      C <- sum(abs(rep(obs.FDRSeg$value, diff(c(obs.FDRSeg$left, obs.FDRSeg$n)) - true_signal)))/length(obs)
      D <- sum(abs(obs.myBinSeg$est.signal - true_signal))/length(obs)
      E <- sum(abs(rep(obs.SMUCE$value, diff(c(obs.SMUCE$leftEnd, length(true_signal)))) - true_signal))/length(obs)
      
      difference_in_means <- difference_in_means + rbind(A, B, C, D, E)
      
    }
    
    foo <- difference_in_changepoints
    a <- c(-1,0,1)
    const <- 1/N
    mean_difference <- rbind(const*sum(foo[1,]* a), const*sum(foo[2,]* a), const*sum(foo[3,]* a), const*sum(foo[4,]* a))
    var_difference <- rbind(const*sum((foo[1,]*a)^2), const*sum((foo[2,]*a)^2), const*sum((foo[3,]*a)^2), const*sum((foo[4,]*a)^2))
    
    bar <- difference_in_means
    Mu_difference <- bar * rep(const, length(bar))
    
    interim <- cbind(difference_in_changepoints, mean_difference, var_difference, Mu_difference)
    
    store <- c(store, interim)
  }
  names(store) <- modelList
  return(store)
}
