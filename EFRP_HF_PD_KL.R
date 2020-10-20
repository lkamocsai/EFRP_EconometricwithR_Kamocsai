rm(list = ls()) # clear all variables from the workspace
graphics.off()  # clear all plots

#function for store the parameters of DGP

#ARc      = coef of AR
#MAc      = coef of MA
#sd       = standard deviation
#tsLength = length of time series
#tsNumber = number of time series

dgp_parameters <- function(ARc, MAc, sd, tsLength, tsNumber)
{
  
  x <- list("ARc" = ARc ,"MAc" = MAc ,"sd" = sd,"tsLength" = tsLength,"tsNumber" = tsNumber)
  return(x)
}

####dgpParams <- dgp_parameters("ARc" = c(0.2, 0.3), "MAc" = c(0.4, 0.5), "sd" = 0.2, "tsLength" = 10, "tsNumber" = 2)####


#function for generating timeseries


arimas <- function(dgpParams)
{
  result = matrix(0,dgpParams$tsNumber,dgpParams$tsLength)
  for (i in (1:dgpParams$tsNumber)){
    result[i,1:dgpParams$tsLength]=arima.sim(n = dgpParams$tsLength, list(order=c(length(dgpParams$ARc),0,length(dgpParams$MAc)),
                                   ar = dgpParams$ARc, ma = dgpParams$MAc), sd = dgpParams$sd)
  }
  return(result)
}

#####test = arimas(dgpParams)#####

#function for generating time series of 4x4 grid

#n <- c(5,10,15,20)                       
#sigma <- c(0.01,0.1,0.02,0.2)
#ARc <- c(0.2, 0.3)                       #####?????????????
#MAc <- c(0.4, 0.5)
#tsNumber <- 2


timeseries_grid <- function(n, sigma, ARc, MAc, tsNumber){
  result <- matrix(list(), length(n), length(sigma))
  
  
  for (i in (1:length(n))){
    for (j in (1:length(sigma))) {
      dgpParams <- dgp_parameters("ARc" = ARc, "MAc" = MAc, "sd" = sigma[j], "tsLength" = n[i], "tsNumber" = tsNumber)
      result[[i,j]] = arimas(dgpParams)
      # first index is for tsLength,  second index is for standard deviation
    }
  }
  return(list(result=result,n=n,sigma=sigma))
}

####x=timeseries_grid(n,sigma,ARc,MAc,tsNumber)#####


#3 sets of different complexity
#case No1

n <- c(10,50,100,200)
sigma <- c(0.01,0.1,1,2)
ARc1 <- c(0.2, 0.4)
MAc1 <- c(0.1, 0.2)
tsNumber <- 100       #not 1000, because lack of time (1000 is running for too long)

finalTimeSeries1 = timeseries_grid(n,sigma,ARc1,MAc1,tsNumber)

#case No2

ARc2 <- c(0.3, 0.05)
MAc2 <- c(0.4, 0.2, 0.5)

finalTimeSeries2 = timeseries_grid(n,sigma,ARc2,MAc2,tsNumber)

#case No3

ARc3 <- c(-0.1, 0.3, 0.1)
MAc3 <- c(0.3, 0.1, 0.2)

finalTimeSeries3 = timeseries_grid(n,sigma,ARc3,MAc3,tsNumber)



#function for fitting models 

fit_model_on_series <- function(inputSeries, orderToFit){
  
  tmp = ts(inputSeries)
  result = matrix(0,dim(tmp)[1],2)
  
  for (i in (1:dim(tmp)[1])){
    
    dataToFit = tmp[i,1:dim(tmp)[2]]
    model <- arima(dataToFit, order=orderToFit, method="ML")
    BIC=AIC(model,k = log(length(dataToFit)))
    
    result[i,1] = model$aic
    result[i,2] = BIC
  }
  return(result)
}


#function for fitting on all models

fit_model_on_grid <- function(dataGrid,orderToFit){
  
  result <- matrix(list(), length(dataGrid$n), length(dataGrid$sigma))
  
  for(i in (1:length(dataGrid$n))){
    for(j in (1:length(dataGrid$sigma))){
      result[i,j] = list(fit_model_on_series(dataGrid$result[i,j][[1]],orderToFit))
      cat("Elvégezve:",((4*(i-1)+j)/16*100), "%\n")
    }
  }
  return(result)
}



#fitting for the first parameter setting
#case No1

orderToFit1 = c(1,0,1)
model11 = fit_model_on_grid(finalTimeSeries1, orderToFit1)
orderToFit2 = c(2,0,2)
model21 = fit_model_on_grid(finalTimeSeries1, orderToFit2)
orderToFit3 = c(4,0,4)
model31 = fit_model_on_grid(finalTimeSeries1, orderToFit3)

#fitting for the first parameter setting
#case No2

orderToFit1 = c(1,0,1)
model12 = fit_model_on_grid(finalTimeSeries2, orderToFit1)
orderToFit2 = c(2,0,3)
model22 = fit_model_on_grid(finalTimeSeries2, orderToFit2)
orderToFit3 = c(3,0,4)
model32 = fit_model_on_grid(finalTimeSeries2, orderToFit3)

#fitting for the first parameter setting
#case No2

orderToFit1 = c(1,0,1)
model13 = fit_model_on_grid(finalTimeSeries3, orderToFit1)
orderToFit2 = c(3,0,3)
model23 = fit_model_on_grid(finalTimeSeries3, orderToFit2)
orderToFit3 = c(4,0,4)
model33 = fit_model_on_grid(finalTimeSeries3, orderToFit3)


#function for choosing the best model from 3 different

choose_best_model <- function(mRes1, mRes2, mRes3){
  
  bestModel = mRes1
  
  for (k in (1:dim(mRes1)[2])) {
    
    for (j in (1:dim(mRes1)[1])){
      
      for (i in (1:(length(mRes1[1,1][[1]])/2))){
        bestModel[j,k][[1]][i,1] = which.min(c(mRes1[j,k][[1]][i,1],mRes2[j,k][[1]][i,1],mRes3[j,k][[1]][i,1]))
        bestModel[j,k][[1]][i,2] = which.min(c(mRes1[j,k][[1]][i,2],mRes2[j,k][[1]][i,2],mRes3[j,k][[1]][i,2]))
      }
    }
  }
  return(bestModel)
}

results1 = choose_best_model(model11,model21,model31)
results2 = choose_best_model(model12,model22,model32)
results3 = choose_best_model(model13,model23,model33)

#function for counting for each model how many times were they the best

count_best_model <- function(bestModel){
  
  res = matrix(list(),dim(bestModel)[1],dim(bestModel)[2])
  
  for (k in (1:dim(bestModel)[2])) {
    
    for (j in (1:dim(bestModel)[1])){

        db1AIC = sum(bestModel[j,k][[1]][1:(length(bestModel[1,1][[1]])/2),1]==1)
        db2AIC = sum(bestModel[j,k][[1]][1:(length(bestModel[1,1][[1]])/2),1]==2)
        db3AIC = sum(bestModel[j,k][[1]][1:(length(bestModel[1,1][[1]])/2),1]==3)
        db1BIC = sum(bestModel[j,k][[1]][1:(length(bestModel[1,1][[1]])/2),2]==1)
        db2BIC = sum(bestModel[j,k][[1]][1:(length(bestModel[1,1][[1]])/2),2]==2)
        db3BIC = sum(bestModel[j,k][[1]][1:(length(bestModel[1,1][[1]])/2),2]==3)
        res[j,k] = list(c(db1AIC,db2AIC,db3AIC,db1BIC,db2BIC,db3BIC))
    }
  }
  return(res)
  
}

AIC_BIC1 = count_best_model(results1)  #e.i. AIC_BIC[1,1]  gives 6 numbers: first 3 numbers: 1./2./3. model were the best x times by the lowest AIC
                                                                          #second 3 numbers: 1./2./3. model were the best x times by the lowest BIC

AIC_BIC2 = count_best_model(results2)
AIC_BIC3 = count_best_model(results3)

