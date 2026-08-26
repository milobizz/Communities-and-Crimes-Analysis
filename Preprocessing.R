library(tidyverse)
library(rsample)
library(FNN)
library(igraph)
library(datasets)
library(maps)
library(usmap)
library(GGally)
library(paletteer)
library(tree)
library(pROC)
library(yardstick)
library(rpart)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
#Lettura dataset originale non normalizzato
data <- read.csv("crimes.csv")
glimpse(data)
apply(data,2,function(x) sum(is.na(x)))
#sono presenti molti NA

# Preprocessing -----------------------------------------------------------

data2 <- data[,c(-3,-4,-(104:120), -(124:129))] 
#tolte quelle relative alla polizia e a county
data2 <- data2[,-c(14, 24, 31, 45, 53, 55, 83, 85, 86, 87, 89, 90)]
#tolte quelle ridondanti
data2 <- data2[data2$state!="AK",]
#tolta alaska
data2["2006","OtherPerCap"]=mean(na.omit(data2$OtherPerCap))
#unico NA al di fuori delle ultime variabili
#lo si riempie con la media delle altre

# KNN ---------------------------------------------------------------------
#capire se il knn da una buona accuracy

data_mascherato <- na.omit(data2)
sum(is.na(data_mascherato$rapesPerPop))
length(unique(data_mascherato$state))
states <- unique(data_mascherato$state)
best_K <- rep(NA,43)
best_RMSE <- rep(NA,43)
'for(j in 1:43){
  y <- data_mascherato[,c("rapes","rapesPerPop")]
  y[data_mascherato$state==states[j],"rapes"] <- NA
  y[data_mascherato$state==states[j],"rapesPerPop"] <- NA
  x <- data_mascherato[,-c(1,2,3,93:110)]
  rmse_arr <- rep(NA,100)
  for(i in 1:100){
    pred <- knn.reg(train = x[data_mascherato$state!=states[j],], test = x[data_mascherato$state==states[j],], y = y[data_mascherato$state!=states[j],"rapesPerPop"], k = i)$pred
    rmse <- sqrt(mean((data_mascherato[data_mascherato$state==states[j],"rapesPerPop"] - pred)^2))
    rmse_arr[i] = rmse
  }
  best_RMSE[j] <- rmse_arr[which.min(rmse_arr)]
  best_K[j] <- which.min(rmse_arr)
}
plot(best_K,best_RMSE)'
k <- 8

# Ripartizione in 4 aree geografiche -------------------------------------------

# legenda:
# 1 = Midwest
# 2 = South
# 3 = Northeast
# 4 = West

#midwest
rip<- ifelse(data2[, "state"] %in% list("IL","IN","IA", "KS", "MI", "MN",
                                        "MO", "NE", "ND", "OH", "SD", "WI"),1,4)
#south
rip <- ifelse(data2[, "state"] %in% list("AL","AR","DE", "FL", "GA", "KY",
                                         "LA", "MD", "MS", "NC", "OK", "SC", 
                                         "TN", "TX", "VA", "WV","DC"),2,rip)
#northeast
rip <- ifelse(data2[, "state"] %in% list("CT","ME","MA", "NH", "NJ", "NY",
                                         "PA", "RI", "VT"),3,rip)

data2 <- cbind(data2, rip) 
#ora nel dataset c'è una colonna che identifica la ripartizione
#geografica di appartenenza
midwest <- ifelse(data2$rip==1,1,0)
south <- ifelse(data2$rip==2,1,0)
northeast <- ifelse(data2$rip==3,1,0)
west <- ifelse(data2$rip==4,1,0)
#si creano 4 dummy da usare nel KNN
data2 <- cbind(data2,midwest,south,northeast,west)
crimes <- colnames(data2[,c(94,96,98,100,102,104,106,108)])
#riempimento degli NA con stime KNN
for(i in 1:8){
  x_train <- data2[!is.na(data2[,crimes[i]]),-c(1,2,3,93:111)]
  Y_train <- data2[!is.na(data2[,crimes[i]]),crimes[i]]
  x_test <- data2[is.na(data2[,crimes[i]]),-c(1,2,3,93:111)]
  if(dim(x_test)[1]!=0){
    kn <- knn.reg(x_train,x_test,Y_train,k)$pred 
    data2[is.na(data2[,crimes[i]]),crimes[i]] <- round(kn,2) #riempiamo dove c'erano gli NA
  }
}
data2[,109] <- apply(data2[,c(94,96,98,100)],1,sum)
data2[,110] <- apply(data2[,c(102,104,106,108)],1,sum)
data2 <- data2[,-c(93,95,97,99,101,103,105,107)]
#tolti i valori assoluti dei crimini
sum(is.na(data2)) #non ci sono più NA
write.csv(data2,"DataNoNa.csv",row.names=F)
