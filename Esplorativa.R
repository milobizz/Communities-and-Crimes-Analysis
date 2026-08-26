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
library(patchwork)
data <- read.csv("DataNoNa.csv")


# Esplorativa univariata --------------------------------------------------
hist(data$ViolentCrimesPerPop,nclass=100)
# categoria 1 : variabili socio demografiche -------------------------------------------------------------

#composizione demografica 
comp_demografica_usa <- function(data) {
  comp_pop <- rep(NA, 4)
  for (i in 6:9) {
    comp_pop[i-5] = sum(data[,i]*data$population/100) #trovo numero di persone di quella razza
  }
  tot <- sum(data[,6:9]*data$population/100)
  rbind(round(comp_pop/tot, 2),colnames(data[,6:9]))
}

comp_demografica_usa(data)
#il 68% è di razza bianca, il 16% nero, il 12% ispanico e il 4% asiatico

ggpairs(data[,10:13], aes(fill="orange", alpha=0.5))
#ci sono variabili correlate perché fasce di età 12-21, 12-29 e 16-24 si sovrappongono
#si opta per tenere solo 12-29 


# categoria 2 - variabili socio-economiche-------------------------------------------------------------------

#distribuzione del reddito pro capite a livello globale nelle varie etnie
df_long <- reshape2::melt(data[, c("whitePerCap", "blackPerCap",
                                   "indianPerCap", "AsianPerCap",
                                   "OtherPerCap","HispPerCap")],
                          variable.name = "race",
                          value.name = "income_per_capita")
ggplot(df_long) +
  geom_boxplot(aes(x="race",y = income_per_capita, fill=race))+
  coord_flip(ylim = c(0, 1e+05))+ #messo limite per zoomare sulle distribuzioni, c'erano outliers che schiacciavano grafico
  scale_fill_brewer(palette = "Set2")

tot_pop_usa <- sum(data$population)
tot_pop_usa
tot_pctUrban <- sum(data[,"pctUrban"]/100*data$population) / tot_pop_usa
tot_pctUrban
tot_pctPoverty <- sum(data[,"PctPopUnderPov"]/100*data$population) / tot_pop_usa
tot_pctPoverty

by(data[,"PctPopUnderPov"], data$rip, mean)

#reddito medio all interno delle aree geografiche
by(data$medIncome, data$rip, mean)
#gap tra etnia caucasica e afro-americana costante in tutte le aree
by(data[,"whitePerCap"]-data[,"blackPerCap"], data$rip, mean)

cor(data[,16:21])
cor.test(data[,16], data[,19])
#pctWWage e pctWSocSec correlazione significativa = -0.90 (p-value < 2.2e-16)

by(data$medIncome, data$rip, mean) #stati northeast hanno reddito mediano più alto

by(data$PopDens, data$rip, mean)

# categoria 3 - istruzione e lavoro -------------------------------------------------------------

by(data$PctLess9thGrade*data$population/100, data$rip, sum)/tot_pop_usa
by(data$PctNotHSGrad*data$population/100, data$rip, sum)/tot_pop_usa
by(data$PctBSorMore*data$population/100, data$rip, sum)/tot_pop_usa

sum(data$PctUnemployed*data$population/100)/tot_pop_usa
by(data$PctUnemployed*data$population/100, data$rip, sum)/tot_pop_usa

sum(data$PctEmploy*data$population/100)/tot_pop_usa

cor(data[,33:38]) #variabili correlate tra loro



# categoria 4 - famiglia, divorzio ----------------------------------------
cor(data[,39:49]) 
p4 = data[, 39:49]
GGally::ggpairs(p4, aes(color = as.factor(data$rip)))

by(p4$PersPerFam, data$rip, summary)
by(p4$PersPerFam, data$rip, sd)

#plot mappa US 

#divorzi maschili
divmale = by(data$MalePctDivorce*data$population/100, data$state, sum)/tot_pop_usa
state_data = data.frame(state = names(divmale), divmale = divmale)

plot_usmap(data = state_data, values = "divmale", regions = "states",exclude = c("HI","AK")) +
  scale_fill_continuous(low = "white", high = "red", name = "Divmaschi") +
  theme(legend.position = "right")

#divorzi femminili
divfem = by(data$FemalePctDiv*data$population/100, data$state, sum)/tot_pop_usa
state_data = data.frame(state_data, divfem = divfem)

plot_usmap(data = state_data, values = "divfem", regions = "states",exclude = c("HI","AK")) +
  scale_fill_continuous(low = "white", high = "red", name = "Divdonne") +
  theme(legend.position = "right")

#figli con 2 genitori
twopar = by(data$PctKids2Par*data$population/100, data$state, sum)/tot_pop_usa
state_data = data.frame(state_data, twopar = twopar)

plot_usmap(data = state_data, values = "twopar", regions = "states",exclude = c("HI","AK")) +
  scale_fill_continuous(low = "white", high = "red", name = "Figli 2 genitori") +
  theme(legend.position = "right")

#figli con genitori non sposati
nomar = by(data$PctKidsBornNeverMar*data$population/100, data$state, sum)/tot_pop_usa
state_data = data.frame(state_data, nomar = nomar)

plot_usmap(data = state_data, values = "nomar", regions = "states",exclude = c("HI","AK")) +
  scale_fill_continuous(low = "white", high = "red", name = "Figli genitori\n non sposati") +
  theme(legend.position = "right")

# categoria 5 - Immigrazione e Lingua ------------------------------------

# Grafico per la variabile più rappresentativa del gruppo popolazione
# ovvero la % della popolazione immigrata negli ultimi 10 anni:

# Popolazione totale degli USA
tot_pop_usa = sum(data$population)

# Calcolo per ogni stato il numero assoluto di immigrati recenti.
# (data$PctRecImmig10*data$population) = numero di immigrati recenti
# in quella comunità.
# Sommiamo i valori raggruppando per stato del numero totale di immigrati
# recenti in una determinata comunità: 
c<-by(data$PctRecImmig10*data$population, data$state, sum)/100


quota_stato = data.frame(
  state = names(c),
  quota_media = as.numeric(c) / tot_pop_usa
)

# Mappa dell'immigrazione negli ultimi 10 anni
plot_usmap(data = quota_stato, values = "quota_media", regions = "states",exclude = c("HI","AK")) +
  scale_fill_gradient(
    low  = "white",
    high = "red",
    name = "Quota immigrati\nsulla popolazione"
  )+theme(legend.position = "right")

# categoria 6 - costo della vita, qualità e condizioni abitative -------------------------------------------------------------
cor(data[,60:66]) #alta correlazione

cor(data[,67:77])

cor(data$NumStreet,data$NumInShelters)
cor.test(data$NumStreet,data$NumInShelters)

sum(data$PctLargHouseOccup*data$population/100)/tot_pop_usa

by(data$PctForeignBorn*data$population/100, data$rip, sum)/tot_pop_usa
sum(data$PctForeignBorn*data$population/100) / tot_pop_usa

cor(data[,85:89])

# categoria 7 - crimini --------------------------------------------------
medie_crimini <- apply(data[,93:102],2,mean)
Violent <- data.frame("Violent"=medie_crimini[1:4])
NonViolent <- data.frame("NonViolent"=medie_crimini[5:8])
# Violent
Violent_long <- data.frame(
  tipo  = rownames(Violent),
  media = Violent$Violent
)

p1 <- ggplot(Violent_long, aes(fill=tipo,x = 1, y = media)) +
  geom_bar(stat = "identity",position="fill") +
  labs(title = "Violent Crimes", x = "",y="")+
  scale_fill_brewer(palette = "Spectral")

# NonViolent
NonViolent_long <- data.frame(
  tipo  = rownames(NonViolent),
  media = NonViolent$NonViolent
)

p2 <- ggplot(NonViolent_long, aes(fill=tipo,x = 1, y = media)) +
  geom_bar(stat = "identity",position="fill") +
  labs(title = "Non-Violent Crimes", x = "",y="")+
  scale_fill_brewer(palette = "Spectral")

p1+p2
sum(ifelse(data$murdPerPop==0,1,0))
sum(ifelse(data$rapesPerPop==0,1,0))
sum(ifelse(data$robbbPerPop==0,1,0))
sum(ifelse(data$assaultPerPop==0,1,0))
sum(ifelse(data$burglPerPop==0,1,0))
sum(ifelse(data$larcPerPop==0,1,0))
sum(ifelse(data$autoTheftPerPop==0,1,0))
sum(ifelse(data$arsonsPerPop==0,1,0))
#molti 0 nei crimini violenti, giustifica certe bande così piccole
#323 0 negli arsons, giustifica la banda così piccola

# Esplorativa bivariata ---------------------------------------------------

#correlazioni tra il tasso di criminalità (ViolentCrimesPerPop) e variabili
#selezionate da 5 macroaree: reddito, livello di istruzione, immigrazione, caratteristiche familiari, povertà.

#reddito 
cor.test(data$ViolentCrimesPerPop, data$medIncome, method="spearman")
ggplot(data)+
  geom_point(aes(medIncome,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

#istruzione
cor.test(data$ViolentCrimesPerPop, data$PctLess9thGrade, method="spearman")
ggplot(data)+
  geom_point(aes(PctLess9thGrade,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

cor.test(data$ViolentCrimesPerPop, data$PctNotHSGrad, method="spearman")
ggplot(data)+
  geom_point(aes(PctNotHSGrad,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

#immigrazione
cor(data$ViolentCrimesPerPop,data[,50:59])

#famiglia
cor(data$ViolentCrimesPerPop,data[,39:49])

cor.test(data$ViolentCrimesPerPop, data$PctKids2Par)
ggplot(data)+
  geom_point(aes(PctKids2Par,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

cor.test(data$ViolentCrimesPerPop, data$PctFam2Par)
ggplot(data)+
  geom_point(aes(PctFam2Par,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

cor.test(data$ViolentCrimesPerPop, data$PctKidsBornNeverMar)
ggplot(data)+
  geom_point(aes(PctKidsBornNeverMar,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

#povertà
cor.test(data$ViolentCrimesPerPop, data$pctWPubAsst)
ggplot(data)+
  geom_point(aes(pctWPubAsst,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

cor.test(data$ViolentCrimesPerPop, data$PctPopUnderPov)
ggplot(data)+
  geom_point(aes(PctPopUnderPov,ViolentCrimesPerPop,colour = factor(rip)))+
  labs(
    colour = "ripartizione\n geografica",
  )+
  scale_colour_discrete(labels = c("1" = "mid-west", "2" = "south", "3" = "north-east", "4" = "west"))

