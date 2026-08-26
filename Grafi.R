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
data <- read.csv("DataNoNa.csv")

# GRAFO -------------------------------------------------------------------
states <- matrix(unique(data$state))
#matrice di adiacenze tra stati
mat <- matrix(NA, nrow = 47, ncol = 9)
mat[1, 1:5] <- c("AL","MS","TN","GA","FL")
mat[2, 1:7] <- c("AR","MO","TN","MS","LA","TX","OK")
mat[3, 1:6] <- c("AZ","CA","NV","UT","CO","NM")
mat[4, 1:4] <- c("CA","OR","NV","AZ")
mat[5, 1:8] <- c("CO","WY","NE","KS","OK","NM","AZ","UT")
mat[6, 1:4] <- c("CT","NY","MA","RI")
mat[7, 1:3] <- c("DC","MD","VA")
mat[8, 1:4] <- c("DE","MD","PA","NJ")
mat[9, 1:3] <- c("FL","AL","GA")
mat[10, 1:6] <- c("GA","FL","AL","TN","NC","SC")
mat[11, 1:7] <- c("IA","MN","WI","IL","MO","NE","SD")
mat[12, 1:7] <- c("ID","MT","WY","UT","NV","OR","WA")
mat[13, 1:6] <- c("IL","IN","KY","MO","IA","WI")
mat[14, 1:5] <- c("IN","MI","OH","KY","IL")
mat[15, 1:5] <- c("KS","NE","MO","OK","CO")
mat[16, 1:8] <- c("KY","IN","OH","WV","VA","TN","MO","IL")
mat[17, 1:4] <- c("LA","TX","AR","MS")
mat[18, 1:6] <- c("MA","RI","CT","NY","NH","VT")
mat[19, 1:6] <- c("MD","VA","WV","PA","DC","DE")
mat[20, 1:2] <- c("ME","NH")
mat[21, 1:4] <- c("MI","WI","IN","OH")
mat[22, 1:5] <- c("MN","WI","IA","SD","ND")
mat[23, 1:9] <- c("MO","IA","IL","KY","TN","AR","OK","KS","NE")
mat[24, 1:5] <- c("MS","LA","AR","TN","AL")
mat[25, 1:5] <- c("NC","VA","TN","GA","SC")
mat[26, 1:4] <- c("ND","MN","SD","MT")
mat[27, 1:4] <- c("NH","VT","ME","MA")
mat[28, 1:4] <- c("NJ","DE","PA","NY")
mat[29, 1:6] <- c("NM","AZ","UT","CO","OK","TX")
mat[30, 1:6] <- c("NV","ID","UT","AZ","CA","OR")
mat[31, 1:6] <- c("NY","NJ","PA","VT","MA","CT")
mat[32, 1:6] <- c("OH","PA","WV","KY","IN","MI")
mat[33, 1:7] <- c("OK","KS","MO","AR","TX","NM","CO")
mat[34, 1:5] <- c("OR","CA","NV","ID","WA")
mat[35, 1:7] <- c("PA","NY","NJ","DE","MD","WV","OH")
mat[36, 1:3] <- c("RI","CT","MA")
mat[37, 1:3] <- c("SC","GA","NC")
mat[38, 1:7] <- c("SD","ND","MN","IA","NE","WY","MT")
mat[39, 1:9] <- c("TN","KY","VA","NC","GA","AL","MS","AR","MO")
mat[40, 1:5] <- c("TX","NM","OK","AR","LA")
mat[41, 1:7] <- c("UT","ID","WY","CO","NM","AZ","NV")
mat[42, 1:7] <- c("VA","NC","TN","KY","WV","MD","DC")
mat[43, 1:4] <- c("VT","NY","NH","MA")
mat[44, 1:3] <- c("WA","ID","OR")
mat[45, 1:5] <- c("WI","MI","MN","IA","IL")
mat[46, 1:6] <- c("WV","OH","PA","MD","VA","KY")
mat[47, 1:7] <- c("WY","MT","SD","NE","CO","UT","ID")
states=sort(states)
adj <- matrix(0,47,47)
colnames(adj) <- states
rownames(adj) <- states
#matrice di adiacenza con 0 e 1
for(i in 1:47){
  for(j in 1:47){
    if(states[j] %in% mat[i,-1]){
      adj[i,j] = 1
    }
  }
}
#matrice di distanza euclidea tra stati
data_states <- by(data[,-c(1,2,3)],data$state,colMeans)
line <- c(names(data_states)[1],data_states[[1]])
dist <- line
for(i in 2:47){
  line <- c(names(data_states)[i],data_states[[i]])
  dist <- rbind(dist,line)
}
dist <- data.frame(dist)
rownames(dist) <- dist$V1
dist <- dist[,-c(1,101:106)]
dist <- apply(dist,2,as.numeric)
dist <- scale(dist)
dist <- as.matrix(dist(dist))
rownames(dist) <- sort(unique(data$state))
colnames(dist) <- sort(unique(data$state))
dist <- dist[sort(rownames(dist)), sort(colnames(dist))]
dist

#grafo con dimensione dei nodi in base al tasso di crimini violenti
#e colore degli archi in base alla distanza euclidea
tassi <- data[,c(2,101)]
tassi[,2]<-as.numeric(tassi[,2])
g <- graph_from_adjacency_matrix(adj*dist, mode = "undirected", diag = FALSE,weighted = T)
#posizioni nello spazio degli stati
coords_df <- data.frame(state.center, row.names = state.abb)
coords_df["DC", ] <- c(-77.01, 38.91)
coords_df["NV", ] <- c(-115, 40)
#si spostano leggermente alcuni stati per vederli meglio
lay <- as.matrix(coords_df[V(g)$name, c("x", "y")])
tassigr <- by(tassi$ViolentCrimesPerPop,tassi$state,sum)
tassi2 <- by(data$ViolentCrimesPerPop,data$state,sum)
size <- 150+500*(tassigr-min(tassi2))/(max(tassi2)-min(tassi2))
state_region <- data.frame(
  state = c("alabama","alaska","arizona","arkansas","california","colorado","connecticut","district of columbia","delaware","florida","georgia","hawaii","iowa","idaho","illinois","indiana","kansas","kentucky","louisiana","massachusetts","maryland","maine","michigan","minnesota","missouri","mississippi","montana","north carolina","north dakota","nebraska","new hampshire","new jersey","new mexico","nevada","new york","ohio","oklahoma","oregon","pennsylvania","rhode island","south carolina","south dakota","tennessee","texas","utah","virginia","vermont","washington","wisconsin","west virginia","wyoming"),
  region = c("South","West","West","South","West","West","Northeast","South","South","South","South","West","Midwest","West","Midwest","Midwest","Midwest","South","South","Northeast","South","Northeast","Midwest","Midwest","Midwest","South","West","South","Midwest","Midwest","Northeast","Northeast","West","West","Northeast","Midwest","South","West","Northeast","Northeast","South","Midwest","South","South","West","South","Northeast","West","Midwest","South","West"),
  stringsAsFactors = FALSE
)
state_region <- state_region[-c(2,12),]
# Ordine esatto degli stati in map("state")
map_names  <- map("state", plot = FALSE)$names
base_names <- sub(":.*", "", map_names)  # rimuove ":main", ":north" ecc.
# Allinea fill_cols all'ordine della mappa
region_cols <- c(Northeast="#BCE4D8", Midwest="#E2E4E7", South="#C4D8F3", West="#DDF1D7")
state_region_named <- setNames(state_region$region, state_region$state)
fill_cols <- region_cols[state_region_named[base_names]]

# Plot
xlim <- range(lay[, 1], na.rm = TRUE) + c(-2, 2)
ylim <- range(lay[, 2], na.rm = TRUE) + c(-2, 2)
plot(NA, NA, xlim = xlim, ylim = ylim,
     xlab = "", ylab = "", axes = FALSE, asp = 1)
map("state", interior = TRUE, add = TRUE,
    fill = TRUE, lwd = 1, col = fill_cols)
# Archi
w        <- E(g)$weight
w_scaled <- round(1+99 * (w - min(w)) / (max(w) - min(w)))
pal      <- colorRampPalette(c("#A50021","#EAC0BD","#FFFFFF"))(100)
edge_cols <- pal[w_scaled]
plot(g,
     layout       = lay,
     rescale      = FALSE,
     add          = TRUE,
     vertex.label = V(g)$name,
     vertex.size  = size,
     edge.color   = edge_cols,
     edge.width   = 3,
     vertex.color = "white",
     vertex.shape="circle",
     vertex.label.color="black")


# stati del Nord-Est
states_zoom <- c("WV", "VA", "DC", "MD", "DE", "PA", "NJ",
                 "NY", "CT", "RI", "MA", "NH", "VT", "ME")

# sottografo solo con gli stati scelti
g_ne <- induced_subgraph(g, vids = V(g)[name %in% states_zoom])
# coordinate solo per questi stati
coords_ne <- matrix(c(
  -80.4549, 38.5976,  # WV
  -78.6569, 37.4316,  # VA
  -77.0369, 38.9072,  # DC
  -76.7908, 39.0639,  # MD
  -75.5071, 38.9108,  # DE
  -77.1945, 41.2033,  # PA
  -74.4057, 40.0583,  # NJ
  -74.9481, 42.9134,  # NY
  -72.7554, 41.6032,  # CT
  -71.4774, 41.5801,  # RI
  -70.5, 42.2589,  # MA
  -71.3, 43.8,  # NH
  -72.5778, 44.75,  # VT
  -69.3819, 45.3706   # ME
), ncol = 2, byrow = TRUE)
rownames(coords_ne) <- states_zoom
colnames(coords_ne) <- c("x", "y")
# riordino coordinate secondo l'ordine dei nodi nel sottografo
lay_ne <- coords_ne[V(g_ne)$name, ]
# tassi per dimensione nodo
#tassigr
v_ne <- tassigr[V(g_ne)$name]
size_ne <- 50 + 100 * (v_ne - min(v_ne, na.rm = TRUE)) /
  (max(v_ne, na.rm = TRUE) - min(v_ne, na.rm = TRUE))
# colori

xlim <- range(lay_ne[,1], na.rm = TRUE) + c(-1.5, 1.5)
ylim <- range(lay_ne[,2], na.rm = TRUE) + c(-1.5, 1.5)
plot(
  NA, NA,
  xlim = xlim,
  ylim = ylim,
  xlab = "", ylab = "",
  axes = FALSE,
  asp = 1
)
map("state", interior = TRUE, add = TRUE,
    fill = TRUE, lwd = 1, col = fill_cols)
plot(g_ne,
     layout       = lay_ne,
     rescale      = FALSE,
     add          = TRUE,
     vertex.label = V(g_ne)$name,
     vertex.size  = size_ne,
     edge.color   = edge_cols[c(
       23, 26, 28, 61, 68, 29, 24, 62, 82, 30,
       65, 83, 89, 25, 63, 27, 66, 64, 81, 90,
       67, 95, 99
     )],
     edge.width   = 3,
     vertex.color = "white",
     vertex.shape="circle",
     vertex.label.color="black")

#divisione binaria violenti o non violenti

xlim <- range(lay[, 1], na.rm = TRUE) + c(-2, 2)
ylim <- range(lay[, 2], na.rm = TRUE) + c(-2, 2)

plot(NA, NA, xlim = xlim, ylim = ylim,
     xlab = "", ylab = "", axes = FALSE, asp = 1)
map("state", interior = TRUE, add = TRUE,
    fill = TRUE, lwd = 1, col = fill_cols)

state_risk <- by(data$ViolentCrimesPerPop,data$state, mean)
state_col <- ifelse(state_risk<500,"#6ae866","#e86f66")
plot(g,
     layout       = lay,
     rescale      = FALSE,
     add          = TRUE,
     vertex.label = V(g)$name,
     vertex.size  = size,
     edge.color   = "black",
     edge.width   = 2,
     vertex.color = state_col,
     vertex.shape="circle",
     vertex.label.color="black"
     )

             