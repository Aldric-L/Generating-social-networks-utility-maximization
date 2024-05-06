# doc : 
# https://programminghistorian.org/en/lessons/temporal-network-analysis-with-r
# https://kateto.net/network-visualization


library(ndtv)
library(network)
library(networkDynamic)
library(animation)
library(ggplot2)

rData <- read.csv("~/Documents/SMS/DataViz/SMS-Save-Edges-235184180.csv")
rDataU <- read.csv("~/Documents/SMS/DataViz/SMS-Save-Utility-235184180.csv")
rDataV <- read.csv("~/Documents/SMS/DataViz/SMS-Save-Vertices-235184180.csv", header=FALSE)
rDataV <- rDataV[-c(1),]
names(rDataV) <- c("round","agentid","gamma","isgreedy", "meandist", "vardist", "maxdist", "P")
rDataV <- rDataV[rDataV$round == 0,]
rDataV$agentid <- as.numeric(rDataV$agentid) + 1
list_P <- list()
for(i in 1:nrow(rDataV)){
  list_P[[i]] <- na.omit(as.numeric(strsplit(rDataV$P[i], "")[[1]]))
}
rDataV$mp <- NA
for(i in 1:nrow(rDataV)){
  rDataV$mp[i] <- mean(list_P[[i]])
}

rData <- rData #[rData$accepted == 1,]
W <- rData$old_weight - rData$new_weight
d <-  data.frame(w=W,r=rData$round,a=rData$accepted)
d <- d[!d$r == 0,]
ggplot() +
  geom_point(aes(x=r,y=w), data = d[d$a == 0,], alpha = 0.2) +
  geom_point(aes(x=r,y=w), col = 'red', data = d[d$a == 1,], alpha = 0.2)


Data <- prepareData(rData)

net <- networkDynamic(edge.spells=Data, 
                      create.TEAs = TRUE,
                      edge.TEA.names = c('w'))


##
Data <- list_data[[1]]$Data
net <- list_data[[1]]$net
dfU <- list_data[[1]]$dfU
dfV <- list_data[[1]]$dfV
dfV <- dfV[dfV$round == unique(dfV$round[2]),]
dfD <- list_data[[1]]$dfD
##

list_P <- list()
for(i in 1:nrow(dfV)){
  list_P[[i]] <- na.omit(as.numeric(strsplit(dfV$P[i], "")[[1]]))
}

l_degree <- data.frame()
num_id <- rDataV$agentid
for(i in 1:max(Data$onset)){
  d <- data.frame(id = num_id, degree = degree(network.extract(net, at=i)), round = i)
  l_degree <- rbind(l_degree,d)
}
min <- min(l_degree$degree)
max <- max(l_degree$degree)

wf <- Data[Data$terminus == 'Inf',]$w 
plot(density(wf))
sum(wf<0.01)/length(wf)

vertex_size_function <- function(slice) {
  degrees <- (degree(slice) - min)/(max - min) + 0.5
  degrees
}
library(viridis)
color_palette <- viridis(1000, option = 'B')

edge_color_function <- function(slice) {
  color_index <- cut((slice %e% "w")^(1/2), breaks = length(color_palette), include.lowest = TRUE)
  color_palette[color_index]
}

## Check negative U
neg <- rDataU[rDataU$round == 999 & rDataU$utility < 0, ]$agentid + 1 
col_neg <- ifelse(1:100 %in% neg, 0, 1)
###

## Check affinity
Data[Data$terminus == Inf,]

indiv = 12
affinity <- as.numeric(lapply(list_P, function(x){sum(x*list_P[[indiv]])/100}))
affinity[indiv] <- 0.2
af <- (affinity - min(affinity)) / (max(affinity) - min(affinity)) + 0.001
color_palette_aff <- viridis(1000, option = 'B')
color <- color_palette_aff[floor(af*1000)]
color[indiv] <- 0
##

compressed_data <- compress_netdata(Data,10)
compressed_net <- networkDynamic(edge.spells=compressed_data, 
                      create.TEAs = TRUE,
                      edge.TEA.names = c('w'))

render.d3movie(compressed_net, usearrows = F, displaylabels = F, 
               render.par=list(tween.frames=2),
               edge.col =  edge_color_function,
               edge.lwd = function(slice)(slice %e% "w" + 0.3),
               vertex.cex = vertex_size_function,
               vertex.col = color,
               output.mode=c('HTML'),script.type='embedded',filename='cc.html')

## Check
d <- Data[(Data$head == indiv | Data$tail == indiv) & Data$terminus == Inf,]
d$aff <- NA
for(i in 1:nrow(d)){d$aff[i] <- sum(list_P[[d$tail[i]]]*list_P[[d$head[i]]])/100}
print(d)
