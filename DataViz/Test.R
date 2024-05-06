install.packages("igraph") 
install.packages("network") 
install.packages("sna")
install.packages("ggraph")
install.packages("visNetwork")
install.packages("threejs")
install.packages("networkD3")
install.packages("ndtv")

library(ndtv)
library(network)
library(networkDynamic)
library(animation)

l_degree <- data.frame()
num_id <- unique(rDataU$agentid)+1
for(i in 750){
  d <- data.frame(id = num_id, degree = degree(network.extract(net, at=i)),
                  U = rDataU$utility[rDataU$round == i], round = i)
  l_degree <- rbind(l_degree,d)
}

reg <- lm(utility~degree, data = dfF)

ggplot(data = dfF[dfF$round == max(dfF$round),]) +
  geom_point(aes(x = degree, y = utility), position = 'jitter') +
  geom_abline(slope = reg$coefficients[2], intercept = reg$coefficients[1], col = 'red')
