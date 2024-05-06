library(viridis)
library(shiny)
library(ggplot2)

list_a <- seq(0,1,0.005)

U <- function(a,G,gamma, delta){
  10*sum(a*G) - sum(a^gamma/(1-a^gamma)) - sum(a)^delta
}

Best_response_g1g2 <- function(list_a, gamma, delta, g1, g2) { #Best response given g1 and g2
  G <- c(g1,g2)
  grid <- as.data.frame(expand.grid(list_a,list_a))
  names(grid)[1:2] <- c("a1","a2")
  grid$U <- NA
  for(i in 1:nrow(grid)){
    grid$U[i] <- U(grid[i,c(1,2)],G,gamma,delta)
  }
  grid$U[grid$U<0] <- 0
  max_row <- grid[which.max(grid$U), ]
  max_row
}

Best_response_g1 <- function(list_a, gamma, delta, g1) { #Best response given g1
  g2 <- seq(0.1,0.5,0.01)
  a1 <- numeric(length(g2))
  a2 <- numeric(length(g2))
  U <- numeric(length(g2))
  for(it in 1:length(g2)){
    print(g2[it])
    r <- Best_response_g1g2(list_a, gamma, delta, g1, g2[it])
    a1[it] <- r$a1
    a2[it] <- r$a2
    U[it] <- r$U
  }
 list(a1 = a1, a2 = a2, U = U, g1 = g1, g2 = g2)
}

b2 <- Best_response_g1(list_a, 10, 2, 0.25)
b1 <- Best_response_g1(list_a, 10, 1, 0.25)
b4 <- Best_response_g1(list_a, 10, 4, 0.25)


data_grid2 = data.frame(a1=b2$a1, a2=b2$a2, U=b2$U, g2=b2$g2)
data_grid1 = data.frame(a1=b1$a1, a2=b1$a2, U=b1$U, g2=b1$g2)
data_grid4 = data.frame(a1=b4$a1, a2=b4$a2, U=b4$U, g2=b4$g2)

data_text = data.frame(x=c(0.04,0.43, 0.8), y=0.9,label = c("d=4","d=2","d=1"))

ggplot() +
  xlim(c(0,1)) +
  ylim(c(0,1)) +
  geom_line(aes(x = a1, y = a2, color = g2), data = data_grid2) +
  geom_line(aes(x = a1, y = a2, color = g2), data = data_grid1) +
  geom_line(aes(x = a1, y = a2, color = g2), data = data_grid4) +
  scale_color_viridis(option = "turbo") +
  geom_abline(intercept = 0, slope = 1, lty = 2) +
  theme_light() +
  geom_text(aes(x=x,y=y,label=label), data = data_text) +
  labs(color = "Affinity with N2", x = "weight with N1", y = "weight with N2")


##

U_conv <- dfU$utility[dfU$round == 999]
mean_p <- dfV$mp
nb_degree <- dfD$degree
ggplot()+
  geom_point(aes(mean_p,U_conv,col=log(nb_degree+1))) +
  scale_color_gradientn(colours = viridis(100))

data.frame(id = dfV$agentid, nb_degree = nb_degree, U_conv = U_conv)

p <- "ggplot()"
for(i in 1:nrow(dfU[sample(unique(dfU$agentid), 100),])){
  p <- paste0(p, "+ geom_line(aes(x=round, y=utility), col =", i,", data = dfU[dfU$agentid == ", i, "-1,])")
}
labs(color = "Affinity with N2", x = "weight with N1", y = "weight with N2")


library(dplyr)

meanU_round <- dfU %>%
  group_by(round) %>%
  summarise(mean_value = sum(utility))

eval(parse(text = p))

