# Set the path to the Outputs folder
outputs_path <- "~/Documents/SMS/DataViz/NewOutputs"

# Get a list of sub-folders (A, B, C, D) in the Outputs folder
subfolders <- list.files(outputs_path, full.names = TRUE)
subfolders <- subfolders[c(4,1,2,3)]
# Create empty lists to store data frames
l1 <- list()
l2 <- list()
l3 <- list()
l4 <- list()

# Loop through each sub-folder
for (s in 1:4) {
  subfolder = subfolders[s]
  
  print(s)
  # Get the list of files matching the pattern
  matching_files <- list.files(subfolder, pattern = "SMS-Save-V", full.names = TRUE)
  
  # Loop through each matching file
  for (i in seq_along(matching_files)) {
    # Read the data frame from the file
    df <- read.csv(matching_files[i], header=FALSE)  # You may need to adjust this based on your actual file format
    df <- df[-c(1),]
    names(df) <- c("round","agentid","gamma","isgreedy", "meandist", "vardist", "maxdist", "P")
    # Assign the data frame to the corresponding list and name
    eval(parse(text = paste0("l", s, "[[i]] = df")))
  }
}

# Access the data frames using l1, l2, l3, l4
# For example:

for(s in 1:4){
  cmd = paste0("meandist_", s, " <- numeric()")
  eval(parse(text = cmd))
  
  cmd = paste0("for(i in 1:length(l",s,")){meandist_",s," <- append(meandist_", s, ", l", s, "[[i]]$meandist)}")
  eval(parse(text = cmd))
  
  cmd = paste0("meandist_",s," = as.numeric(meandist_",s,")")
  eval(parse(text = cmd))
  
  cmd = paste0("meandist_",s," = meandist_",s,"[meandist_",s,"> 0.1]")
  eval(parse(text = cmd))
}

d1 <- data.frame(x = 50, y = meandist_1)
d2 <- data.frame(x = 100, y = meandist_2)
d3 <- data.frame(x = 250, y = meandist_3)
d4 <- data.frame(x = 500, y = meandist_4)
d <- rbind(d1,d2,d3,d4)

dstat = data.frame(x=c(50,100,250,500), y=c(mean(meandist_1), mean(meandist_2), mean(meandist_3), mean(meandist_4)))
reg = lm(y ~ log(x), dstat)
xx = seq(10, 550, 1)
yy = predict(reg, data.frame(x=xx))
  
ggplot(aes(x=x,y=y), data = d) +
  geom_point(alpha = 0.1) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "red") +
  geom_line(aes(x=x,y=y), data = data.frame(x=xx,y=yy)) +
  labs(x = 'Number of agent', y = "Mean distance") +
  theme_light()



#############################################

# Set the path to the Outputs folder
outputs_path <- "/Users/yannkerzreho/Documents/SMS/DataViz/NewOutputs"

# Get a list of sub-folders (A, B, C, D) in the Outputs folder
subfolders <- list.files(outputs_path, full.names = TRUE)
subfolders <- subfolders[c(4,1,2,3)]
# Create empty lists to store data frames
lU1 <- list()
lU2 <- list()
lU3 <- list()
lU4 <- list()

for (s in 1:4) {
  subfolder = subfolders[s]
  
  print(s)
  # Get the list of files matching the pattern
  matching_files <- list.files(subfolder, pattern = "SMS-Save-U", full.names = TRUE)
  
  # Loop through each matching file
  for (i in seq_along(matching_files)) {
    # Read the data frame from the file
    df <- read.csv(matching_files[i], header=T)  
    df$size = s
    # Assign the data frame to the corresponding list and name
    eval(parse(text = paste0("lU", s, "[[i]] = df")))
  }
}

DU <- data.frame()
for(s in 1:4) {
  cmd = paste0("for(i in 1:length(lU",s,")){DU <- rbind(DU,lU", s, "[[i]])}")
  eval(parse(text = cmd))
}

DU <- subset(DU, (round == 200 & size == 1) |
                        (round == 500 & size == 2) |
                        (round == 750 & size == 3) |
                        (round == 1000 & size == 4))


#######################################################################################################################################
library(stringr)
list_data <- list()
subfolder = subfolders[1]
# Get the list of files matching the pattern
matching_files_U <- list.files(subfolder, pattern = "SMS-Save-U", full.names = TRUE)
matching_files_V <- list.files(subfolder, pattern = "SMS-Save-V", full.names = TRUE)
matching_files_E <- list.files(subfolder, pattern = "SMS-Save-E", full.names = TRUE)

for(i in seq_along(matching_files_U)) {
  print(paste0("Treating the ", i, " data sets..."))
  # Read the data frame from the file
  dfU <- read.csv(matching_files_U[i], header=T)  
  dfU$agentid <- dfU$agentid + 1
  dfV <- read.csv(matching_files_V[i], header=T, numerals = 'no.loss') #No loss allow to keep the exact P
  dfV$agentid <- dfV$agentid + 1
  dfV$mp <- str_count(as.character(dfV$P), "1")/100
  dfE <- read.csv(matching_files_E[i], header=T)
  
  Data <- prepareData(dfE)
  
  print(paste0("Doing the graph..."))
  net <- networkDynamic(edge.spells=Data, 
                        create.TEAs = TRUE,
                        edge.TEA.names = c('w'))
  
  num_id <- unique(dfU$agentid) + 1
  dfD <- data.frame(id = num_id, degree = degree(network.extract(net, at=1000)))
  
  print(paste0("Log..."))
  eval(parse(text = paste0("list_data[[", i, "]] = list(dfU=dfU, dfV=dfV, Data=Data, net=net, dfD=dfD, dfE=dfE)")))
}

treatedDF <- data.frame()
for(i in 1:length(list_data)){ #length(list_data)
  print(i)
  du <- list_data[[i]]$dfU 
  du <- du[du$round == 999,]
  dv <- list_data[[i]]$dfV 
  dv <- dv[dv$round == 0,]
  dd <- list_data[[i]]$dfD
  
  du <- du[order(du$agentid),]
  dv <- dv[order(dv$agentid),]
  dd <- dd[order(dd$id),]
  
  data <- data.frame("Data_set" = i, "agentid" = du$agentid, "Utility" = du$utility,
                     "Gamma" = dv$gamma, "Mean_Dist" = dv$meandist, "Mean_P" = dv$mp,
                     "Degree" = dd$degree)
  treatedDF <- rbind(treatedDF, data)
}

corrplot::corrplot.mixed(cor(data[-c(1,2)]), order = 'AOE')

heatmap(cor(data[-c(1)]))

corrplot::corrplot(cor(data[-c(1)]))

GGally::ggpairs(data[-c(1,2)])


rd <- list_data[[1]]$dfE
W <- rd$old_weight - rd$new_weight
d <-  data.frame(w=W,r=rd$round,a=rd$accepted)
d <- d[!d$r == 0,]
ggplot() +
  geom_point(aes(x=r,y=w), data = d[d$a == 0,], alpha = 0.2) +
  geom_point(aes(x=r,y=w), col = 'red', data = d[d$a == 1,], alpha = 0.2)
