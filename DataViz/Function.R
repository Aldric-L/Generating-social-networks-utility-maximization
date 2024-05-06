prepareData <- function(rData){
  from_to <- data.frame(tail = apply(rData[,c(2,3)], 1, min), head = apply(rData[,c(2,3)], 1, max))
  rData$tail <- from_to$tail + 1
  rData$head <- from_to$head + 1
  rData <- rData[ !(rData$round == 0 & rData$new_weight == 0), ]
  rData <- rData[rData$accepted == 1,]
  max_id <- max(rData$tail,rData$head)
  list_possible_edge <- list()
  k <- 1
  for(i in 2:max_id){
    for(j in 1:(i-1)){
      list_possible_edge[[k]] <- c(j,i)
      k <- k+1
    }
  }
  
  Data <- data.frame()
  for(p in 1:length(list_possible_edge)){
    tail = list_possible_edge[[p]][1]
    head = list_possible_edge[[p]][2]
    d <- rData[rData$tail == tail & 
                 rData$head == head, ]
    nd <- nrow(d)
    res <- data.frame()
    if(nd > 1){
      d <- d[order(d$round),]
      for(i in 1:(nd-1)){
        if(d$new_weight[i]>0) res <- rbind(res, c(d$round[i], d$round[i+1],
                                                  tail, head, d$new_weight[i]))
      }
      if(d$new_weight[nd]>0) res <- rbind(res, c(d$round[nd], Inf,
                                                 tail, head, d$new_weight[nd]))
    }
    if(nd == 1 & d$new_weight[1] > 0) res <- rbind(res, c(d$round[1], Inf,
                                                          tail, head, d$new_weight[1]))
    if(nrow(res)>0) names(res) <- c("onset", "terminus", "tail", "head", "w")
    Data <- rbind(Data, res)
  }
  Data
}



compress_netdata <- function(data, by = 5){
  tmin <- min(data$onset) + 1
  tmax <- max(data$onset)
  tseq <- seq(tmin, tmax, by)
  res <- data.frame()
  for(i in 1:(length(tseq)-1)){
    d <- data[data$onset %in% seq(tseq[i], tseq[i+1]), ]
    d <- d[d$terminus > tseq[i+1], ]
    d$onset <- i
    d$terminus <- floor(d$terminus / by)
    res <- rbind(res, d)
  }
  res
}

C <- function(x,gamma){
  x^gamma/(1-x^gamma)
}



SMSanalyse <- function(outputs_path, compression = 0, graph = FALSE){
  
  ## Find the files
  files_utility <- list.files(outputs_path, pattern = "SMS-Save-Utility", full.names = TRUE)
  if(length(files_utility) != 1) {
    stop("You should have a SMS-Save-Utility file")
  }
  files_edges <- list.files(outputs_path, pattern = "SMS-Save-Edges", full.names = TRUE)
  if(length(files_edges) != 1) {
    stop("You should have a SMS-Save-Edges file")
  }
  files_cluster <- list.files(outputs_path, pattern = "SMS-Save-Clustering", full.names = TRUE)
  if(length(files_edges) != 1) {
    stop("You should have a SMS-Save-Clustering file")
  }
  files_vertices <- list.files(outputs_path, pattern = "SMS-Save-Vertices", full.names = TRUE)
  if(length(files_edges) != 1) {
    stop("You should have a SMS-Save-Vertices file")
  }
  
  ## Preprocess the data
  print("Preprocess the data")
  dfE <- read.csv(files_edges)
  data <- prepareData(dfE)
  dfU <- read.csv(files_utility)
  dfU$agentid <- dfU$agentid + 1
  dfC <- read.csv(files_cluster)
  dfV <- read.csv(files_vertices, header=FALSE)
  dfV <- dfV[-c(1),]
  names(dfV) <- c("round","agentid","gamma","isgreedy", "meandist", "vardist", "maxdist", "P")
  dfV$agentid <- as.numeric(dfV$agentid) + 1
  dfV <- dfV[dfV$round == 0,]
  list_P <- list()
  for(i in 1:nrow(dfV)){
    list_P[[i]] <- na.omit(as.numeric(strsplit(dfV$P[i], "")[[1]]))
  }
  dfV$mp <- NA
  for(i in 1:nrow(dfV)){
    dfV$mp[i] <- mean(list_P[[i]])
  }
  
  if(compression > 0){
    data <- compress_netdata(data, compression)
  }
  
  ## Creating the net and analyse it
  print("Creating the Net")
  net <- networkDynamic::networkDynamic(edge.spells=data, 
                                        create.TEAs = TRUE,
                                        edge.TEA.names = c('w'))
  
  dfD <- data.frame()
  num_id <- dfV$agentid
  for(i in 1:max(Data$onset)){
    d <- data.frame(agentid = num_id, degree = degree(network.extract(net, at=i)), round = i)
    dfD <- rbind(dfD,d)
  }
  
  dfU <- merge(dfD, dfU, by = c("agentid", "round"))
  dfU <- merge(dfU, dfV[,c('agentid', 'gamma', 'mp')])
  
  if(graph == TRUE){
    ## Creating the render
    print("Creating the render")
    min_degree <- min(dfD$degree)
    max_degree <- max(dfD$degree)
    vertex_size_function <- function(slice) {
      degrees <- (degree(slice) - min)/(max - min) + 0.5
      degrees
    }
    color_palette <- viridis::viridis(1000, option = 'B')
    edge_color_function <- function(slice) {
      color_index <- cut((slice %e% "w")^(1/2), breaks = length(color_palette), include.lowest = TRUE)
      color_palette[color_index]
    }
    
    ## Check affinity
    
    indiv = sample(dfU[(dfU$round == max(dfU$round)) & (dfU$degree > 5),]$agentid, 1)
    affinity <- as.numeric(lapply(list_P, function(x){sum(x*list_P[[indiv]])/100}))
    affinity[indiv] <- 0.2
    af <- (affinity - min(affinity)) / (max(affinity) - min(affinity)) + 0.001
    color_palette_aff <- viridis(1000, option = 'B')
    color <- color_palette_aff[floor(af*1000)]
    color[indiv] <- 0
    ##
    
    render.d3movie(net, usearrows = F, displaylabels = F, 
                   render.par=list(tween.frames=2),
                   edge.col =  edge_color_function,
                   edge.lwd = function(slice)(slice %e% "w" + 0.3),
                   vertex.cex = vertex_size_function,
                   vertex.col = color,
                   output.mode=c('HTML'),script.type='embedded',filename='graph.html')
    
  }
  list(dfU=dfU, net=net)
}

outputs_path = "/Users/yannkerzreho/Documents/SMS/DataViz/out"
analy <- SMSanalyse(outputs_path, compression = 20, graph = T)
