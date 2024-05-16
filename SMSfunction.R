# Titre du fichier : SMSfunction.R

# Auteur : Yann Kerzreho
# Date : 16 mai 2024

# Description : Function to analyze data from SMS simulations

# Exemple:
data_SMSTable <- SMStaTable('~/Outputs')
comp <- compare_runs(c(13, 10, 11, 12, 14), "prop_isolated", data_SMSTable)
comp$plt
comp$ttest

outputs_path = "~/Outputs/2 200i-1000r nogreedy /sim_3971"
dt <- SMS_fulldata(outputs_path, compression = 10, graph = T)


# Core Functions ###############################################################

# dtAgent agregate all individual information from simulation in a folder which
# have folder_path as path. Individual information are compute with AgentData.
dtAgent <- function(folder_path){
  files <- list.files(folder_path, full.names = TRUE)
  dtAgent <- data.frame()
  for(path in files){
    print(path)
    dt <- AgentData(path)
    dtAgent <- rbind(dtAgent, dt)
  }
  dtAgent
}

# SMStaTable computes summary statistics for a different batch of simulations
# folder_path is a folder containing one folder per batch of simulation. A 
# batch can store as many simulation as you want. It return table, the
# SMStaTable_static for each batch and a summary using summarize. 
SMStaTable = function(folder_path){
  #Return all runs
  runs_folder <- list.files(folder_path, pattern = "i-", full.names = TRUE)
  list_table <- list()
  summary <- data.frame()
  #For each folder creat summary stat
  for(i in seq_along(runs_folder)){
    list_table[[i]] = SMStaTable_static(runs_folder[i])
    sum <- summarize(list_table[[i]], index = i)
    #Add some param info
    summary <- rbind(summary, sum)
  }
  list(summary = summary, table = list_table)
}

# compare_runs use a data_SMSTable to compare the variable col_name across
# batch selected by index_run. It perform a series of t.test with perform_t_test
# and return a plot.
compare_runs <- function(index_run, col_name, data_SMSTable){
  n <-  length(index_run)
  data_values <- data.frame()
  list_values <- list()
  for(i in 1:n){
    list_values[[i]] <- data_SMSTable$table[[index_run[i]]][col_name]
    data_values <- rbind(data_values, 
                         data.frame('run' = index_run[i], 
                                    'value' = list_values[[i]]))
  }
  names(data_values) = c('run', 'value')
  data_values$run <- as.factor(data_values$run)
  plt <- ggplot(data = data_values) +
    geom_jitter(aes(x = run, y = value), height = 0, width = 0.2) +
    stat_summary(aes(x = run, y = value), fun = mean, geom = "point", 
                 color = "red", shape = 16, size = 4) +
    stat_summary(aes(x = run, y = value), fun = quantile, 
                 fun.args = list(probs = 0.1), geom = "point", 
                 color = "blue", shape = 16, size = 3) +
    stat_summary(aes(x = run, y = value), fun = quantile, 
                 fun.args = list(probs = 0.9), geom = "point", 
                 color = "blue", shape = 16, size = 3)
  
  ttest <- perform_t_test(list_values, index_run)
  list(plt=plt, ttest=ttest, dt=data_values)
}

# SMS_fulldata takes in input a outputs_path directing to a simulation folder
# and return a data.frame for the edge, vertices, cluster, utiliy, a list
# of personality trait and the networkDynamic (net) as well as the data for
# computing it (data). compression indicates the number of rounds aggregated 
# together and graph whether a render.d3movie should be saved as graph.html in
# your directory.
SMS_fulldata <- function(outputs_path, compression = 0, graph = FALSE){
  library(ndtv)
  library(network)
  library(networkDynamic)
  library(animation)
  ## Find the files, we need a folder as outputs_path
  files_utility <- list.files(outputs_path, pattern = "SMS-Save-Utility",
                              full.names = TRUE)
  if(length(files_utility) != 1) {
    stop("You should have a SMS-Save-Utility file")
  }
  files_edges <- list.files(outputs_path, pattern = "SMS-Save-Edges",
                            full.names = TRUE)
  if(length(files_edges) != 1) {
    stop("You should have a SMS-Save-Edges file")
  }
  files_cluster <- list.files(outputs_path, pattern = "SMS-Save-Clustering",
                              full.names = TRUE)
  if(length(files_edges) != 1) {
    stop("You should have a SMS-Save-Clustering file")
  }
  files_vertices <- list.files(outputs_path, pattern = "SMS-Save-Vertices",
                               full.names = TRUE)
  if(length(files_edges) != 1) {
    stop("You should have a SMS-Save-Vertices file")
  }
  
  ## Preprocess the data
  print("Preprocess the data")
  dfE <- read.csv(files_edges)
  data <- prepareData(dfE) # Transform files_edges to a format networkDynamic 
  # compatible
  dfU <- read.csv(files_utility)
  dfU$agentid <- dfU$agentid + 1 # Normalise ID from 1, not 0
  dfC <- read.csv(files_cluster)
  dfV <- read.csv(files_vertices, header=FALSE) # P is a sequence of 1 and 0,
  # header=FALSE allow to import as string and keep all the P even if one starts
  # with 0s
  dfV <- dfV[-c(1),]
  names(dfV) <- c("round","agentid","gamma","isgreedy", "meandist", "vardist", 
                  "maxdist", "P")
  dfV$agentid <- as.numeric(dfV$agentid) + 1
  dfV <- dfV[dfV$round == 0,]
  list_P <- list() # Creat a vector of personality for each individuals
  for(i in 1:nrow(dfV)){
    list_P[[i]] <- na.omit(as.numeric(strsplit(dfV$P[i], "")[[1]]))
  }
  dfV$mp <- NA
  for(i in 1:nrow(dfV)){
    dfV$mp[i] <- mean(list_P[[i]])
  }
  
  if(compression > 0){ # If ask compress the data, merge some round
    data_comp <- compress_netdata(data, compression)
    dfU$round <- (dfU$round+compression-1)/compression
  }
  row.names(data) <- NULL
  
  ## Creating the net and analyse it
  print("Creating the Net")
  net <- networkDynamic::networkDynamic(edge.spells=data_comp, 
                                        create.TEAs = TRUE,
                                        edge.TEA.names = c('w'))
  
  dfD <- data.frame()
  num_id <- dfV$agentid
  for(i in 1:max(data$onset)){
    d <- data.frame(agentid = num_id, 
                    degree = degree(network.extract(net, at=i)), 
                    round = as.numeric(i))
    dfD <- rbind(dfD,d) # Attention si igraph est installé cela peut poser prblm
  }
  
  dfU <- merge(dfD, dfU, by = c("agentid", "round"))
  dfU <- merge(dfU, dfV[,c('agentid', 'gamma', 'mp')])
  
  if(graph == TRUE){
    library(viridis)
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
      color_index <- cut((slice %e% "w")^(1/2), breaks = length(color_palette), 
                         include.lowest = TRUE)
      color_palette[color_index]
    }
    
    ## Check affinity
    indiv = sample(dfU[(dfU$round == max(dfU$round)) & 
                         (dfU$degree > 5),]$agentid, 1)
    affinity <- as.numeric(lapply(list_P, 
                                  function(x){sum(x*list_P[[indiv]])/100}))
    affinity[indiv] <- 0.2
    af <- (affinity - min(affinity)) / (max(affinity) - min(affinity)) + 0.001
    color_palette_aff <- viridis(1000, option = 'B')
    color <- color_palette_aff[floor(af*1000)]
    color[indiv] <- 0
    
    render.d3movie(net, usearrows = F, displaylabels = F, 
                   render.par=list(tween.frames=2),
                   edge.col =  edge_color_function,
                   edge.lwd = function(slice)(slice %e% "w" + 0.3),
                   vertex.cex = vertex_size_function,
                   vertex.col = color,
                   output.mode=c('HTML'),script.type='embedded',
                   filename='graph.html')
  }
  detach("ndtv")
  detach("network")
  detach("networkDynamic")  
  detach("animation")
  list(dfE=dfE, dfV=dfV, dfC=dfC, dfU=dfU, list_P=list_P, net=net, data=data)
}

# Functional Functions #########################################################

# SMStaTable_static compute statistics with compute_stat for each simulation
# in the folder directed by folder_path and return a data.frame. Only terminus
# network is used: the log named SMS-Save-AdjacencyMatrix
SMStaTable_static = function(folder_path){
  files <- list.files(folder_path, full.names = TRUE)
  stat <- data.frame()
  net <- list()
  for(i in seq_along(files)){
    files_SM <- list.files(files[i], pattern = "SMS-Save-AdjacencyMatrix",
                           full.names = TRUE)
    SM <- read.csv(files_SM, header = F)
    network <- graph_from_adjacency_matrix(as.matrix(SM), weighted = TRUE,
                                           mode = "undirected")
    res <- compute_stat(network)
    net[[i]] <- network
    stat <- rbind(stat, res)
  }
  data.frame(stat)
}

# perform_t_test takes a data frame df, each row will be compared with a t.test,
# index_run is a vector containing the name of the variables in df.
# It return the list of compared variables, the statistic and the p-value
perform_t_test <- function(df, index_run) {
  combinations <- combn(length(df), 2) # Toutes les combinaisons de colonnes
  result <- data.frame(matrix(NA, ncol = 4, nrow = ncol(combinations)))
  colnames(result) <- c("Variable_1", "Variable_2", "t_value", "p_value")
  
  for (i in 1:ncol(combinations)) {
    var1 <- df[[combinations[1, i]]]
    var2 <- df[[combinations[2, i]]]
    test <- t.test(var1, var2)
    result[i, 1] <- index_run[combinations[1, i]]
    result[i, 2] <- index_run[combinations[2, i]]
    result[i, 3] <- test$statistic
    result[i, 4] <- test$p.value
  }
  
  return(result)
}

# summarize creats summary statistics from a SMStaTable_static
summarize <- function(df, index = 0) {
  col_names <- names(df)
  res <- data.frame(index = index)
  val_names <- c("index")
  for(i in seq_along(df)){
    val_names <- c(val_names, paste0(col_names[i],"_1"),
                   paste0(col_names[i],"_m"),
                   paste0(col_names[i],"_9"))
    qt <- quantile(t(df[i]), probs = c(0.1, 0.9))
    res <- data.frame(res, qt[1], mean(t(df[i])), qt[2])
  }
  names(res) <- val_names
  res
}

# compute_stat compute statistics from a static network, return a list with 
# respectively: n - the number of vertice; m - the number of edge; 
# dmean - the mean of degree without isolated; dmin, dmax - the min and max
# degree; meanDist - the mean distance of the graph; gini - the gini coefficient
# compute over degree distribution; gamma - the parameter of the power low fit;
# pval - the p-value (KS) of the godness of fit; rho - assortativity coefficient
# compute with degree; cmean - the mean cluster/transitivity coefficient;
# prop_isolated - the proportion of isolated vertice
compute_stat = function(network){
  degree <- igraph::degree(network)
  n <- length(degree)
  m <- length(igraph::E(network))
  dmin = min(degree[degree > 0])
  dmax = max(degree)
  dmean = mean(degree[degree > 0])
  
  cumulative_degree <- cumsum(sort(degree, decreasing = F)) / sum(degree)
  lorenz_curve <- c(0, cumulative_degree)
  gini <- DescTools::Gini(degree)
  
  # Analyze power law
  pl <- igraph::fit_power_law(degree)
  # Assortativity
  rho = igraph::assortativity_degree(network)
  # Diameter and radius 
  ecc = igraph::eccentricity(network)
  diam = max(ecc)
  rad = min(ecc)
  # Density
  edge_density = igraph::edge_density(network)
  # Clustering
  c = igraph::transitivity(network, type = "localundirected", isolates = 'NaN')
  cmean = mean(c, na.rm = T)
  # similarity?
  # Distance
  meanDist = igraph::mean_distance(network)
  
  prop_isolated = sum(degree == 0)/n
  
  list(n=n, m=m, dmin=dmin, dmean=dmean, dmax=dmax, meanDist=meanDist, 
       gini=gini, gamma = pl$alpha, pval = pl$KS.p, rho=rho, cmean = cmean, 
       prop_isolated = prop_isolated)
}

# prepareData transform a log SMS-Save-Edges to a networkDynamic format
prepareData <- function(rData){
  from_to <- data.frame(tail = apply(rData[,c(2,3)], 1, min),
                        head = apply(rData[,c(2,3)], 1, max))
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
    if(nd == 1 & d$new_weight[1] > 0) res <- rbind(res, 
                                                   c(d$round[1], Inf, tail, 
                                                     head, d$new_weight[1]))
    if(nrow(res)>0) names(res) <- c("onset", "terminus", "tail", "head", "w")
    Data <- rbind(Data, res)
  }
  Data
}

# compress_netdata compress a prepareData by aggregating rounds
compress_netdata <- function(data, by = 5){
  tmin <- min(data$onset) + 1
  tmax <- max(data$onset) + 1
  tseq <- seq(tmin, tmax, by)
  res <- data.frame()
  for(i in 1:(length(tseq)-1)){
    d <- data[data$onset %in% seq(tseq[i], tseq[i+1]-1), ]
    d <- d[d$terminus > tseq[i+1], ]
    d$onset <- i
    d$terminus <- floor(d$terminus / by)
    res <- rbind(res, d)
  }
  res
}

# AgentData takes a path of a simulation folder to creat a data.frame that
# summarize individuals at the terminus. It return the agent id, the utility
# the gamma of the individual, mp the proportion of 1 in its personality vector, 
# and the eccentricity, cluster coefficient and the degree.
AgentData <- function(simulation_path){
  files_SM <- list.files(simulation_path, pattern = "SMS-Save-AdjacencyMatrix",
                         full.names = TRUE)
  files_Ut <- list.files(simulation_path, pattern = "SMS-Save-Utility",
                         full.names = TRUE)
  files_Vr <- list.files(simulation_path, pattern = "SMS-Save-Vertices", 
                         full.names = TRUE)
  
  dfU <- read.csv(files_Ut)
  dfU <- dfU[dfU$round == max(dfU$round),]
  dfU$agentid <- dfU$agentid + 1 # Normalise ID from 1, not 0
  
  dfV <- read.csv(files_Vr, header=FALSE) # P is a sequence of 1 and 0, 
  # header=FALSE allow to import as string and keep all the P even if one starts
  # with 0s
  dfV <- dfV[-c(1),]
  names(dfV) <- c("round","agentid","gamma","isgreedy", "meandist", "vardist",
                  "maxdist", "P")
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
  
  SM <- read.csv(files_SM, header = F)
  network <- igraph::graph_from_adjacency_matrix(as.matrix(SM), weighted = TRUE,
                                         mode = "undirected")
  list_degree <- igraph::degree(network)
  list_cluster <- igraph::transitivity(network, type = "localundirected",
                               isolates = 'NaN')
  list_eccentricity <- igraph::eccentricity(network)
  
  data.frame(agentid = dfU$agentid, utility = dfU$utility, 
             gamma = as.numeric(dfV$gamma), mp = dfV$mp, 
             eccentricity = list_eccentricity, 
             coefCl = list_cluster, degree = list_degree)
}

