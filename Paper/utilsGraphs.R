###############################################################################
# SOCIAL NETWORKS AS A UTILITY MAXIMIZATION DECISION
# Stat. processing by Aldric Labarthe @ ENS Paris-Saclay
# Licensed under GNU GPL v3 - 2O24
###############################################################################

compute_global_clustering_coefficient <- function(adj_matrix) {
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  return(transitivity(g, type = "global"))
}

importSimulatedGraphs <- function(true_adj_matrix, true_comp_matrix, localpath, clearing=T, mode="optimal", minPrecisionTh=0.05, digitRound=4){
  library(igraph)
  library(ineq)
  library(e1071)
  library(rlist)
  
  results <- list()
  
  if (is.matrix(true_adj_matrix)){
    g <- graph_from_adjacency_matrix(true_adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    degrees <- degree(g)
    if (is.matrix(true_comp_matrix))
      results <- list.append(results, list(id=0, g=g, compmatrix=true_comp_matrix, adjmatrix=true_adj_matrix, degrees=degrees, utilities=NULL))
    else 
      results <- list.append(results, list(id=0, g=g, compmatrix=NULL, adjmatrix=true_adj_matrix, degrees=degrees, utilities=NULL))
    minPrecision <- min(c(true_adj_matrix)[c(true_adj_matrix)!=0])*minPrecisionTh
  }else {
    minPrecision <- 0
  }
  
  dirs <- list.files(path = localpath, pattern = "^sim_", full.names = TRUE)
  
  print(minPrecision)
  for (dir in dirs) {
    id <- as.integer(gsub("^.*sim_([0-9]+)$", "\\1", dir))
    
    pattern <- ""
    pattern_utility <- ""
    if (mode == "optimal"){
      pattern <- "OptimalGraph"
      pattern_utility <- "OptimalUtility"
    }else {
      pattern <- "AdjacencyMatrix"
      pattern_utility <- "Utility"
    }
    file_path <- file.path(dir, sprintf(paste0("SMS-Save-", pattern, "-%d.csv"), id))
    file_path_utility <- file.path(dir, sprintf(paste0("SMS-Save-", pattern_utility, "-%d.csv"), id))
    
    if (!file.exists(file_path)){
      files <- list.files(dir, pattern = paste0("^SMS-Save-", pattern, ".*\\.csv$"), full.names = TRUE)
      
      if (length(files) == 0) {
        stop(paste0("No files found with the pattern '", pattern, "' in 'sim_", id ,"'"))
      }
      
      file_path <- files[1]
    }
    
    if (!file.exists(file_path_utility)){
      files <- list.files(dir, pattern = paste0("^SMS-Save-", pattern_utility, ".*\\.csv$"), full.names = TRUE)
      
      if (length(files) == 0) {
        stop(paste0("No files found with the pattern '", pattern_utility, "' in 'sim_", id ,"'"))
      }
      
      file_path_utility <- files[1]
    }
    
    
    adj_matrix <- fread(file_path, header = FALSE)
    if (pattern == "OptimalGraph")
      adj_matrix <- as.matrix(adj_matrix)[c(1:ncol(adj_matrix)),]
    else 
      adj_matrix <- as.matrix(adj_matrix)[c((nrow(adj_matrix)-ncol(adj_matrix)+1):nrow(adj_matrix)),]
    
    if (pattern == "OptimalGraph"){
      utilityVector <- fread(file_path_utility, header = FALSE)
      utilityVector <- as.numeric(utilityVector$V1)
      comp_matrix <- fread(file_path, header = FALSE)
      comp_matrix <- as.matrix(comp_matrix)[c((ncol(adj_matrix)+1):(2*ncol(adj_matrix))),]
    }else{ 
      utilityVector <- fread(file_path_utility, header = TRUE)
      utilityVector <- as.numeric(utilityVector[utilityVector$round==max(utilityVector$round), ]$utility)
      comp_matrix <- NULL
    }
    
    if (clearing){
      adj_matrix <- ifelse(adj_matrix < minPrecision, 0, adj_matrix)
      adj_matrix <- round(adj_matrix, digits=digitRound)
    }
    
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
    degrees <- degree(g)
    
    results <- list.append(results, list(id=id, g=g, compmatrix=comp_matrix, adjmatrix=adj_matrix, degrees=degrees, utilities=utilityVector))
  }
  
  return(results)
}

importAndCompareSimulatedGraphs <- function(true_adj_matrix, true_comp_matrix, localpath, clearing=T, mode="optimal", minPrecisionTh=0.05, digitRound=4){
  library(igraph)
  library(ineq)
  library(e1071)
  library(data.table)
  
  compute_mse <- function(original_matrix, approximated_matrix) {
    return(mean((original_matrix - approximated_matrix)^2))
  }
  
  graphs <- importSimulatedGraphs(true_adj_matrix, true_comp_matrix, localpath, clearing, mode, minPrecisionTh, digitRound)
  return(compareSimulatedGraphs(graphs))
}

compareSimulatedGraphs <- function(importedGraphs){
  library(igraph)
  library(ineq)
  library(e1071)
  library(data.table)
  
  compute_mse <- function(original_matrix, approximated_matrix) {
    return(mean((original_matrix - approximated_matrix)^2))
  }
  
  graphs <- importedGraphs
  
  results <- data.table(id = integer(), mse = numeric(),
                        rels = numeric(),
                        countDelRel=numeric(), countNewRel=numeric(),countWeightChanges=numeric(),
                        global_clustering_coefficient = numeric(), 
                        global_gini_coefficient=numeric(), global_density=numeric(),
                        mean_distance=numeric(), assortativity=numeric(),
                        deg_mean=numeric(), deg_sec_cmom=numeric(), deg_thd_cmom=numeric(),
                        deg_kurtosis=numeric(), deg_entropy=numeric(),
                        meanComp=numeric(), sdComp=numeric(), meanCorComp=numeric(), sdCorComp=numeric(),
                        utilityMean=numeric(), utilitySd=numeric(), utilityMin=numeric(), utilityMax=numeric(), utilitySum=numeric())
  
  if (graphs[[1]]$id == 0){
    g <- graphs[[1]]$g
    true_comp_matrix <- graphs[[1]]$compmatrix
    true_adj_matrix <- graphs[[1]]$adjmatrix
    
    global_clustering_coefficient <- transitivity(g, type = "global")
    global_density <- edge_density(g)
    mean_distance <- mean_distance(g)
    assortativity <- assortativity_degree(g)
    global_gini_coefficient <- ineq::Gini(E(g)$weight)
    degrees <- degree(g)
    deg_mean <- mean(degrees)
    deg_sec_cmom <- mean((degrees - deg_mean)^2)
    deg_thd_cmom <- mean((degrees - deg_mean)^3)
    deg_kurtosis <- kurtosis(degrees)
    degree_prob <- table(degrees) / length(degrees)
    deg_entropy <- -sum(degree_prob * log(degree_prob))
    
    if (is.matrix(true_comp_matrix)){
      meanComp <- mean(true_comp_matrix)
      sdComp <- sd(true_comp_matrix)
      corComp <- cor(true_comp_matrix)
      meanCorComp <- mean(corComp)
      sdCorComp <- sd(corComp)
    }else {
      meanComp <- NA
      sdComp <- NA
      corComp <- NA
      meanCorComp <- NA
      sdCorComp <- NA
    }
    
    
    results <- rbind(results, data.table(id = 0, mse = NA, rels = sum(true_adj_matrix >0),
                                         countDelRel = NA, countNewRel=NA,countWeightChanges = NA,
                                         global_clustering_coefficient = global_clustering_coefficient, 
                                         global_gini_coefficient=global_gini_coefficient,
                                         assortativity=assortativity,
                                         mean_distance=mean_distance,
                                         global_density=global_density,
                                         deg_mean=deg_mean,deg_sec_cmom=deg_sec_cmom,deg_thd_cmom=deg_thd_cmom,
                                         deg_kurtosis=deg_kurtosis, deg_entropy=deg_entropy,
                                         meanComp=meanComp, sdComp=sdComp, meanCorComp=meanCorComp, sdCorComp=sdCorComp,
                                         utilityMean=NA, utilitySd=NA, utilityMin=NA, utilityMax=NA, utilitySum=NA))
    graphs <- graphs[-1]
  }else {
    true_comp_matrix <- NULL
    true_adj_matrix <- NULL
  }
  
  
  for (graph in graphs) {
    id <- graph$id
    adj_matrix <- graph$adjmatrix
    
    mse <- if (is.matrix(true_adj_matrix)) compute_mse(true_adj_matrix, adj_matrix) else NA
    countDelRel <- if (is.matrix(true_adj_matrix)) sum((true_adj_matrix > 0 & adj_matrix == 0)) else NA
    countNewRel <- if (is.matrix(true_adj_matrix)) sum((true_adj_matrix == 0 & adj_matrix > 0) ) else NA
    countWeightChanges <- if (is.matrix(true_adj_matrix))sum((round(true_adj_matrix,digits = 2) != round(adj_matrix,digits = 2)) & true_adj_matrix > 0 & adj_matrix > 0) else NA
    
    g <- graph$g
    global_clustering_coefficient <- transitivity(g, type = "global")
    global_density <- edge_density(g)
    mean_distance <- mean_distance(g)
    assortativity <- assortativity_degree(g)
    global_gini_coefficient <- ineq::Gini(E(g)$weight)
    degrees <- graph$degrees
    deg_mean <- mean(degrees)
    deg_sec_cmom <- mean((degrees - deg_mean)^2)
    deg_thd_cmom <- mean((degrees - deg_mean)^3)
    deg_kurtosis <- kurtosis(degrees)
    degree_prob <- table(degrees) / length(degrees)
    deg_entropy <- -sum(degree_prob * log(degree_prob))
    
    if (is.matrix(graph$compmatrix)){
      comp_matrix <- graph$compmatrix
      meanComp <- mean(comp_matrix)
      sdComp <- sd(comp_matrix)
      corComp <- cor(comp_matrix)
      meanCorComp <- mean(corComp)
      sdCorComp <- sd(corComp)
    }else {
      meanComp <- NA
      sdComp <- NA
      corComp <- NA
      meanCorComp <- NA
      sdCorComp <- NA
    }
    
    utilityVector <- graph$utilities
    results <- rbind(results, data.table(id = id, mse = mse, rels = sum(adj_matrix >0),
                                         countDelRel = countDelRel, countNewRel=countNewRel,
                                         countWeightChanges = countWeightChanges, 
                                         global_clustering_coefficient = global_clustering_coefficient, 
                                         global_gini_coefficient=global_gini_coefficient,
                                         assortativity=assortativity,
                                         mean_distance=mean_distance,
                                         global_density=global_density,
                                         deg_mean=deg_mean,deg_sec_cmom=deg_sec_cmom,deg_thd_cmom=deg_thd_cmom,
                                         deg_kurtosis=deg_kurtosis, deg_entropy=deg_entropy,
                                         meanComp=meanComp, sdComp=sdComp, meanCorComp=meanCorComp, sdCorComp=sdCorComp,
                                         utilityMean=mean(utilityVector), utilitySd=sd(utilityVector), 
                                         utilityMin=min(utilityVector), utilityMax=max(utilityVector),
                                         utilitySum=sum(utilityVector)))
  }
  
  return(results)
}

buildGraphFromStrangeCSV <- function(filename){
  library(data.table)
  library(igraph)
  
  edges <- fread(filename, header = FALSE)
  
  # Check if there are at least three columns for weights
  if (ncol(edges) >= 3) {
    edges <- edges[, .(V1, V2, V3)]
  } else {
    edges <- edges[, V3 := 1]  # Assign a default weight of 1 if not present
    edges <- edges[, .(V1, V2, V3)]
  }
  
  edges[, 3] <- edges[, 3] / (max(edges[, 3])+0.1*max(edges[, 3]))
  
  g <- graph_from_data_frame(edges, directed = FALSE)
  
  return(list(igraph=g, adj_mat=as_adjacency_matrix(g, attr = "V3", sparse = FALSE)))
}

visGraph <- function(igraphobj){
  library(visNetwork)
  vis_g <- toVisNetworkData(igraphobj)
  visNetwork(nodes = vis_g$nodes, edges = vis_g$edges) %>%
    visEdges(arrows = "to") %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)
}

checksForSolidity <- function(adjmat, N, verticesToKill=1){
  library(igraph)
  library(ineq)
  library(data.table)
  
  I <- nrow(adjmat)
  results <- data.table(id = integer(), global_clustering_coefficient = numeric(), 
                        global_gini_coefficient=numeric(), global_density=numeric(),
                        mean_distance=numeric(), assortativity=numeric())
  
  g <- graph_from_adjacency_matrix(adjmat, mode = "undirected", weighted = TRUE, diag = FALSE)
  global_clustering_coefficient <- transitivity(g, type = "global")
  global_density <- edge_density(g)
  mean_distance <- mean_distance(g)
  assortativity <- assortativity_degree(g)
  global_gini_coefficient <- ineq::Gini(E(g)$weight)
  
  results <- rbind(results, data.table(id = 0,
                                       global_clustering_coefficient = global_clustering_coefficient, 
                                       global_gini_coefficient=global_gini_coefficient,
                                       assortativity=assortativity,
                                       mean_distance=mean_distance,
                                       global_density=global_density))
  
  for (n in 1:N){
    localAdjMat <- adjmat
    for (i in 1:verticesToKill){
      row <- 0
      col <- 0
      while (row == col || as.numeric(adjmat[row,col] == 0)){
        row <- round(runif(1, min = 1, max = I)[1], digits=0)
        col <- round(runif(1, min = 1, max = I)[1], digits=0)
      }
      localAdjMat[row,col] <- 0
      localAdjMat[col,row] <- 0
    }
    
    g <- graph_from_adjacency_matrix(localAdjMat, mode = "undirected", weighted = TRUE, diag = FALSE)
    global_clustering_coefficient <- transitivity(g, type = "global")
    global_density <- edge_density(g)
    mean_distance <- mean_distance(g)
    assortativity <- assortativity_degree(g)
    global_gini_coefficient <- ineq::Gini(E(g)$weight)
    
    results <- rbind(results, data.table(id = n,
                                         global_clustering_coefficient = global_clustering_coefficient, 
                                         global_gini_coefficient=global_gini_coefficient,
                                         assortativity=assortativity,
                                         mean_distance=mean_distance,
                                         global_density=global_density))
  }
  return(results)
}

checksForSolidityWRIndividuals <- function(adjmat, N, individualsToKill=1){
  library(igraph)
  library(ineq)
  library(data.table)
  
  I <- nrow(adjmat)
  results <- data.table(id = integer(), global_clustering_coefficient = numeric(), 
                        global_gini_coefficient=numeric(), global_density=numeric(),
                        mean_distance=numeric(), assortativity=numeric())
  
  g <- graph_from_adjacency_matrix(adjmat, mode = "undirected", weighted = TRUE, diag = FALSE)
  global_clustering_coefficient <- transitivity(g, type = "global")
  global_density <- edge_density(g)
  mean_distance <- mean_distance(g)
  assortativity <- assortativity_degree(g)
  global_gini_coefficient <- ineq::Gini(E(g)$weight)
  
  results <- rbind(results, data.table(id = 0,
                                       global_clustering_coefficient = global_clustering_coefficient, 
                                       global_gini_coefficient=global_gini_coefficient,
                                       assortativity=assortativity,
                                       mean_distance=mean_distance,
                                       global_density=global_density))
  
  for (n in 1:N){
    localAdjMat <- adjmat
    to_delete <- round(runif(individualsToKill, min = 1, max = I), digits=0)
    range_vector <- 1:I 
    result_vector <- setdiff(range_vector, to_delete)
    localAdjMat <- adjmat[result_vector, result_vector]
    
    g <- graph_from_adjacency_matrix(localAdjMat, mode = "undirected", weighted = TRUE, diag = FALSE)
    global_clustering_coefficient <- transitivity(g, type = "global")
    global_density <- edge_density(g)
    mean_distance <- mean_distance(g)
    assortativity <- assortativity_degree(g)
    global_gini_coefficient <- ineq::Gini(E(g)$weight)
    
    results <- rbind(results, data.table(id = n,
                                         global_clustering_coefficient = global_clustering_coefficient, 
                                         global_gini_coefficient=global_gini_coefficient,
                                         assortativity=assortativity,
                                         mean_distance=mean_distance,
                                         global_density=global_density))
  }
  return(results)
}

combine_means <- function(...) {
  data_frames <- list(...)
  
  # Check that all data frames have the same columns
  colnames_list <- lapply(data_frames, colnames)
  if (!all(sapply(colnames_list, identical, colnames_list[[1]]))) {
    stop("All data frames must have the same columns")
  }
  
  mean_list <- list()
  
  # Calculate the mean of each column for each data frame
  for (df in data_frames) {
    mean_vector <- colMeans(df, na.rm = TRUE)
    mean_list <- c(mean_list, list(mean_vector))
  }
  
  # Convert the list of means into a data frame
  mean_df <- as.data.frame(mean_list)
  
  # Set the column names of the new data frame to be the names of the original data frames
  mean_df_colnames <- sapply(substitute(list(...))[-1], deparse)
  colnames(mean_df) <- unlist(lapply(strsplit(mean_df_colnames, split = "[[]"), function(x) x[1]))
  
  return(mean_df)
}

combine_sd <- function(...) {
  data_frames <- list(...)
  
  # Check that all data frames have the same columns
  colnames_list <- lapply(data_frames, colnames)
  if (!all(sapply(colnames_list, identical, colnames_list[[1]]))) {
    stop("All data frames must have the same columns")
  }
  
  # Initialize an empty list to store standard deviation vectors
  sd_list <- list()
  
  # Calculate the standard deviation of each column for each data frame
  for (df in data_frames) {
    sd_vector <- sapply(df, sd, na.rm = TRUE)
    sd_list <- c(sd_list, list(sd_vector))
  }
  
  # Convert the list of standard deviations into a data frame
  sd_df <- as.data.frame(sd_list)
  
  # Set the column names of the new data frame to be the names of the original data frames
  sd_df_colnames <- sapply(substitute(list(...))[-1], deparse)
  colnames(sd_df) <- unlist(lapply(strsplit(sd_df_colnames, split = "[[]"), function(x) x[1]))
  
  return(sd_df)
}

generate_latex_table <- function(means, sds, ignoreCols=c()) {
  means <- round(means, 2)
  sds <- round(sds, 3)
  
  row_names <- rownames(means)
  col_names <- colnames(means)
  
  latex_table <- "\\begin{tabular}{l" 
  latex_table <- paste0(latex_table, paste(rep("c", ncol(means)), collapse = ""), "}\n")
  
  latex_table <- paste0(latex_table, "\\midrule\\midrule\n")
  header <- paste0(" & ", paste(col_names, collapse = " & "), " \\\\")
  latex_table <- paste0(latex_table, header, "\n\\midrule\n")
  
  for (i in 1:nrow(means)) {
    row <- row_names[i]
    if (!(row %in% ignoreCols)){
      row_data <- paste(means[i,], sep = "", collapse = " & ")
      row_line <- paste(row, "&", row_data, "\\\\", sep = " ")
      latex_table <- paste0(latex_table, row_line, "\n")
      row_data <- paste("{\\scriptsize(", sds[i,], ")}", sep = "", collapse = " & ")
      row_line <- paste("", "&", row_data, "\\\\", sep = " ")
      latex_table <- paste0(latex_table, row_line, "\n")
    }
  }
  
  latex_table <- paste0(latex_table, "\\midrule\\midrule\\end{tabular}")
  return(latex_table)
}

