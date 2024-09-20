###############################################################################
# SOCIAL NETWORKS AS A UTILITY MAXIMIZATION DECISION
# Stat. processing by Aldric Labarthe @ ENS Paris-Saclay
# Licensed under GNU GPL v3 - 2O24
###############################################################################

rebuildCompMatrix <- function(adjMatrix, kappa, gamma, delta){
  I <- nrow(adjMatrix)
  M <- matrix(0, nrow=I, ncol=I)
  
  sum_rows <- rowSums(adjMatrix)
  sum_cols <- colSums(adjMatrix)
  
  for (i in 1:I) {
    for (j in i:I) {
      if (i != j) { 
        A_ij <- adjMatrix[i, j]
        
        term1 <- (gamma * A_ij^(gamma - 1)) / (1 - A_ij^gamma)^2
        term2 <- (gamma * A_ij^(gamma - 1)) / (1 - A_ij^gamma)^2
        term3 <- delta * sum_rows[i]^(delta - 1)
        term4 <- delta * sum_cols[j]^(delta - 1)
        
        M[i, j] <- (1 / (2 * kappa)) * (term1 + term2 + term3 + term4)
        M[j, i] <- (1 / (2 * kappa)) * (term1 + term2 + term3 + term4)
      }else {
        M[i, j] <-0
      }
    }
  }
  return (M)
}

checkLoveTriangle <- function(compMatrix){
  single_attempts = 0
  single_successes = 0
  graphsize = sqrt(length(compMatrix))
  issues = 0
  for (i in 1:graphsize){
    check = 0
    attempts = 0
    for (m in 1:graphsize){
      for (k in 1:graphsize){
        if(m!=k && k!= i && m != i && compMatrix[m,k] < compMatrix[i,k] + compMatrix[m,i]){
          check = check+1
          single_successes = single_successes+1
        }
        if(m!=k && k!= i && m != i){
          attempts = attempts+1
          single_attempts = single_attempts+1
        }
      }
    }
    if (check <1)
      issues = issues+1
  }
  print(paste("Issues=", issues, "ProbaTotal=", issues/(graphsize), "ProbaIndiv",single_successes/single_attempts))
  
}

insertWNinCompMatrix <- function(compMatrix, sd=0.1){
  I <- nrow(compMatrix)
  for (i in 1:I) {
    for (j in i:I) {
      if (i != j) {
        compMatrix[i, j] <- compMatrix[i, j] + rnorm(1, mean=0, sd=sd)[1]
        compMatrix[j, i] <- compMatrix[i, j]
      }

    }
  }
  return (compMatrix)
}
# DEPRECATED 
# assessDeviationInPopulation <- function(localpath){
#   # Load necessary libraries
#   library(igraph)
#   library(ineq)
#   library(data.table)
#   
#  
#   # Initialize a dataframe to store results
#   results <- data.table(id = integer(), global_clustering_coefficient = numeric(), 
#                         global_gini_coefficient=numeric(), global_density=numeric(),
#                         mean_distance=numeric(), assortativity=numeric())
#   
#   # Get the list of directories
#   dirs <- list.files(path = localpath, pattern = "^sim_", full.names = TRUE)
#   
#   #print(minPrecision)
#   # Loop through each directory
#   for (dir in dirs) {
#     # Get the id from the directory name
#     id <- as.integer(gsub("^.*sim_([0-9]+)$", "\\1", dir))
#     
#     # Construct the file path
#     file_path <- file.path(dir, sprintf("SMS-Save-OptimalGraph-%d.csv", id))
#     
#     # Check if the file exists
#     if (file.exists(file_path)) {
#       # Read the adjacency matrix from the file
#       adj_matrix <- fread(file_path, header = FALSE)
#       adj_matrix <- as.matrix(adj_matrix)[c(1:ncol(adj_matrix)),]
# 
#       # Compute the global clustering coefficient
#       g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
#       global_clustering_coefficient <- transitivity(g, type = "global")
#       global_density <- edge_density(g)
#       mean_distance <- mean_distance(g)
#       assortativity <- assortativity_degree(g)
#       global_gini_coefficient <- ineq::Gini(E(g)$weight)
#       
#       # Add the results to the dataframe
#       results <- rbind(results, data.table(id = id, 
#                                            global_clustering_coefficient = global_clustering_coefficient, 
#                                            global_gini_coefficient=global_gini_coefficient,
#                                            assortativity=assortativity,
#                                            mean_distance=mean_distance,
#                                            global_density=global_density))
#     }
#   }
#   
#   return(results)
# }

buildCompMatrix <- function(pmatrix){
  I <- nrow(pmatrix)
  dimP <- ncol(pmatrix)
  compMat <- matrix(nrow=I, ncol=I)
  for (i in 1:I) {
    for (j in i:I) {
      compMat[i, j] <- 0
      if (i != j){
        for (k in 1:(dimP)) {
          compMat[i, j] <- compMat[i, j] + pmatrix[i, k] * pmatrix[j, k]
        }
        compMat[i, j] <- compMat[i, j] #* (100 / dimP)
        compMat[j, i] <- compMat[i, j]
      }
    }
  }
  return(compMat)
}

fromPbuildCompMatrix <- function(pmatrixfile, sd=0, components=0){
  pmatrix <- if (!is.matrix(pmatrixfile)) as.matrix(read.csv(pmatrixfile, header = FALSE)) else pmatrixfile
  I <- nrow(pmatrix)
  dimP <- ncol(pmatrix)
  
  if (sd > 0){
    for (i in 1:I) {
      if (components == 0){
        wn <- rnorm(1, mean=0, sd=sd)[1]
        for (j in 1:dimP) {
            pmatrix[i, j] <- pmatrix[i, j] + wn
        }
      }
    }
    if (components > 0) {
      new_cols <- matrix(rnorm(I*I, mean = 0, sd = sd), nrow = I, ncol = components)
      pmatrix <- cbind(pmatrix, new_cols)
    }
  }
  compMat <- buildCompMatrix(pmatrix)
  return(compMat)
}

fromPbuildBarycentricCompMatrix <- function(pmatrixfile, N) {
  pmatrix <- as.matrix(read.csv(pmatrixfile, header = FALSE))
  # Number of rows and columns in the matrix
  I <- nrow(pmatrix)
  dimP <- ncol(pmatrix)
  
  # Check if N is greater than the number of rows
  if (N > I) 
    stop("N cannot be greater than the number of rows in pmatrix.")
  
  newpmatrix <- pmatrix
  # Select N random rows to replace
  rows_to_replace <- sample(1:I, N)
  print(paste("We alter", rows_to_replace))
  
  for (i in rows_to_replace) {
    # Generate random weights that sum to 1
    weights <- runif(I, min = 0, max = 1)
    weights <- weights / sum(weights)
    
    # Create a new row as a barycentric combination of all rows
    new_row <- colSums(pmatrix * weights)
    
    # Replace the ith row with the new row
    newpmatrix[i, ] <- new_row
  }
  
  compMat <- buildCompMatrix(newpmatrix)
  return(compMat)
}

fromPbuildNewBarycentricCompMatrix <- function(pmatrixfile, N, symOrder) {
  pmatrix <- as.matrix(read.csv(pmatrixfile, header = FALSE))
  # Number of rows and columns in the matrix
  I <- nrow(pmatrix)
  dimP <- ncol(pmatrix)
  
  # Check if N is greater than the number of rows
  if (N > I) 
    stop("N cannot be greater than the number of rows in pmatrix.")
  
  newpmatrix <- pmatrix
  # Select N random rows to replace
  rows_to_replace <- sample(1:I, N)
  print(paste("We alter", rows_to_replace))
  
  for (i in rows_to_replace) {
    syms <- if (N < 0.75 * I) sample(setdiff(1:I, rows_to_replace), symOrder) else sample(1:I, symOrder)
    # Generate random weights that sum to 1
    #weights <- runif(I, min = 0, max = 5)
    
    weights <- runif(I, min = 0, max = 0)
    #weights[i] <- 6
    t <- 0
    for (s in syms){
      #weights[s] <- runif(1, min = 0, max = symOrder)[1]
      #weights[s] <- (symOrder)-t
      #t <- t+1
      weights[s] <- 1
    }
    weights <- weights / sum(weights)

    # Create a new row as a barycentric combination of all rows
    new_row <- colSums(pmatrix * weights)
    
    # Replace the ith row with the new row
    newpmatrix[i, ] <- new_row
  }
  
  originalCompMat <- buildCompMatrix(pmatrix)
  compMat <- buildCompMatrix(newpmatrix)
  
  print(paste(mean(compMat), mean(originalCompMat), sd(colMeans(compMat)), sd(colMeans(originalCompMat))))
  
  return(compMat)
}

fromPbuildNewCoupleCompMatrix <- function(pmatrixfile, N, adjmat, homophilyScore=5) {
  pmatrix <- as.matrix(read.csv(pmatrixfile, header = FALSE))
  # Number of rows and columns in the matrix
  I <- nrow(pmatrix)
  dimP <- ncol(pmatrix)
  
  # Check if N is greater than the number of rows
  if (N > I) 
    stop("N cannot be greater than the number of rows in pmatrix.")
  
  newpmatrix <- pmatrix
  originalCompMat <- buildCompMatrix(pmatrix)
  # Select N random rows to replace
  rows_to_replace <- sample(1:I, N)
  print(paste("We alter", rows_to_replace))
  
  alreadySelected <- c()
  alreadySelectedTwice <- c()
  previous_matches <- numeric(I)
  
  for (i in rows_to_replace) {
    firstParent <- if (N < 0.75 * I) sample(setdiff(1:I, rows_to_replace), 1) else sample(1:I, 1)

    #row_means <- adjmat[,firstParent]
    #differences <- adjmat[,firstParent]
    
    #differences[firstParent] <- 0
    #for (a in alreadySelected)
    #  differences[a] = -Inf
    
    #prob <- exp(differences * homophilyScore) / sum(exp(differences * homophilyScore))

    row_means <- apply(originalCompMat, 1, mean)
    differences <- abs(row_means - row_means[firstParent])
    differences[firstParent] <- Inf
    for (a in alreadySelectedTwice)
      differences[a] = Inf
    
    prob <- exp(-differences * homophilyScore) / sum(exp(-differences * homophilyScore))

    # Draw the second row based on the calculated probabilities
    secondParent <- firstParent
    while(secondParent == firstParent || secondParent %in% alreadySelectedTwice || previous_matches[secondParent] == firstParent)
      secondParent <- sample(1:I, 1, prob = prob)
    if (!(secondParent %in% alreadySelected))
      alreadySelected <- c(alreadySelected, secondParent)
    else 
      alreadySelectedTwice <- c(alreadySelectedTwice, secondParent)
    previous_matches[firstParent] <- secondParent
    #secondParent <- which.min(differences)
    
    # row_means <- adjmat[,firstParent]
    # secondParent <- which.max(row_means)
    
    print(paste("Parent 1", row_means[firstParent], "Parent 2", row_means[secondParent]))

    weights <- runif(I, min = 0, max = 0)
    weights[firstParent] <- I*0.5#*0.45
    weights[secondParent] <- I*0.5#*0.45
    
    weights <- weights / sum(weights)
    
    # Create a new row as a barycentric combination of all rows
    new_row <- colSums(pmatrix * weights)
    
    # Replace the ith row with the new row
    newpmatrix[i, ] <- new_row
  }
  
  compMat <- buildCompMatrix(newpmatrix)
  
  print(paste(mean(compMat), mean(originalCompMat), sd(colMeans(compMat)), sd(colMeans(originalCompMat))))

  return(compMat)
}

reduceCompatibilitiesAtThreshold <- function(adj_mat, comp_mat){
  if (!is.matrix(comp_mat))
    comp_mat <- as.matrix(read.csv(comp_mat, header = FALSE))
  I <- nrow(comp_mat)
  for (i in 1:I) {
    for (j in 1:I) {
      if (i != j) {
        if (adj_mat[i, j] == 0){
          comp_mat[i, j] <- max(comp_mat[i, j] - (runif(1, min = 0.01, max = 0.9)[1])*comp_mat[i, j], 0)
          comp_mat[j, i] <- comp_mat[i, j]
        }
      }
    }
  }
  return (comp_mat)
}

shuffleColumns <- function(mat) {
  shuffled_mat <- apply(mat, 2, function(col) sample(col))
  return(shuffled_mat)
}

fromPbuildShuffleCompMatrix <- function(pmatrixfile) {
  pmatrix <- as.matrix(read.csv(pmatrixfile, header = FALSE))
  oldcompMat <- buildCompMatrix(pmatrix)
  newpmatrix <- shuffleColumns(pmatrix)
  compMat <- buildCompMatrix(newpmatrix)
  if (min(compMat) < 0)
    compMat <- compMat - min(compMat)
  compMat <- compMat * max(oldcompMat)  / max(compMat)
  return(compMat)
}

buildInvisibleHandMatrix <- function(matrix) {
  library(Matrix)
  dim <- nrow(matrix)
  
  # Sample dim unique values from the matrix to create perm_vector
  unique_values <- unique(as.vector(matrix))
  if (length(unique_values) < dim) {
    stop("Not enough unique values in the matrix to sample from.")
  }
  perm_vector <- sample(unique_values, dim)
  
  # Step 2: Function to find the nearest value in perm_vector
  nearest_value <- function(x, perm_vector) {
    perm_vector[which.min(abs(perm_vector - x))]
  }
  
  # Step 3: Transform the matrix by replacing each element with the nearest value in perm_vector
  transformed_matrix <- apply(matrix, c(1, 2), nearest_value, perm_vector = perm_vector)
  
  for (i in 1:dim)
    transformed_matrix[i, i] <- 0

  return(transformed_matrix)
}

buildInvisibleHandMatrix2 <- function(matrix) {
  library(Matrix)
  dim <- nrow(matrix)
  
  # We select a column of the matrix to be our candidate 
  origin <- sample.int(dim, 1)[1]
  perm_vector <- matrix[, origin]
  
  # We remove the 0 in the column
  perm_vector <- perm_vector[-origin]
  
  
  insert_value <- function(vec, value, pos) {
    # Combine the parts of the vector with the new value inserted
    if (pos > 1 && pos < length(vec)){
      return(c(vec[1:(pos-1)], value, vec[pos:length(vec)]))
    }else if (pos == 1){
      return(c(value, vec[1:length(vec)]))
    }else {
      return(c(vec[1:length(vec)], value))
    }
  }
  
  print(perm_vector)
  final_matrix <- matrix
  print(final_matrix)
  for (i in 1:(dim-1)){
    print(paste0("We work on the ", i, "th column"))
    local_perm_vector <- perm_vector
    if (i > 1){
      to_delete <- final_matrix[1:(i-1),i]
        
      remove_first_occurrence <- function(vec, value) {
          pos <- which(vec == value)[1]  # Find the first occurrence
          if (!is.na(pos)) {
            vec <- vec[-pos]  # Remove the element at that position
          }else{
            stop(paste("Not found!", value, vec))
          }
          return(vec)
      }
      
      for (value in to_delete) {
        local_perm_vector <- remove_first_occurrence(local_perm_vector, value)
      }
    }
    if (length(local_perm_vector) > 1){
      # We generate permutations of this vector 
      permutations <- replicate(min(dim*10000, factorial(length(local_perm_vector)) *10), sample(local_perm_vector))
      
      errors <- 0
      permutations_diff <- abs(permutations - final_matrix[(i+1):dim,i])
      permutations_diff_sums <- colSums(permutations_diff)
      while(dim(permutations)[[2]] > errors){
        final_perm_id <- which.min(permutations_diff_sums)
        error <- FALSE
        if (i > 1){
          for (z in 1:(length(permutations[, final_perm_id]))){
            if (permutations[z, final_perm_id] %in% final_matrix[1:(i-1),i+z-1]){
              error <- TRUE
              break
            }
          }
        }
        if (!error){
          break
        }else {
          errors <- errors+1
          #print("Failure with")
          #print(permutations_local[,final_perm_id] )
          permutations_diff_sums[final_perm_id] <- Inf
        }
      }
      if (error == TRUE)
        stop("Unable to find a candidate matrix, please restart the process")
      
      final_row <- if (i>1) c(final_matrix[1:(i-1),i], 0, permutations[, final_perm_id]) else c(0, permutations[, final_perm_id])
    }else {
      final_row <- c(final_matrix[1:(i-1),i], 0, local_perm_vector)
    }
    final_matrix[, i] <- final_row
    final_matrix[i, ] <- final_row
    #print(final_matrix)
  }
 
  return(final_matrix)
}

buildInvisibleHandMatrix3 <- function(matrix, nb, max_attempts=2000){
  if (dim(matrix)[[2]] != dim(matrix)[[1]])
    stop("We only work with squared matrices")
  
  if (dim(matrix)[[2]] %%2 != 0)
    stop("We can only work with even matrices")
  
  dimMat <- dim(matrix)[[2]]
  
  if (factorial(dimMat) < nb)
    stop("You ask for more matrices that I can build")
  
  nb <- min(factorial(dimMat),nb)
  originalVector <- 2:dimMat
 
  library(permute)
  library(parallel)
  library(rlist)
  compute_mse <- function(original_matrix, approximated_matrix) {
    return(mean((original_matrix - approximated_matrix)^2))
  }
  
  remove_first_occurrence <- function(vec, value) {
    pos <- which(vec == value)[1]  # Find the first occurrence
    if (!is.na(pos)) {
      vec <- vec[-pos]  # Remove the element at that position
    }else{
      stop(paste("Not found!", value, vec))
    }
    return(vec)
  }
  
  generate_unique_permutations <- function(vec, N) {
    num_cores <- detectCores() - 4
    cl <- makeCluster(num_cores)
    unique_permutations <- vector("list", N)
    
    perm_set <- new.env()
    
    generate_and_check <- function(id, vec, perm_set) {
      while (TRUE) {
        perm <- sample(vec)
        perm_key <- paste(perm, collapse = "-")
        if (!exists(perm_key, envir = perm_set)) {
          assign(perm_key, TRUE, envir = perm_set)
          return(perm)
        }
      }
    }
    
    for (i in 1:N) {
      unique_permutations[[i]] <- parLapply(cl, 1, generate_and_check, vec, perm_set)[[1]]
    }
    
    stopCluster(cl)
    do.call(rbind, unique_permutations)
  }
  first_permutations <- generate_unique_permutations(originalVector, nb)

  final_results <- list()
  
  for (perm_id in 1:nb){
    perm <- c(1, first_permutations[perm_id, ])
    local_matrix_mask <- matrix(, nrow=dimMat, ncol=dimMat)
    local_matrix_mask[1,] <- perm
    local_matrix_mask[,1] <- perm
    
    # We build the antidiagonals until reaching the main antidiagonal
    for (i in 2:(dimMat)){
      local_matrix_mask[i,i] <- 1
      if (i <= dimMat/2){
        for (j in (i+1):(dimMat-i+1)){
          local_matrix_mask[i,j] <- local_matrix_mask[i-1,j+1]
          local_matrix_mask[j,i] <- local_matrix_mask[i,j]
        }
      }
    }
    #print(local_matrix_mask)
    # We fill the lower anti-triangle of the matrix
    for (i in 2:(dimMat-1)){
      print(paste("Working on row", i))
      existing_row <- local_matrix_mask[i, 1:max((dimMat-i+1), i)]
      remaining_possibilities <- setdiff(c(1, originalVector), existing_row)
      if (length(remaining_possibilities) == 1){
        local_matrix_mask[i, dimMat] <- remaining_possibilities[1]
        local_matrix_mask[dimMat, i] <- remaining_possibilities[1]
      }else {
        error <- TRUE
        candidate <- sample(remaining_possibilities, length(remaining_possibilities))
        while(error){
          errs <- 0
            for (z in 1:(length(candidate))){
              if (candidate[z] %in% local_matrix_mask[1:(i-1),(max((dimMat-i+1), i)+z)]){
                errs <- errs+1
                break
              }
            }
          
          if (errs == 0)
            error <- FALSE
          else
            candidate <- sample(remaining_possibilities, length(remaining_possibilities))
        }
        local_matrix_mask[i, (max((dimMat-i+1), i)+1):dimMat] <- candidate
        local_matrix_mask[(max((dimMat-i+1), i)+1):dimMat, i] <- candidate
      }
    }
    
    #We select the permutation of row and columns that better fit our mask
    attempt <- 0
    best_matrix <- matrix(, ncol=dimMat, nrow=dimMat)
    best_mse <- Inf
      
    while(attempt < max_attempts){
      cols <- sample(dimMat)
      matrix_candidate <- matrix[cols, cols]
      donnor_vector <- remove_first_occurrence(matrix_candidate[1,], 0)
      donnor_vector <- c(0, donnor_vector)
      local_matrix <- apply(local_matrix_mask, c(1, 2), function(x){
        return(donnor_vector[x])
      })
      mse <- compute_mse(matrix_candidate, local_matrix)
      if (mse < best_mse){
        best_matrix <- local_matrix
        best_mse <- mse
      }
      attempt <- attempt+1
    }
    final_results <- list.append(final_results, list(mask=local_matrix_mask, candidate=matrix_candidate, result=best_matrix, mse=best_mse))
  }
  return(final_results)
}


buildInvisibleHandMatrix4 <- function(matrix, nb, max_attempts=2000){
  if (dim(matrix)[[2]] != dim(matrix)[[1]])
    stop("We only work with squared matrices")
  
  if (dim(matrix)[[2]] %%2 != 0)
    stop("We can only work with even matrices")
  
  dimMat <- dim(matrix)[[2]]
  
  if (factorial(dimMat) < nb)
    stop("You ask for more matrices that I can build")
  
  nb <- min(factorial(dimMat),nb)
  originalVector <- 2:dimMat
  
  library(permute)
  library(parallel)
  library(rlist)
  compute_mse <- function(original_matrix, approximated_matrix) {
    return(mean((original_matrix - approximated_matrix)^2))
  }
  
  remove_first_occurrence <- function(vec, value) {
    pos <- which(vec == value)[1]  # Find the first occurrence
    if (!is.na(pos)) {
      vec <- vec[-pos]  # Remove the element at that position
    }else{
      stop(paste("Not found!", value, vec))
    }
    return(vec)
  }
  
  generate_unique_permutations <- function(vec, N) {
    num_cores <- detectCores() - 4
    cl <- makeCluster(num_cores)
    unique_permutations <- vector("list", N)
    
    perm_set <- new.env()
    
    generate_and_check <- function(id, vec, perm_set) {
      while (TRUE) {
        perm <- sample(vec)
        perm_key <- paste(perm, collapse = "-")
        if (!exists(perm_key, envir = perm_set)) {
          assign(perm_key, TRUE, envir = perm_set)
          return(perm)
        }
      }
    }
    
    for (i in 1:N) {
      unique_permutations[[i]] <- parLapply(cl, 1, generate_and_check, vec, perm_set)[[1]]
    }
    
    stopCluster(cl)
    do.call(rbind, unique_permutations)
  }
  first_permutations <- generate_unique_permutations(originalVector, nb)
  
  final_results <- list()
  
  for (perm_id in 1:nb){
    perm <- c(1, first_permutations[perm_id, ])
    local_matrix_mask <- matrix(, nrow=dimMat, ncol=dimMat)
    local_matrix_mask[1,] <- perm
    local_matrix_mask[,1] <- perm
    
    # We build the antidiagonals until reaching the main antidiagonal
    for (i in 2:(dimMat)){
      local_matrix_mask[i,i] <- 1
      if (i <= dimMat/2){
        for (j in (i+1):(dimMat-i+1)){
          local_matrix_mask[i,j] <- local_matrix_mask[i-1,j+1]
          local_matrix_mask[j,i] <- local_matrix_mask[i,j]
        }
      }
    }
    # We fill the lower anti-triangle of the matrix
    buildTriangle <- function(local_matrix_mask, dimMat){
      for (i in 2:(dimMat-1)){
        #print(local_matrix_mask)
        print(paste("Working on row", i))
        existing_row <- local_matrix_mask[i, 1:max((dimMat-i+1), i)]
        remaining_possibilities <- setdiff(c(1, originalVector), existing_row)
        if (length(remaining_possibilities) != dimMat-max((dimMat-i+1), i))
          return(NULL)
        if (length(remaining_possibilities) == 1){
          local_matrix_mask[i, dimMat] <- remaining_possibilities[1]
          local_matrix_mask[dimMat, i] <- remaining_possibilities[1]
        }else {
          individual_possibilities <- list()
          for (z in 1:(dimMat-max((dimMat-i+1), i))){
            individual_possibilities <- list.append(individual_possibilities, setdiff(remaining_possibilities, local_matrix_mask[1:(i-1),(max((dimMat-i+1), i)+z)]))
          }
          for (a in 1:length(individual_possibilities)){
            for (z in 1:length(individual_possibilities)){
              if (length(individual_possibilities[[z]]) == 1){
                for (zbis in 1:length(individual_possibilities)){
                  if (z != zbis && length(individual_possibilities[[zbis]]) > 1 && individual_possibilities[[z]][1] %in% individual_possibilities[[zbis]]){
                    individual_possibilities[[zbis]] <- setdiff(individual_possibilities[[zbis]], individual_possibilities[[z]][1])
                  }
                  if (z != zbis && length(individual_possibilities[[zbis]]) == 1 && individual_possibilities[[zbis]][1] == individual_possibilities[[z]][1]){
                    return(NULL)
                  }
                  if (length(individual_possibilities[[zbis]]) == 0)
                    return(NULL)
                }
              }
              if (length(individual_possibilities[[z]]) == 0)
                return(NULL)
            }
          }
          #print(individual_possibilities)
          candidate <- c()
          completed <- FALSE
          attempts <- 0
          while(!completed){
            #print("Building candidate")
            for (z in 1:length(individual_possibilities)){
              if (length(individual_possibilities[[z]]) == 1){
                if (individual_possibilities[[z]][1] %in% candidate){
                  break
                }else {
                  candidate <- c(candidate, individual_possibilities[[z]][1])
                }
              }else {
                #print("Remaining possibilities")
                local_candidate <- setdiff(individual_possibilities[[z]], candidate)
                #print(local_candidate)
                if (length(local_candidate) == 0){
                  break
                }else if (length(local_candidate) == 1){
                  if (local_candidate[1] %in% candidate){
                    break
                  }else {
                    candidate <- c(candidate, local_candidate[1])
                  }
                }else {
                  if (z < length(individual_possibilities) && length(intersect(local_candidate, setdiff(individual_possibilities[[z+1]], candidate))) != length(local_candidate))
                    candidate <- c(candidate, sample(setdiff(local_candidate, intersect(local_candidate, individual_possibilities[[z+1]])), 1)[1])
                  else
                    candidate <- c(candidate, sample(local_candidate, 1)[1])
                }
              }
              #print("Candidate at this stage")
              #print(candidate)
            }
            #print(candidate)
            if (is.null(candidate))
              return(NULL)
            if (length(candidate) == length(individual_possibilities))
              completed <- TRUE
            else
              candidate <- c()
            attempts <- attempts+1
            if (attempts > 20000)
              break
          }
          print(local_matrix_mask)
          if (length(candidate) != length(individual_possibilities))
            return(NULL)
          #print("selected candidate")
          #print(candidate)
          local_matrix_mask[i, (max((dimMat-i+1), i)+1):dimMat] <- candidate
          local_matrix_mask[(max((dimMat-i+1), i)+1):dimMat, i] <- candidate
        }
      }
      return(local_matrix_mask)
    }
    
    final_local_matrix_mask <- NULL
    while (!is.matrix(final_local_matrix_mask)){
      print("Starting process")
      final_local_matrix_mask <- buildTriangle(local_matrix_mask, dimMat)
      break
    }
    print(final_local_matrix_mask)
    stop()
    local_matrix_mask <- final_local_matrix_mask
    
    #We select the permutation of row and columns that better fit our mask
    attempt <- 0
    best_matrix <- matrix(, ncol=dimMat, nrow=dimMat)
    best_mse <- Inf
    
    while(attempt < max_attempts){
      cols <- sample(dimMat)
      matrix_candidate <- matrix[cols, cols]
      donnor_vector <- remove_first_occurrence(matrix_candidate[1,], 0)
      donnor_vector <- c(0, donnor_vector)
      local_matrix <- apply(local_matrix_mask, c(1, 2), function(x){
        return(donnor_vector[x])
      })
      mse <- compute_mse(matrix_candidate, local_matrix)
      if (mse < best_mse){
        best_matrix <- local_matrix
        best_mse <- mse
      }
      attempt <- attempt+1
    }
    final_results <- list.append(final_results, list(mask=local_matrix_mask, candidate=matrix_candidate, result=best_matrix, mse=best_mse))
  }
  return(final_results)
}

makeMatrixSymmetric <- function(matrix_input) {
  # Check if the input is a matrix
  if (!is.matrix(matrix_input)) {
    stop("Input must be a matrix.")
  }
  
  # Make the matrix symmetric
  symmetric_matrix <- (matrix_input + t(matrix_input)) / 2
  
  return(symmetric_matrix)
}

buildCompMatrixFromRelationVector <- function(relationsVector, reduction=0.5){
  remove_first_occurrence <- function(vec, value) {
    pos <- which(vec == value)[1]  # Find the first occurrence
    if (!is.na(pos)) {
      vec <- vec[-pos]  # Remove the element at that position
    }else{
      stop(paste("Not found!", value, vec))
    }
    return(vec)
  }
  
  dim <- length(relationsVector)
  adjMat <- matrix(, nrow=dim, ncol=dim)
  for (row in 1:dim){
    for (col in 1:dim){
      adjMat[row, col] <- 0
    }
  }

  relationsVector <- remove_first_occurrence(relationsVector, 0)
  for (row in 1:dim){
    adjMat[row, row] <- 0
    if (row < dim)
      adjMat[row, (row+1):dim] <- relationsVector[1:(dim-row)]
    if (row > 1)
      adjMat[row, 1:(row-1)] <- relationsVector[(dim-row+1):(dim-1)]
  }
  
  adjMat <- makeMatrixSymmetric(adjMat)
  
  compMat <- rebuildCompMatrix(adjMat, 10, 9, 2)
  for (row in 1:dim){
    for (col in 1:dim){
      if (adjMat[row, col] == 0 && row != col)
        compMat[row, col] <- compMat[row, col]*reduction
    }
  }
  return(list(compMat=compMat, adjMat=adjMat))
}

isAValidIHMatrix <- function(mat) {
  # Check if the input is a valid matrix
  if (!is.matrix(mat)){
    print("Input must be a matrix")
    return(FALSE)
  } 
  
  if (dim(mat)[[1]] != dim(mat)[[2]]){
    print("Matrix must be squarred")
    return(FALSE)
  } 

  for (i in 1:(dim(mat)[[1]]-1)){
    for (j in (i+1):dim(mat)[[2]]){
      if (abs(tmp[i,j] - tmp[j, i]) > 1e-6){
        print(paste("Matrix must be symmetric", i, j, tmp[i,j], tmp[j, i]))
        return(FALSE)
      }
    }
  }
  
  # Get the first row and column as reference vectors
  ref_row <- sort(mat[1, ])
  ref_col <- sort(mat[, 1])
  
  # Check if all rows are permutations of the first row
  for (i in 1:nrow(mat)) {
    if (!all(abs(sort(mat[i, ]) - ref_row) <= 1e-6)) {
      print(paste("Issue with row", i))
      return(FALSE)
    }
  }
  
  # Check if all columns are permutations of the first column
  for (j in 1:ncol(mat)) {
    if (!all(abs(sort(mat[, j]) - ref_col) <= 1e-6)) {
      print(paste("Issue with col", j))
      return(FALSE)
    }
  }
  
  # If all checks passed, return TRUE
  return(TRUE)
}

