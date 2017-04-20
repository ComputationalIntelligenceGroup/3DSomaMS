library(caret)
library(bnlearn)

#' Compute the structure of the BN
#' 
#' Compute the structure of a Dynamical Spatial Bayesian Network. Also apply bootstrap to choose the most significant arcs.
#'
#' @param somas dataframe with the characterization of the somas
#' @param n_boot number of times that bootstrap learn the bayesian structure
#' @param significant_arcs a value between 0 and 1 that denotes a threshold over percentage of times that an arc appeared in the boostrap structures
#' 
#' @return bayesian_model an object of class bn from bnlearn package that representates the structure of the DSBN
#' 
#' @export

BN_structure_learning <- function(somas, nboots, significant_arcs) 
{
  #Read csv file
  #somas <- read.table(path_to_csv, sep = ";", header = T, row.names = 1)
  
  #Column index of the variables
  height_position <- grep("height", colnames(somas))
  phi_position <- grep("^phi", colnames(somas))
  theta_position <- grep("^theta", colnames(somas))
  r_h_position <- grep("r_h", colnames(somas))
  e_position <- grep("^e",colnames(somas))
  PCA_azimuth_position <- grep("PCA_phi", colnames(somas))
  PCA_theta_position <- grep("PCA_theta", colnames(somas))
  w_position <- grep("^w", colnames(somas))
  b_0_position <- grep("b_0", colnames(somas))
  b_1_position <- grep("b_1", colnames(somas))
  b_2_position <- grep("b_2", colnames(somas))

  #Generate the blacklist of variables for the bayesian network. Actual state define the variables in a section.
  #Transition is an arc from the variables of the actual section to the next section of the soma
  #There can be arcs between the variables of the actual state but not between transitions
  #Neither can be arcs from transitions to actual state. Transitions only point out the next state
  dataset <- data.frame()
  state_names <- c("height", "phi", "theta", "r_h", "e", "PCA_phi", "PCA_theta", "w", "b_0", "b_1", "b_2")
  transition_names <- paste0(state_names, "+1")
  black_list <- expand.grid(transition_names, transition_names)
  black_list <- rbind(cbind(as.character(black_list[ ,1]), as.character(black_list[ ,2])), 
                      cbind(as.character(black_list[ ,1]), gsub("\\+1", "", as.character(black_list[ ,2]))))


  #Create a new dataset combining actual state with the next future state. Each state is a section of the soma.
  #Thus, we obtain a dynamic bayesian network applied to spatial variables.
  for(i in 1:(length(height_position) - 2))
  {
    xi <- c(height_position[i], phi_position[i], theta_position[i], r_h_position[i], e_position[i], 
            PCA_azimuth_position[i], PCA_theta_position[i], w_position[i], b_0_position[i], b_1_position[i], 
            b_2_position[i])
    
    xi1 <- c(height_position[i+1], phi_position[i+1], theta_position[i+1], r_h_position[i+1], e_position[i+1],
             PCA_azimuth_position[i+1], PCA_theta_position[i+1], w_position[i+1], b_0_position[i+1],
             b_1_position[i+1], b_2_position[i+1])	
    
    state <- somas[ ,xi]
    colnames(state) <- state_names
    transition <- somas[ ,xi1]
    colnames(transition) <- transition_names
    dataset <- rbind(dataset, cbind(state, transition))
  }
  
  #Sampling index of the dataset to apply bootstrap
  bootstrap <- list()
  for(i in 1:nboots)
  {
    bootstrap[[i]] <- sample(1:nrow(dataset), size = nrow(dataset), replace = T)
  }

  #All posible combinations of variables
  all_combinations <- expand.grid(colnames(dataset), colnames(dataset))
  all_combinations <- all_combinations[-which(all_combinations[ ,1] == all_combinations[ ,2]), ]
  distribution <- rep(0, length = nrow(all_combinations))

  for(i in 1:nboots)
  {
    model <- bnlearn::hc(dataset[bootstrap[[i]], ], blacklist = black_list)
    white_list <- as.data.frame(bnlearn::arcs(model))
    update_positions <- which(paste(all_combinations[ ,1], all_combinations[ ,2]) %in% paste(white_list[ ,1], white_list[ ,2]))
    distribution[update_positions] <- distribution[update_positions] + 1
  }

  #Arcs that were in more than the significant_arcs percentage of learned structures are included in the future model
  #Arcs that were in less than the significant_arcs percentage of learned structures are discarded
  #Thus, only the robust arcs are saved. It is useful when the number of variables is bigger than the number of instances
  white_list <- all_combinations[which(distribution > (nboots * significant_arcs)), ]
  new_black_list <- all_combinations[which(distribution <= (nboots * significant_arcs)), ]
  
  bootstrap_model <- bnlearn::hc(dataset, whitelist = white_list, blacklist = new_black_list)

  #Develop the network repeating the learned structure for each variable
  model_arcs <- data.frame()
  for(i in 1:(length(height_position) - 2))
  {
    graph_arcs <- bnlearn::arcs(bootstrap_model)
    trans_variables_pos <- grep("\\+1", graph_arcs)
    graph_arcs <- gsub("\\+1", i + 1, graph_arcs)
    graph_arcs[setdiff(1:length(graph_arcs), trans_variables_pos)] <- paste0(graph_arcs[setdiff(1:length(graph_arcs), trans_variables_pos)], i)
    model_arcs <- rbind(model_arcs, graph_arcs)
  }

  graph_arcs <- bnlearn::arcs(bootstrap_model)
  graph_arcs <- graph_arcs[which(graph_arcs[ ,2] == "height+1"), ]
  if(dim(graph_arcs)[1] > 0)
  {
    trans_variables_pos <- grep("\\+1", graph_arcs)
    graph_arcs <- gsub("\\+1", length(height_position), graph_arcs)
    graph_arcs[setdiff(1:length(graph_arcs),trans_variables_pos)] <- paste0(graph_arcs[setdiff(1:length(graph_arcs), trans_variables_pos)], length(height_position) - 1)
    model_arcs <- rbind(model_arcs, graph_arcs)
  }

  all_combinations <- expand.grid(colnames(somas), colnames(somas))
  black_list_reloaded <- setdiff(all_combinations, model_arcs)

  bayesian_model <- bnlearn::hc(somas, whitelist = model_arcs, blacklist = black_list_reloaded)

  return(bayesian_model)
}