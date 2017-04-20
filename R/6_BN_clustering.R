#Compute denominator of the expectation step
logSumExp<-function(x)
{
  y = apply(x, 1, max)
  x = sweep(x, 1, y, "-")
  s = y + log(rowSums(exp(x)))
  return(s)
}

#Estimate coefficients of the parents for the maximization step
#A Linear Gaussian distribution is represented as N(x|B_0+B_1*x_p1+B_2*x_p2,sigma)
#B_0+B_1*x_p1+B_2*x_p2 is the mean of the gaussian, B_ is a coefficient 
#and x_p are the column vector of each parent of the actual node
#To compute the mle, the partial derivative of each coefficient is computed giving as result a linear system
estimateCoefficients<-function(data, responsabilities, node, parents)
{
	#Initialize coefficient and independent term matrix
	coefficient_matrix <- matrix(NA, nrow = (length(parents) + 1), ncol = (length(parents) + 1))
	independent_term_matrix <- matrix(rep(NA, length(parents) + 1), ncol = 1)
	
	#The column of the data that corresponds to the actual node is represented as a row to multiply faster
	node_column_values <- matrix(data[ ,node], nrow = 1)
  #The same as before but with the parents nodes. We added a column of 1s to multiply the intercept
	parent_columns_values <- as.matrix(cbind(1, data[ ,parents]))
  
  #Responsabilities (weights) multiplying each parent
	comb_respons_parents <- matrix(sweep(data.matrix(data[ ,parents]), 1, responsabilities, "*"),ncol = length(parents))#Se computa el vector de responsabilidades por los datos del padre que se esta estudiando.
	
  #Compute the ecuation of the partial derivative of the intercept
	coefficient_matrix[1, ] <- matrix(responsabilities, nrow = 1) %*% parent_columns_values
	independent_term_matrix[1] <- node_column_values %*% matrix(responsabilities, ncol = 1)

  
	#Compute the ecuation of the partial derivative of each one of the parents
	for(parent in 2:(length(parents)+1))
	{
		coefficient_matrix[parent, ] <- matrix(comb_respons_parents[ ,parent-1], nrow = 1) %*% parent_columns_values
		independent_term_matrix[parent] <- node_column_values %*% comb_respons_parents[ ,parent-1]
  }
  
	#Solve for the coefficients
	coefs<-solve(coefficient_matrix, independent_term_matrix)
  
	return(coefs)
}

#M-step
maximization<-function(bn, data, responsabilities)
{
	node_order <- bnlearn::node.ordering(bn)
	for(node in node_order)
	{
		parents <- bn[[node]]$parents
		children <- bn[[node]]$children

		if (length(parents) == 0) #If it has not got parents
      {		  
			mean <- (matrix(responsabilities, nrow = 1) %*% matrix(data[,node], ncol = 1)) / sum(responsabilities)
			coefs <- c("(Intercept)" = mean)
			resid <- data[ ,node] - mean
			sd <- sqrt((matrix(responsabilities, nrow = 1) %*% matrix((resid)^2, ncol = 1)) / sum(responsabilities))
			bn[[node]] <- list(coef = as.numeric(coefs), sd = as.numeric(sd))  
      
      }else{ #If there is at least one parent
        
          coefs <- estimateCoefficients(data, responsabilities, node, parents)
          resid <- data.matrix(data)[ ,node] - matrix(coefs, nrow = 1) %*% t(as.matrix(cbind(1, data[ ,parents])))
  			  sd <- sqrt((matrix(responsabilities, nrow = 1) %*% matrix((resid)^2, ncol = 1)) / sum(responsabilities))
          bn[[node]] <- list(coef = as.numeric(coefs), sd = as.numeric(sd))
       }#ELSE
	}#FOR
	return(bn)
}

#E-step
#Compute (alpha_k * N(x|theta_k))/sum^K_i(alpha_i * N(x|theta_i))
expectation <- function(model_parameters, data, priori)
{
  #Compute loglikelihood for each cluster and instance
  log_likelihood <- matrix(NA, nrow = nrow(data), ncol = length(model_parameters))
  for(k in 1:length(model_parameters))
  {
    log_likelihood[ ,k] <- logLik(model_parameters[[k]], data,by.sample = T)
  }
  
  #Compute the weights or responsabilities of each cluster in each instance with logSumExp to avoid underflow
  priori_loglikelihood <- sweep(log_likelihood, 2, log(priori), "+")
  denominator <- logSumExp(priori_loglikelihood)  
  weights <- exp(sweep(priori_loglikelihood, 1, denominator, "-"))
  
  return (list(weights = weights, log_lik = sum(denominator)))
}

#' Clustering Gaussian Bayesian network
#' 
#' Given a Bayesian network structure, apply clustering to the data with the aim to find clusters of somas
#'
#' @param bayesian_structure a bn object from bnlearn package whose nodes are gaussians
#' @param data the same dataset used to compute the bayesian_structure
#' @param clusters array of values where each value indicates the number of clusters for that model  
#' @param initialization a natural number to initialize the algorithm. Use 0 to use kmeans initialization and any other positive number as a seed for a random initialization
#' @param verbose show BIC scores for each model and the optimal number of clusters
#' 
#' @return model_list is the list of parameters that best fit the data according to BIC criteria
#' 
#' @examples 
#' somas <- read.table(file.path(tempdir(),"somaReebParameters.csv"), sep = ";", header = T, row.names = 1)
#' bayesian_structure <- BN_structure_learning(path_to_csv = file.path(tempdir(),"somaReebParameters.csv"), nboots = 200, significant_arcs = 0.95)
#' fittest_model <- BN_clustering(bayesian_structure, data = somas, clusters = c(2, 3, 4), initialization = 2, T)
#' 
#' @export
BN_clustering <- function(bayesian_structure, data, clusters, initialization, verbose)
{
  #Initializataion
  model_parameters <- list()
  model_list <- list()
  priori_list <- list()
  BIC_list <- list()
  weights_list <- list()
    
  #For each of the integers introduced by the user compute the mle for the parameters
  for(cluster in 1:length(clusters))
  {
    #Reboot variables
    old_log <- -Inf
    K <- clusters[cluster]
    data <- data.frame(data)
 
    #Estimate initial parameters. If it is 0 compute k-means if it is any other value with that seed random initialization
    if(initialization == 0 )
    {
      initiation <- kmeans(data, K)$cluster
      }else{
        set.seed(initialization)
        initiation <- round(runif(nrow(data), min = 1, max = K))
    }
  
    #Initialization of the parameters
    priori <- matrix(0, nrow = K, ncol = 1)
    for(i in 1:K)
    {
      model_parameters[[i]] <- bnlearn::bn.fit(bayesian_structure, data[which(initiation == i), ])
      priori[i] <- length(which(initiation == i)) / length(initiation)
    }
  
    e_step <- expectation(model_parameters, data, priori)
    new_log <- sum(e_step$log_lik)
    
    #When a NA value is detected execution stops. When there is not enough data to fit the parameters NA values appears 
    if(is.na(new_log))
    {
      stop("Not enough data in one of the clusters to compute the parameters. Please, change initialization values or reduce the number of clusters.")
    }
    
    #Repeat until convergence
    while(new_log > old_log)
    {
	      #M-step
	      #Maximizacion de los parametros para cada uno de los nodos
	      for(k in 1:K)
	      {
	        model_parameters[[k]] <- maximization(model_parameters[[k]], data, e_step$weights[ ,k])
	      }

	      #Priori maximization
	      priori <- colSums(e_step$weights) / nrow(e_step$weights)
        #end M-step
      
        #E-step
        e_step <- expectation(model_parameters, data, priori)
        #end E-step
      
	      old_log <- new_log
	      new_log <- e_step$log_lik
	    
        if(is.na(new_log))
	      {
	        stop(paste0("Not enough data in one of the clusters to compute the parameters.",
                      "Please, change initialization values or reduce the number of clusters."))
	      }
      
        individualBIC <- matrix(0, nrow = length(model_parameters), ncol = 1)
        for(i in 1:length(model_parameters))
	      {
	        individualBIC[i] <- BIC(model_parameters[[i]], data)
        }
    }

    individual_BIC <- matrix(0, nrow = length(model_parameters), ncol = 1)
    for(i in 1:length(model_parameters))
    {
      individual_BIC[i] <- BIC(model_parameters[[i]], data)
    }
    
    model_list[[cluster]] <- model_parameters
    priori_list[[cluster]] <- priori
    BIC_list[[cluster]] <- sum(individual_BIC)
    weights_list[[cluster]] <- e_step$weights
  }
  
  if(verbose)
  {
    print("BIC values by cluster order:")
    print(unlist(BIC_list))
    print(paste0("Optimal number of clusters = ", clusters[which.min(unlist(BIC_list))]))
  }
  return(list(priori=priori_list[[which.min(unlist(BIC_list))]], weights=weights_list[[which.min(unlist(BIC_list))]], parameters=model_list[[which.min(unlist(BIC_list))]]))
  #return(list(priori=priori_list[[which.min(unlist(BIC_list))]], parameters=model_list[[which.min(unlist(BIC_list))]]))
}
