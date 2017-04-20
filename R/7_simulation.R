#' Simulation of new somas
#' 
#' Simulation of new somas of the clustering of bayesian networks
#'
#' @param model a bn.fit object from bnlearn package whose nodes are gaussians. It must be the result of the BN_clustering function
#' @param data the same dataset used to compute the bayesian_structure
#' @param cluster_number a number between 1 and the number of clusters of model. When it is NULL, new somas are generated from the priori distribution. When is a number between 1 and the number of clusters, simulate new somas from the selected cluster.
#' @param number_of_new_instances a natural number denoting the number of new somas sampled from simulation
#' @param seed a natural number used as seed during the sampling, by default it is 1
#' 
#' @return new_somas a dataset similar to data with the simulated somas
#' 
#' 
#' @export
simulate_new_somas <- function(model, data, cluster_number = NULL, number_of_new_instances = 100, seed = 1)
{
  new_somas<-data.frame()
  if(is.null(cluster_number)) #New somas from priori distribution
  {
    set.seed(seed)
    sampling_cluster <- sample(length(model$priori), number_of_new_instances, replace = T, prob = model$priori)#Sample clusters according to priori distribution
    soma_classes <- as.numeric(table(sampling_cluster)) #Count the number of instances for each cluster
    
    #Simulate new somas from each cluster according to the number of instances computed in the previous sampling
    for(i in 1:length(soma_classes))
    {
      
      new_somas <- rbind(new_somas, cluster_simulation(model, data, i, soma_classes[i], seed))
    
    } 
 
    }else{ #Sample from the selected cluster
    if(cluster_number %in% c(1:length(model$priori)))
    {
      
      new_somas <- cluster_simulation(model, data, cluster_number, number_of_new_instances, seed)
      
      }else{ 
        
        stop(paste0("The selected cluster number does not exist. The model has ", length(model), " clusters"))
      
      }
  }
  
  return(new_somas)
}


#' Simulation of new somas from a selected cluster
#' 
#' Simulation of new somas from the choosen cluster of the clustering of bayesian networks
#'
#' @param model a bn.fit object from bnlearn package whose nodes are gaussians. It must be the result of the BN_clustering function
#' @param data the same dataset used to compute the bayesian_structure
#' @param cluster_number a natural number between 1 and the number of clusters of model. Simulate new somas from the selected cluster.
#' @param number_of_new_instances a natural number denoting the number of new somas sampled from simulation
#' @param seed a natural number used as seed during the sampling, by default it is 1
#' 
#' @return new_somas a dataset similar to data with the simulated somas of the selected cluster
cluster_simulation <- function(model, data, cluster_number, number_of_new_instances, seed = 1)
{
  set.seed(seed)
  
  parameters <- model$parameters[[cluster_number]]
  nodeOrder <- node.ordering(parameters)
  
  new_somas <- data.frame(matrix(NA, nrow = number_of_new_instances, ncol = length(nodeOrder)))
  colnames(new_somas) <- names(parameters)

  classes <- apply(model$weights, 1, which.max)
  maximum <- apply(data[which(classes == cluster_number), nodeOrder], 2, max)
  minimum <- apply(data[which(classes == cluster_number), nodeOrder], 2, min)
  
  for(i in 1:length(nodeOrder))
  {
    node.parameters <- parameters[[nodeOrder[i]]]
    node.parents <- parameters[[nodeOrder[i]]]$parents
    
    mean <- as.numeric(node.parameters$coefficients %*% rbind(1, t(data.matrix(new_somas[ ,node.parents]))))
    sd <- node.parameters$sd
    new_somas[ ,nodeOrder[i]] <- rnorm(number_of_new_instances, mean, sd)
  }
  return(new_somas)
}


#' Three-dimensional representation of the simulated somas
#' 
#' Three-dimensional renderization from a dataset where the simulated somas were saved
#'
#' @param new_somas a dataset with the simulated somas obtained in simulate_new_somas function
#' @param index index of the row of the new soma to renderize
#' 
#' @return curve_points a n x 3 dataframe where each row represents a point in the cartesian three-dimensional space
#' 
#' @export
simulation_3D_mesh<-function(new_somas, index)
{   
  
  curve_points<-data.frame(x = 0, y = 0, z = 0)
  cartesianEllipses<-list()
  cartesianEllipses[[1]]<-c(0,0,0)
  matrixOfFaces <- data.frame()
  
	#Column index of the variables
  height_position <- grep("height", colnames(new_somas))
  phi_position <- grep("^phi", colnames(new_somas))
	theta_position <- grep("^theta", colnames(new_somas))
	r_h_position <- grep("r_h", colnames(new_somas))
	e_position <- grep("^e",colnames(new_somas))
	PCA_azimuth_position <- grep("PCA_phi", colnames(new_somas))
	PCA_theta_position <- grep("PCA_theta", colnames(new_somas))
	w_position <- grep("^w", colnames(new_somas))
	b_position <- grep("b_", colnames(new_somas))
	
	vector_length <- as.numeric(new_somas[index, height_position])
	phi <- as.numeric(c(0, new_somas[index, phi_position]))
	theta <- as.numeric(c(pi/2, new_somas[index, theta_position]))
	r_h <- as.numeric(new_somas[index, r_h_position])
	eccentricity <- as.numeric(new_somas[index, e_position])
	perp_vectors_sph <- matrix(as.numeric(new_somas[index, c(PCA_azimuth_position, PCA_theta_position)]), ncol = 2)
	ellipse_angles <- as.numeric(new_somas[index, w_position])
	curve_coefficients <- matrix(as.numeric(new_somas[index, c(b_position)]), ncol = 3)
	
  matrix_skeleton <- matrix(0, nrow = 7, ncol = 3)
  matrix_skeleton <- spherical2Cartesian(phi, theta, vector_length)
  
  skeleton <- apply(matrix_skeleton, 2, cumsum)

  major_axis <- r_h * (vector_length[1:(length(vector_length)-1)] + vector_length[2:length(vector_length)])
  minor_axis <- major_axis * eccentricity
  
  perp_PCA_matrix <- spherical2Cartesian(perp_vectors_sph[ ,1], perp_vectors_sph[ ,2], 1)
    
  for (i in 1:nrow(perp_PCA_matrix))
  {
    u <- vcrossp(perp_PCA_matrix[i, ], c(0, 0, 1))
    u <- u / as.numeric(normv(u))
    rotation_angle <- acos((perp_PCA_matrix[i,] %*% c(0, 0, 1)))
    rotationMatrix <- rotation_matrix(u,rotation_angle)
  
    centers <- rotationMatrix %*% matrix(skeleton[i, ], ncol = 1)
    
    fit <- list()
    fit$angle <- ellipse_angles[i]
    fit$maj <- major_axis[i]
    fit$min <- minor_axis[i]
    fit$center <- centers[1:2]
    
    points_ellipse <- get.ellipse(fit)
    x <- points_ellipse[ ,1]
    y <- points_ellipse[ ,2]
    z <- centers[3]
    
    elipse <- cbind(x, y, z)
#     u <- vcrossp(c(cos(ellipse_angles[i]), sin(ellipse_angles[i]), 0), c(1, 0, 0));
#     u <- u / as.numeric(norm(u))
#     rotationZ <- rotation_matrix(u, ellipse_angles[i])
#   
#     rotated_ellipse <- cbind(x, y, matrix(0, nrow = length(x), ncol = 1))
#     rotated_ellipse <- t(rotationZ %*% t(rotated_ellipse))
#  
#     z <- curve_coefficients[i, 1] + curve_coefficients[i,2] * (rotated_ellipse[ ,1]^2) + curve_coefficients[i,3] * (rotated_ellipse[,2]^2)
#  
#     rotated_ellipse[,3] <- z
#   
#     elipse <- t(t(rotationZ) %*% t(rotated_ellipse))
#     
    recovered_ellipse <- t(t(rotationMatrix) %*% t(elipse))
    curve_points <- rbind(curve_points, data.frame(x = recovered_ellipse[,1], y = recovered_ellipse[,2], z = recovered_ellipse[,3]))
    cartesianEllipses[[i+1]]<-recovered_ellipse
  }

  curve_points <- rbind(curve_points,skeleton[nrow(skeleton),])
  
  
  points<-rbind(c(0,0,0),cartesianEllipses[[2]])
  faces<-cbind(1,seq(2,(nrow(points)-1),by=1),seq(3,nrow(points),by=1))
  matrixOfFaces<-rbind(matrixOfFaces,faces)

  for(i in 2:(length(cartesianEllipses)-2))
  {
    points<-matrix(rbind(cartesianEllipses[[i]],cartesianEllipses[[i+1]]),ncol=3,byrow=F)

    initialPosition<-which.min(cartesianEllipses[[i]][,3])
    
    apex<-seq(initialPosition,length=nrow(cartesianEllipses[[i]])/2,by=2)
    apex[apex>nrow(cartesianEllipses[[i]])]<-(apex[apex>nrow(cartesianEllipses[[i]])]%%(nrow(cartesianEllipses[[i]])+1))+1
    
    positionMatch<-which.min(sqrt(rowSums(sweep(cartesianEllipses[[i+1]],2,cartesianEllipses[[i]][initialPosition,],"-")^2)))+nrow(cartesianEllipses[[i]])
    base<-seq(positionMatch,length=nrow(cartesianEllipses[[i]])/2,by=2)
    
    faces1<-cbind(apex,base-1,base)
    faces1[faces1>nrow(points)]<-(faces1[faces1>nrow(points)]%%(nrow(points)+1))+361
    
    faces2<-cbind(apex,base,base+1)
    faces2[faces2>nrow(points)]<-(faces2[faces2>nrow(points)]%%(nrow(points)+1))+361
    
    positionMatch<-positionMatch+1
    initialPosition<-initialPosition+1
    apex<-seq(positionMatch,length=nrow(cartesianEllipses[[i]])/2,by=2)
    apex[apex>nrow(points)]<-(apex[apex>nrow(points)]%%(nrow(points)+1))+361
    
    base<-seq(initialPosition,length=nrow(cartesianEllipses[[i]])/2,by=2)
    
    prefaces<-cbind(base-1,base)
    prefaces[prefaces>nrow(cartesianEllipses[[i]])]<-(prefaces[prefaces>nrow(cartesianEllipses[[i]])]%%(nrow(cartesianEllipses[[i]])+1))+1
    faces3<-cbind(apex,prefaces)
    
    
    prefaces<-cbind(base+1,base)
    prefaces[prefaces>nrow(cartesianEllipses[[i]])]<-(prefaces[prefaces>nrow(cartesianEllipses[[i]])]%%(nrow(cartesianEllipses[[i]])+1))+1
    faces4<-cbind(apex,prefaces)
    
    
    faces<-rbind(faces1,faces2,faces3,faces4)
    colnames(faces)<-NULL
    matrixOfFaces<-rbind(matrixOfFaces,(faces+1+(360*(i-2))))
  }
  
  points<-rbind(cartesianEllipses[[7]],skeleton[nrow(skeleton),])
  faces<-cbind(seq(1,(nrow(points)-1),by=1),seq(2,nrow(points),by=1),361)
  matrixOfFaces<-rbind(matrixOfFaces,(faces+1+(360*(i-2+1))))
  
  points<-rbind(skeleton[1,], cartesianEllipses[[2]])
  faces<-cbind(seq(1,(nrow(points)-1),by=1),seq(2,nrow(points),by=1),361)
  matrixOfFaces<-rbind(matrixOfFaces,(faces+1+(360*(i-2+1))))
  
  reoriented <- vcgClean(tmesh3d( t(cbind(curve_points,1)),t(matrixOfFaces) ),sel=7)
  new_vertices <- vcgSample(reoriented,type="mc",SampleNum = 1E4)
  
  write_PLY_only_points( new_vertices, file.path( tempdir(), "vertices_non_poisson" ))
  execute_meshlab_script(paste0(file.path( tempdir(), "vertices_non_poisson" ), ".ply"), paste0(file.path( tempdir(), "vertices_poisson" ), ".ply"), "poisson_reconstruction_with_normals")
  poisson_output_mesh <- vcgPlyRead(paste0(file.path( tempdir(), "vertices_poisson" ), ".ply"))

  return(list(vertices = t(poisson_output_mesh$vb[1:3,]), faces = t(poisson_output_mesh$it) ))
}