library(Rvcg)
library(Morpho)
library(rgl)
library(misc3d)
library(geometry)
library(pracma)
library(matrixStats)
library(tools)
library(bnlearn)

source("./R/Geometric_operators.R")

#' Characterize soma with ray tracing reeb graph 
#' 
#' Compute the morphological characterization of the soma. First, cast rays from the centroid of the soma to the surface,
#' defining level curves. Next, orientate soma and adjust curves. Finally compute parameters according to the curves.
#'
#' @param path_somas_repaired path to the folder where segmented somas were placed. They have to be PLY files.
#' @param path_csv_file path to a csv file where the values of the characterization will be saved 
#' 
#' @return None
#' @examples 
#' soma_caracterization(path_somas_repaired = file.path(tempdir(),"reconstructed"), path_csv_file = file.path(tempdir(),"somaReebParameters.csv"))
#' 
#' @export

soma_caracterization <- function(path_somas_repaired, path_csv_file)
{
  PLY_files <- list.files(path_somas_repaired)[which(file_ext(list.files(path_somas_repaired)) == "ply")]
  soma_names <- unlist(strsplit(basename(PLY_files), "\\."))
  soma_names <- soma_names[which(soma_names != "ply")]
  dataset <- data.frame()
  
  #For each soma file compute parameters and save in the csv file
  for(file in PLY_files)
  {
    #Read soma from PLY file
    soma_PLY <- vcgImport(paste(path_somas_repaired, file, sep = "/"), F, F, F)
    
    #Remove isolate pieces and smooth the mesh
    smoothed_mesh <- vcgSmooth(vcgIsolated(soma_PLY, facenum = NULL, diameter = NULL), type = c("laplacian"), iteration = 50)
    
    #Compute centroid and place the centroid in the origin
    center_of_mass <- colSums(vert2points(smoothed_mesh)) / nrow(vert2points(smoothed_mesh))
    smoothed_mesh <- translate3d(smoothed_mesh, -center_of_mass[1], -center_of_mass[2], -center_of_mass[3])
    
    vertices_soma <- vert2points(smoothed_mesh)
    
    #Compute Feret diameter to obtain the 
    convex_hull <- convhulln(vertices_soma, "FA")$hull
    convex_vertices <- vertices_soma[sort(unique(c(convex_hull))), ]
    radius <- max(apply(convex_vertices, 1, normv))

    #Generate points in the sphere to cast rays to them. Each curve is composed by 3000 points (theta).
    #We generate 6 curves and 2 extreme points (phi). 
    set.seed(1)
    theta <- seq(0, 2*pi, length = 3000)
    phi <- seq(pi, 0, length = 8)
    
    #Normal vectors denote the direction of the rays.
    malla <- smoothed_mesh
    malla$vb[1:3, ] <- t(vertices_soma)
    ellipse <- list()
    list_of_ellipses <- list()
    ellipses_centroids <- matrix(NA, ncol = 3, nrow = length(phi))
    for(i in 1:length(phi))
    {
      #Cartesian points from the spherical coordinates
      x <- radius * sin(theta) * sin(phi[i])
      y <- radius * cos(phi[i])
      z <- radius * cos(theta) * sin(phi[i])
      
      #Cartesian points denote the normals
      normals <- cbind(x,y,z)
      
      #Define origin for each ray. As always is the same point, we generate a zero matrix
      vertices <- matrix(rep(c(0, 0, 0), 3000), ncol = 3000)
      
      x$vb <- vertices
      x$normals <- t(normals)
      #Change class to mesh3d to execute vcgRaySearch
      class(x) <- "mesh3d"
      
      #Cast rays
      raytracer <- vcgRaySearch(x,malla)
      ellipse[[i]] <- raytracer
      #Save each vertex of the curve
      list_of_ellipses[[i]] <- vert2points(raytracer)
      #Centroid of each curve
      ellipses_centroids[i, ] <- colSums(list_of_ellipses[[i]]) / nrow(list_of_ellipses[[i]])
    }
    
    #Place lowest ellipse in the origin
    insertion_point<-ellipses_centroids[1, ]
    for(i in 1:length(list_of_ellipses))
    {
      list_of_ellipses[[i]] <- translate3d(list_of_ellipses[[i]], -insertion_point[1], -insertion_point[2], -insertion_point[3])  
    }
    
 
    #Also shift the soma mesh to the same place
    smoothed_mesh <- translate3d(smoothed_mesh, -insertion_point[1], -insertion_point[2], -insertion_point[3])                             
    ellipses_centroids <- translate3d(ellipses_centroids, -insertion_point[1], -insertion_point[2], -insertion_point[3]) 
    
    number_of_curves<-length(list_of_ellipses) - 2
    
    perp_vectors_sph <- matrix(0, nrow = number_of_curves, ncol = 2) #Perpendicular vector to curve in spherical coords.
    ellipse_angles <- matrix(0, nrow = number_of_curves, ncol = 1)     #Angle of the elipse in the XY plane.     
    eccentricity <- matrix(0, nrow = number_of_curves, ncol = 1)      #Relation between semiaxes.
    curve_coefficients <- matrix(0, number_of_curves, ncol = 3) #Coefficients of z=a_0+a_1*x^2+a_2*y^2.
    curve_centroids <- matrix(0, nrow = length(list_of_ellipses), ncol = 3)      #Centroids of the curves.
    perp_PCA_matrix <- matrix(0, nrow = number_of_curves, ncol = 3)    #Perpendicular vector to curve in cartesian coords.
    ellipse_angles_vectors <- matrix(0, nrow = number_of_curves, ncol = 3) #Vector representation of the angle of the ellipse.
    semiaxis <- matrix(0, nrow = number_of_curves, ncol = 1)           #One of the semiaxis of the ellipse.
                                  
    #Compute parameters of each curve.
    for (i in 2:(length(list_of_ellipses) - 1))
    {
      PCA <- prcomp(list_of_ellipses[[i]]) #Compute main axis of the curve
      perp_PCA_matrix[i-1, ] <- PCA$rotation[ ,3]#Perpendicular vector to the fittest plane to the curve
          
      #Rotate curve to place the normal to the plane parallel to the Z axis, that is, the curve will be in the XY plane
      u <- vcrossp(perp_PCA_matrix[i-1, ],c(0, 0, 1))
      u <- u / as.numeric(normv(u))
      theta <- acos((perp_PCA_matrix[i-1,] %*% c(0, 0, 1))) 
      rotationMatrix <- rotation_matrix(u, theta)
                                  
      rotated_curve <- t(rotationMatrix %*% t(list_of_ellipses[[i]]))
      
      #Check that the ellipse is over the XY plane because in the other case (when it is under the plane) 
      #rotation would be in the opposite direction
      if(min(rotated_curve[ ,3]) < 0) 
      {
        perp_PCA_matrix[i-1, ] <- -perp_PCA_matrix[i-1, ]
        u <- vcrossp(perp_PCA_matrix[i-1, ],c(0, 0, 1))
        u <- u / as.numeric(normv(u))
        theta <- acos((perp_PCA_matrix[i-1, ] %*% c(0, 0, 1)))
        rotationMatrix <- rotation_matrix(u,theta)
                   
        rotated_curve <- t(rotationMatrix %*% t(list_of_ellipses[[i]]))
      }
    
      #Compute ellipse approximation to the curve
      fit_ellipse <- fit.ellipse(rotated_curve[ ,1], rotated_curve[ ,2])
    
      #Rotate angle and centroid of the ellipse to place it in the original space (before PCA rotation)
      ellipse_angles_vectors[i-1, ] <- as.numeric(t(rotationMatrix) %*% as.numeric(c(cos(fit_ellipse$angle), sin(fit_ellipse$angle), 0)))
      ellipses_centroids[i, ] <-  as.numeric(t(rotationMatrix) %*% as.numeric(c(fit_ellipse$center, mean(rotated_curve[ ,3]))))
    
      eccentricity[i-1] <- fit_ellipse$minor / fit_ellipse$major
      semiaxis[i-1] <- fit_ellipse$major
      ellipse_angles[i-1] <- fit_ellipse$angle
    
      #Place the angle of the ellipse parallel to X axis. This is needed to compute the approximation of the curve.
      #The angle of the ellipse must be always the same.
      u <- vcrossp(c(cos(fit_ellipse$angle), sin(fit_ellipse$angle), 0), c(1, 0, 0))
      u <- u / as.numeric(normv(u))
      rotationZ <- rotation_matrix(u,fit_ellipse$angle)
      rotated_curve <- t(rotationZ %*% t(rotated_curve))
                      
      x <- rotated_curve[ ,1]
      y <- rotated_curve[ ,2]
      z <- rotated_curve[ ,3]
      
      #Adjust curve with a cuadratic regression
      data <- data.frame(x, y, z)
      fitModel <- lm(z ~ I(x^2) + I(y^2), data = data)
      curve_coefficients[i-1, ] <- fitModel$coefficients
    }

    ellipses_centroids[length(list_of_ellipses), ] <- list_of_ellipses[[length(list_of_ellipses)]][1, ]

    #Rotate the soma to align the vector between de lowest point and
    #the centroid of the first curve with the Z axis. Also, rotate the
    #curves and the centroids.
    #u=rotation axis, theta=rotation angle
    unit_vector <- ellipses_centroids[2, ] / normv(ellipses_centroids[2, ]);
    u <- vcrossp(unit_vector,c(0, 0, 1));
    u <- u / as.numeric(normv(u));
    theta <- acos((unit_vector %*% c(0, 0, 1)));  

    rotationMatrix <- rotation_matrix(u, theta);

    #Rotation of all the geometric measures
    for (i in 1:length(list_of_ellipses))
    {
      list_of_ellipses[[i]] <- t(rotationMatrix %*% t(list_of_ellipses[[i]]))
    }
    ellipses_centroids <- t(rotationMatrix %*% t(ellipses_centroids))
    perp_PCA_matrix <- t(rotationMatrix %*% t(perp_PCA_matrix))
    ellipse_angles_vectors <- t(rotationMatrix %*% t(ellipse_angles_vectors))

    #Rotate the soma to place the upper point of the soma parallel to the
    #X axis. It is rotated around the Z axis.
    unit_vector <- c(ellipses_centroids[8, 1], ellipses_centroids[8, 2], 0) / as.numeric(normv(c(ellipses_centroids[8, 1], ellipses_centroids[8, 2], 0)))

    u <- vcrossp(unit_vector,c(1, 0, 0))
    u <- u / as.numeric(normv(u))
    theta <- acos((unit_vector %*% c(1, 0, 0)))

    rotationZ <- rotation_matrix(u, theta)

    #Rotate all the geometric measures
    for (i in 1:length(list_of_ellipses))
    {
      
      list_of_ellipses[[i]] <- t(rotationZ %*% t(list_of_ellipses[[i]]))
    }
    
    ellipses_centroids <- t(rotationZ %*% t(ellipses_centroids))
    perp_PCA_matrix <- t(rotationZ %*% t(perp_PCA_matrix))
    ellipse_angles_vectors <- t(rotationZ %*% t(ellipse_angles_vectors))

    #Compute current angle of the ellipses after all the transformations
    for (i in 1:nrow(perp_PCA_matrix))
    {
      u <- vcrossp(perp_PCA_matrix[i, ], c(0, 0, 1))
      u <- u / as.numeric(normv(u))
      theta <- acos((perp_PCA_matrix[i,] %*% c(0, 0, 1)))
      rotationMatrix <- rotation_matrix(u, theta)
      rotated_angle_matrix <- t(rotationMatrix %*% matrix(ellipse_angles_vectors[i, ], ncol = 1))
      ellipse_angles[i] <- cart2pol(rotated_angle_matrix[1], rotated_angle_matrix[2])[1]
    }
    
    perp_vectors_sph <- cartesian2Spherical(perp_PCA_matrix[ ,1], perp_PCA_matrix[ ,2], perp_PCA_matrix[ ,3])[ ,1:2]

    #Compute the parameters needed to build the skeleton of the spine.
    #|h|, phi, theta, alpha 
    vector_length <- matrix(0, nrow = nrow(ellipses_centroids) - 1, ncol = 1)
    phi <- matrix(0, nrow = nrow(ellipses_centroids) - 2, ncol = 1)
    theta <- matrix(0, nrow = nrow(ellipses_centroids) - 2, ncol = 1)
  
    section_vector <- diff(ellipses_centroids)
    vector_length <- sqrt(rowSums(section_vector^2))
    angles_section <- cartesian2Spherical(section_vector[ ,1], section_vector[ ,2], section_vector[ ,3])
    phi <- angles_section[ ,1]
    theta <- angles_section[ ,2]

    r_h <- matrix(0, nrow = length(list_of_ellipses) - 2, ncol = 1)
    ellipse_area <- matrix(0, nrow = length(list_of_ellipses) - 2, ncol = 1)
    r_h <- semiaxis / (vector_length[1:(length(vector_length) - 1)] + vector_length[2:(length(vector_length))])

    instance_row <- data.frame(t(c(vector_length, phi[2:length(phi)], theta[2:length(theta)], r_h,eccentricity, perp_vectors_sph, ellipse_angles, curve_coefficients)))
    rownames(instance_row)<-file
    dataset<-rbind(dataset,instance_row)
    
  }
    
    if(!file.exists(path_csv_file))
    {
      #Names of the height variables
      height_names <- list()
      for (i in 1:length(vector_length))
      {
        height_names[[i]] <- paste0("height", i)
      }
      
      #Names of all the other variables
      variable_names <- c("phi", "theta", "r_h", "e", "PCA_phi", "PCA_theta", "w", "b_0", "b_1", "b_2")
      variable_names_matrix <- matrix(rep(unlist(variable_names), number_of_curves),ncol = number_of_curves)
      for(i in 1:length(variable_names))
      {
          for(j in 1:number_of_curves)
          {
            variable_names_matrix[i,j] <- paste0(variable_names_matrix[i,j], j)
          }
      }
      column_names <- c(unlist(height_names), as.character(t(variable_names_matrix)))
      colnames(dataset) <- column_names
      
      #Save dataset
      write.table(dataset, file = path_csv_file, sep = ";", col.names = NA)
    }else{
      #If the file exist, append the new values at the end of the file
      write.table(dataset, file = path_csv_file, sep = ";", append = T, col.names = F, row.names = T)
  }
}