#' Segment somas
#' 
#' Apply Gaussian mixture models to segment soma from dendrites
#'
#' @param directory path to the folder where PLY files of the somas are saved
#' @param output_shape_diameter path to the folder where somas with shape diameter will be saved   
#' @param broken_mesh path to the folder where somas will be saved after that the GMM threshold removes the dendrites
#' @param output_poisson_reconstruction path to the folder where somas will be saved after the mesh closing process
#' @param final_result path to the folder where final segmented somas will be saved after the isolated pieces are removed
#' @return None
#' 
#' @examples  
#' segment_somas(directory = file.path(tempdir(), "poisson_reconstruction_ao"), output_shape_diameter = file.path(tempdir(), "sdf"), broken_mesh = file.path(tempdir(), "broken_mesh_sdf"), output_poisson_reconstruction = file.path(tempdir(), "poisson_reconstruction_sdf"), final_result = file.path(tempdir(), "final_result"))
#' 
#' @export
segment_somas <- function(directory, output_shape_diameter, broken_mesh, output_poisson_reconstruction, final_result) {
  if (!file.exists(file.path(directory))){
    stop("The reading directory does not exist")
  }
  
  if (!file.exists(file.path(output_shape_diameter))){
    if (output_shape_diameter == "") {
      stop("The shape diameter funcion directory is not defined")
    }else{
      warning("The shape diameter function directory does not exist, it was created automatically")
      dir.create(output_shape_diameter, showWarning = TRUE, recursive = TRUE)
    }
  }
  
  if (!file.exists(file.path(broken_mesh))){
    if (broken_mesh == "") {
      stop("The broken mesh directory is not defined")
    }else{
      warning("The broken mesh directory does not exist, it was created automatically")
      dir.create(broken_mesh, showWarning = TRUE, recursive = TRUE)
    }
  }
  
  if (!file.exists(file.path(output_poisson_reconstruction))){
    if (output_poisson_reconstruction == "") {
      stop("The poisson reconstruction directory is not defined")
    }else{
      warning("The poisson reconstruction directory does not exist, it was created automatically")
      dir.create(output_poisson_reconstruction, showWarning = TRUE, recursive = TRUE)
    }
  }  
  
  if (!file.exists(file.path(final_result))){
    if (final_result == "") {
      stop("The final result directory is not defined")
    }else{
      warning("The final result directory does not exist, it was created automatically")
      dir.create(final_result, showWarning = TRUE, recursive = TRUE)
    }
  }  
  
  PLY_files <- list.files(directory)[which(file_ext(list.files(directory)) == "ply")]
    
  for (file in PLY_files) {
    soma_name <- file_path_sans_ext(file)
        
    output <- execute_meshlab_script(input = file.path(directory, file), output = file.path(output_shape_diameter, soma_name), meshlab_script = "shape_diameter")
        
    soma <- read_binary_PLY(output)

    #First GMM
    quality_classification <- Mclust(soma$quality, G = 2)
    soma_class <- which.min(c(mean(c(apply(soma$vertices[quality_classification$classification == 1, ], 2, sd))), mean(c(apply(soma$vertices[quality_classification$classification == 2, ], 2, sd)))))
        
    quality_sec_gauss  <- soma$quality [quality_classification$classification == soma_class]
    vertices_sec_gauss <- soma$vertices[quality_classification$classification == soma_class, ]
    
    #Second GMM
    quality_classification <- Mclust(quality_sec_gauss, G = 2)
    write_PLY_only_points(vertices_sec_gauss[quality_classification$classification == which.max(table(quality_classification$classification)), ], paste0(file.path(broken_mesh, soma_name), ".ply"))

    execute_meshlab_script(paste0(file.path(broken_mesh,soma_name), ".ply"), paste0(file.path(output_poisson_reconstruction, soma_name), ".ply"), "poisson_reconstruction_with_normals")

    poisson_output_mesh <- vcgPlyRead(paste0(file.path(output_poisson_reconstruction, soma_name), ".ply"))
    
    #Poisson reconstruction sometimes generates isolated mesh pieces that must be removed
    cleaned_mesh <- vcgIsolated(poisson_output_mesh, silent=TRUE)
          
    vcgPlyWrite(cleaned_mesh, file.path(final_result, soma_name))
  }
}


#' Extract single soma from data
#' 
#' Apply Gaussian mixture models to segment soma from dendrites
#'
#' @param soma A filepath or the structure returned by read_VRML or read_PLY representing damaged soma mesh @seealso \code{read_VRML}
#' @param sdf.path folder where intermediate files with shape diameter function data are saved
#' @param broken_sdf.path Folder where intermediate files with broken somas after GMM threshold are saved
#' @param reconstruction.path Folder where reconstructed file is stored
#' @param name Name of the soma file
#' 
#' @return Extracted soma 
#' 
#' @useDynLib SomaMS
#' 
#' @export
segment.soma3d <- function(soma, sdf.path = tempdir(), broken_sdf.path = tempdir(),
                           reconstruction.path = tempdir(), name="soma"){
  
  ###
  # INPUT SANITIZATION
  ###
  # Validate repaired
  stopifnot( length(name)>0 )
  
  # Validate tempdirs
  stopifnot( utils.validateDir(sdf.path, pre="Shape Diameter Function measures", createDir = T) )
  stopifnot( utils.validateDir(broken_sdf.path, pre="Broken SDF Meshes", createDir = T) )
  stopifnot( utils.validateDir(reconstruction.path, pre="Reconstructed Meshes", createDir = T) )
  
  # If soma is already a ply path keeps it, if its a vrml read and validates  and if is a list validates.
  soma.path <- utils.somaToPLY(soma, paste("original",name,sep="_"), moveToLocation = F)
  stopifnot(!is.null(soma.path))
  
  ####
  
  ###
  #  CALL SDF
  ###
  soma.sdf.path <- file.path(sdf.path,paste0("sdf_",name,".ply"))
  if( execute_meshlab_script(input = soma.path, 
                             output = soma.sdf.path, 
                             meshlab_script = "shape_diameter") != 0)
    stop("Failure while executing SDF meshlab script")
  # Read soma.sdf
  soma.sdf <- read_binary_PLY(soma.sdf.path)
  
  #First GMM
  quality_classification <- Mclust(soma.sdf$quality, G = 2)
  soma_class <- which.min(c(mean(c(apply(soma.sdf$vertices[quality_classification$classification == 1, ], 2, sd))), 
                            mean(c(apply(soma.sdf$vertices[quality_classification$classification == 2, ], 2, sd)))))
  
  quality_sec_gauss  <- soma.sdf$quality [quality_classification$classification == soma_class]
  vertices_sec_gauss <- soma.sdf$vertices[quality_classification$classification == soma_class, ]
  
  #Second GMM
  quality_classification <- Mclust(quality_sec_gauss, G = 2)
  soma.broken_sdf.path <- file.path(broken_sdf.path,paste0("brokensdf_",name))
  
  write_PLY_only_points(
    vertices_sec_gauss[quality_classification$classification == which.max(table(quality_classification$classification)), ], 
    soma.broken_sdf.path)
  soma.broken_sdf.path <- paste0(soma.broken_sdf.path,".ply")
  
  ###
  #  CALL POISSON RECONSTRUCTION
  ###
  soma.reconstruction.path <- file.path(broken_sdf.path,paste0("reconstructed_",name,".ply"))
  if( execute_meshlab_script(input = soma.broken_sdf.path, 
                         output = soma.reconstruction.path, 
                         meshlab_script = "poisson_reconstruction_with_normals") != 0 )
    stop("Failure while executing SDF meshlab script")
  
  #Poisson reconstruction sometimes generates isolated mesh pieces that must be removed
  poisson_output_mesh <- vcgPlyRead(soma.reconstruction.path)
  poisson_output_mesh <- vcgIsolated(poisson_output_mesh, silent=TRUE)
  
  # Transform to soma3d
  return(list( vertices= t(poisson_output_mesh$vb[1:3,]), faces = t(poisson_output_mesh$it),
               normals = t(poisson_output_mesh$normals[1:3,]) ))
}

