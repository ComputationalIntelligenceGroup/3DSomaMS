#' Read VRMLs of somas 
#' 
#' Read all the VRMLs of somas in a folder and export them to PLY format in other folder
#'
#' @param read_directory path to the folder where VRML files are placed
#' @param write_directory path to the folder where PLY files will be saved    
#' 
#' @return None
#' 
#' @examples 
#' convert_somas_to_PLY(read_directory = system.file("test/VRMLs",package="SomaMS"), write_directory = file.path(tempdir(),"PLYs"), parallel = TRUE)
#' 
#' @export
convert_somas_to_PLY <- function(read_directory, write_directory, parallel = TRUE) {
    if (!file.exists(file.path(read_directory))){
      stop("The reading directory does not exist")
    }
  
    if (!file.exists(file.path(write_directory))){
      if (write_directory == "") {
        stop("The writing directory is not defined")
      }else{
        warning("The writing directory does not exist, it was created automatically")
        dir.create(write_directory, showWarning = TRUE, recursive = TRUE)
      }
    }
    
    VRML_files <- list.files(file.path(read_directory))[which(file_ext(list.files(read_directory)) == "vrml")]
    
    if (length(VRML_files)==0){
      warning('There are not VRMLs in the reading directory')
    }
    
    if(!parallel)
    {
      ncores <- 1
      cl <- NULL
    }else{
      ncores <- detectCores()
      if(is.na(ncores)){
        ncores <- 1
        cl <- NULL
      }else{
        ncores <- ncores - 1
        cl <- makeCluster(ncores, type = "SOCK")
        registerDoSNOW(cl)
      }
    }
    
    if(!parallel){
      foreach (file = VRML_files) %do% {
        # Number 1 denotes that it must be choosen that mesh with more vertices. It is useful if there are more
        # than one mesh save in the same file. p.e: 1 mesh with soma+dendrites and 1 mesh with the soma
        # segmented by an expert.
        soma <- read_VRML(file.path(read_directory, file))
        
        soma_name <- file_path_sans_ext(file)
        soma_name <- gsub(" ", "_", soma_name)
        
        vcgPlyWrite(tmesh3d(rbind(t(soma$vertices), 1), t(soma$faces)), file.path(write_directory, soma_name))
      }
    }else{
      foreach (file = VRML_files, .packages=c("Rvcg","Morpho","tools", "rgl"), .export = "read_VRML") %dopar% {
        # Number 1 denotes that it must be choosen that mesh with more vertices. It is useful if there are more
        # than one mesh save in the same file. p.e: 1 mesh with soma+dendrites and 1 mesh with the soma
        # segmented by an expert.
        soma <- read_VRML(file.path(read_directory, file))
        
        soma_name <- file_path_sans_ext(file)
        soma_name <- gsub(" ", "_", soma_name)
        
        vcgPlyWrite(tmesh3d(rbind(t(soma$vertices), 1), t(soma$faces)), file.path(write_directory, soma_name))
      }
    }
    
    if(parallel)
    {
      stopCluster(cl)
    }
    
    return(0)
} 


