#options(SomaMS.meshlabserver.path="/home/luis/package_sources/MESHLAB_133/meshlab_working/meshlab/distrib/meshlabserver")

# Get meshlabserver path
meslabserverpath <- function(){
  
  # Check for a given path
  path <- getOption("SomaMS.meshlabserver.path")
  
  # Check that the file exists and
  if(!is.null(path)){
    if( file.exists(path) )
      return(path)
    else
      warning(sprintf("Meshlab server not found in %s. Falling back to default path",path))
  }

  # Default behaviour: meshlabserver in path
  if(.Platform$OS.type == "unix") {
    return("meshlabserver")
  }
  else {
    return ("meshlabserver.exe")
  }
}

# Compute ambient occlusion through MeshLab.
execute_meshlab_script <- function(input = "", output = "", meshlab_script) {
  if (input == "") {
    input <- tk_choose.files(default = "", caption = "Select a file")
  }
  if (output == "") {
    output <- tclvalue(tkgetSaveFile())
  }
  
  # Force the file to be have .ply extension
  if (file_ext(output) != "ply") {
    output <- paste0(output, ".ply")
  }
  
  # Write the script to run with cmd. NOTE: A path with one or more spaces does not work. Example: 'C:/my
  # soma.off' does not work. 'C:/mysoma.off' is OK.
  
  # Path to MeshLab script
  if(meshlab_script=="ambient_occlusion"){
    script <- file.path(system.file("meshlab_scripts",package="SomaMS"),"ambientOcclusion.mlx")
  }else if(meshlab_script=="shape_diameter"){
    script <- file.path(system.file("meshlab_scripts",package="SomaMS"),"shapeDiameter.mlx")
  }else if(meshlab_script=="poisson_reconstruction"){
    script <- file.path(system.file("meshlab_scripts",package="SomaMS"),"poisson_reconstruction.mlx")
  }else if(meshlab_script=="poisson_reconstruction_with_normals"){
    script <- script <- file.path(system.file("meshlab_scripts",package="SomaMS"),"poisson_reconstruction_with_normals.mlx")
  }
  
  # shquote
  if(.Platform$OS.type == "unix") {
    type <- "sh"
  }
  else {
    type <- "cmd"
  }
  
  # Run the script
  #commandOrder <- paste0(meslabserverpath()," -i '", input, "' -o '", output, "' -m vq ", "-s '", script,"'")
  return( system2(meslabserverpath(),c("-i", shQuote(input,type=type),
                               "-o", shQuote(output,type=type),
                               "-m vq",
                               "-s", shQuote(script,type=type) ),
                                stderr = "" , stdout = "") )
}

