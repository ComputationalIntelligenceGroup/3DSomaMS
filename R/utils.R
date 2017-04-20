utils.validateSoma <- function(soma, check.quality = F){
  
  # This can be done shorter in one expresison, but legibility is a thing.
  # 
  
  # Check soma null values
  if( is.null(soma)  ||  is.null(soma$vertices) || is.null(soma$faces) )
    return (FALSE)
  # Check sublist elements
  else if( !is.matrix(soma$vertices) ||  !is.matrix(soma$faces) || ncol(soma$vertices)!=3 || ncol(soma$faces)!= 3  )
    return (FALSE)
  # Check that values are numeric
  else if( !is.numeric(soma$vertices) || !is.numeric(soma$faces) )
    return (FALSE)
  # Check quality parameter
  else if( check.quality && ( is.null(soma$quality) || !is.numeric(soma$quality) ) )
    return (FALSE)
  else
    return (TRUE)
}

utils.validateDir <- function(path, pre="", createDir = F ){
  
  # Check null and empty
  if( length(path)==0 )
    return (FALSE)
  
  # Check that given paths exist (tmpdir always exist and is unique per R session)
  if (!dir.exists(path) ) {
    if(createDir){
      warning( sprintf("%s Path \"%s\" does not exist. Creating empty folder",pre,path ) )
      # Stop if failure creating ao.path
      return(dir.create(path, showWarnings = T, recursive = T ))
    }
    else{
      return(FALSE)
    }
  }
  else{
    return(TRUE)
  } 
}

utils.somaToPLY <- function(soma, filename, dir=tempdir(), moveToLocation = F ){
  
  # Remove filename extension (if it is ply)
  filename <- sub("\\.ply$","", filename, ignore.case = T)
  
  # Is already a file -> Check if it is ply
  if( is.character(soma) ){
    # Is ply
    if( tolower(file_ext(soma)) == "ply"){
      if(moveToLocation){
        
        if(file.rename(soma,file.path(dir,paste0(filename,".ply"))))
          return (file.path(dir,paste0(filename,".ply")))
        else
          stop("Error moving file to destination")
      }
      else
        return ( soma )
    }
    # Is vrml
    else if( tolower(file_ext(soma)) == "vrml"){
      soma <- read_VRML(soma)
    }
    else
      return (NULL)
  }
  
  # Is a soma (or we have used read_VRML)
  if( is.list(soma) ){
    
    stopifnot(utils.validateSoma(soma, check.quality  = F))
    
    # Dump to temp PLY file
    if( vcgPlyWrite(tmesh3d(rbind(t(soma$vertices), 1), t(soma$faces)), binary = T, file.path(dir,filename)  ) != 0)
      stop("Error while writing PLY temporal file")
    
    ## update name name (.ply is added by vcgPlyWrite)
    return( file.path(dir,paste0(filename,".ply")) )
  }
  else
    return(NULL)
}

utils.loadsoma <- function(soma){
  
  # Is already a file -> Check if it is ply
  if( is.character(soma) ){
    # Is ply
    if( tolower(file_ext(soma)) == "ply"){
      
       vcgmesh <- vcgImport(soma)
       return(list(vertices=t(vcgmesh$vb[1:3,]), faces=t(vcgmesh$it)))
    }
    # Is vrml
    else if( tolower(file_ext(soma)) == "vrml"){
      return (read_VRML(soma))
    }
    else
      return (NULL)
  }
  
  # Is a soma (or we have used read_VRML)
  if( is.list(soma) ){
    
    stopifnot(utils.validateSoma(soma, check.quality  = F))
    return(soma)
  }
  else
    return(NULL)
}