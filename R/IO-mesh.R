#' Read binary PLY file
#' 
#' Read vertices, faces and color or quality of the vertices of a mesh from a PLY file in little endian 1.0 format
#'
#' @param file_name path to the PLY file
#' 
#' @return a list of three components of the mesh, the vertices, the faces and the color or quality of the vertices
#' @export
read_binary_PLY <- function(filePath) {
    # Call to C++ code that reads the binary file and write it to a temp file
    path <- read_PLY(filePath, file.path(tempdir(),"temp.dat") )
    data <- scan(file.path(tempdir(),"temp.dat"))
    
    # Read vertices, faces and quality
    vertex <- matrix(data[1:(4 * path[[1]])], ncol = 4, byrow = TRUE)
    faces <- matrix(data[((4 * path[[1]]) + 1):length(data)], ncol = 4, byrow = TRUE)
    
    faces <- faces[, 2:4]
    quality <- vertex[ ,4]
    vertices <- vertex[ ,1:3]
    
    # Add 1 to the faces because index in R starts at 1
    faces <- faces + 1
    
    # Remove temp file
    file.remove( file.path(tempdir(),"temp.dat")  )
    return(list(vertices = vertices, faces = faces, quality = quality))
}

#' Write only the vertices of a mesh in a PLY file
#' 
#' Write only the vertices of a mesh in a PLY file so it is saved a point cloud
#'
#' @param vertices matrix Nx3 where each row represents the coordinates of a vertex of the mesh
#' @param file_name path and name of the PLY file where point cloud is going to be saved
#' 
#' @return None
write_PLY_only_points <- function(vertices, file_name = "")
{
  if(file_name==""){file_name <- tclvalue(tkgetSaveFile())}
  
  file_name <- sub("\\.ply$","", file_name, ignore.case = T)
  
  # Write as matrix
  return ( vcgPlyWrite(vertices, filename = file_name, binary = T) )
  
}

#' Read neuron from a VRML file
#' 
#' Read vertices, faces of a mesh from a VRML file. Adapted from package geomorph
#'
#' @param file_name path to the PLY file
#' 
#' @return a list of two components of the mesh, the vertices and the faces
#' @export
read_VRML <- function(file_name) {
    # Read each caracter one by one.
    b <- scan(file_name, what = "character", skip = 0, sep = "", quiet = T)
    
    cos <- which(b == "Coordinate")  # 'Coordinate' denotes that the definition of vertices starts.
    tri <- which(b == "coordIndex")  # 'coordIndex' denotes that the definition of faces starts.
    zz <- triang <- NULL  # Data frames of vertices and faces.
    co <- c()  # Positions of 'b' that identifies the beginning of the vertices
    en <- c()  # Positions of 'b' that identifies the end of the vertices
    
    # Loop to find the start and the beginning of each block of vertices.
    for (i in 1:length(cos)) {
        co[i] <- cos[i] + 4  # Beginning of the block is 4 character after 'Coordinate'
        en[i] <- which(b[co[i]:length(b)] == "]")[1] - 1  #']' denotes that block is finished so that it should be removed the last character
    }
    
    # Select that with more vertices
    positionCoords <- which(en == sort(en, T)[1])
    
    
    # Remove incorrect coordinates
    zerocoord <- NULL
    # Select vertices that only belongs to the neuron.
    newrange <- b[co[positionCoords]:length(b)]
    # Check that there are vertices
    if (en[positionCoords] != 0) {
        # Selected values are disposed in a matrix
        z <- as.matrix(newrange[1:en[positionCoords]])
        # There must be 3 columns if not there is a definition error
        if (nrow(z)%%3 == 0) {
            # Matrix is disposed in 3 columns
            zz.temp <- matrix(as.numeric(unlist(strsplit(z, split = ","))), ncol = 3, byrow = T)
            zerocoord <- which(as.character(zz.temp[, 1]) == as.character(0) & as.character(zz.temp[, 2]) == 
                as.character(0) & as.character(zz.temp[, 3]) == as.character(0))
            # Wrong positions are removed
            if (length(zerocoord) > 0) {
                zz.temp <- zz.temp[-zerocoord, ]
            }
            
            zz <- rbind(zz, zz.temp)
        }
    }
    
    # Similar to the vertices definition applied to the faces
    co.tri <- c()
    en.tri <- c()
    for (i in 1:length(tri)) {
        co.tri[i] <- tri[i] + 2
        en.tri[i] <- which(b[co.tri[i]:length(b)] == "]")[1] - 1
    }
    position <- which(en.tri == sort(en.tri, T)[1])
    newrange.tri <- b[co.tri[position]:length(b)]
    if (en.tri[position] != 0) {
        tria <- as.matrix(newrange.tri[1:en.tri[position]])
        
        if (nrow(tria)%%4 == 0) {
            triang.temp <- matrix(as.numeric(unlist(strsplit(tria, split = ","))), ncol = 4, byrow = T)[, 
                -4]
            if (i == 1) {
                triang <- triang.temp + 1
            } else {
                triang <- rbind(triang, triang.temp + 1)
            }
        }
    }
    colnames(zz) <- c("x", "y", "z")
    
    return(list(vertices = zz, faces = triang))
} 

