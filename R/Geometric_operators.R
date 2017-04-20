
#Se recibe un conjunto de puntos que conforman aproximadamente una elipse. Los puntos son transladados al plano XY donde se definen eje mayor y menor.
#A continuacion se pasa 1 de los puntos del eje mayor y 1 punto del eje menor a coordenadas sphericas. Se pretende representar una elipse en base a sus diagonales.
ellipseSphPoints<-function(polygon)
{
	#Se calcula el centroide del poligono
	elipse_center_of_mass <- colSums(polygon) / nrow(polygon)
	
	#Calculo del PCA para obtener la rotacion al plano XY. Rotacion al plano XY
	component_analysis <- prcomp(polygon)
	rotated_ellipse <- rotate3d(polygon, matrix = component_analysis$rotation)

	#Se cogen 2 puntos. Uno para cada elipse
	point_over_diameter_X <- polygon[which.max(rotated_ellipse[ ,1]), ]
	point_over_diameter_Y <- polygon[which.max(rotated_ellipse[ ,2]), ]

	#Se calcula el vector entre el centroide del poligono y los puntos sobre los ejes que se calcularon. 
	#De esta forma el centroide del poligono pasa a ser el centro del sistema de coordenadas esfericas.
	v_point_X <- point_over_diameter_X - elipse_center_of_mass
	v_point_Y <- point_over_diameter_Y - elipse_center_of_mass
	
	#Finalmente se almacenan las coordenadas esfericas para ambos puntos respecto al centro de masas de la elipse
	spherical_X <- cartesian2Spherical(c(v_point_X))
	spherical_Y <- cartesian2Spherical(c(v_point_Y))
	return(c(spherical_X, spherical_Y))
}


ellipseCartPoints<-function(ellipse_diameters, ellipse_center_of_mass, number_of_points = 1000)
{
	#Se almacenan las coordenadas spherical con la mitad del diametro de la elipse
	spherical_coordinates_X <- c(ellipse_diameters[1], ellipse_diameters[2], ellipse_diameters[3])
	spherical_coordinates_Y <- c(ellipse_diameters[4], ellipse_diameters[5], ellipse_diameters[6])
	
	#Se pasan las coordenadas de esfericas a cartesianas     
	cartesian_coordinatesX <- spherical2Cartesian(c(spherical_coordinates_X)) 	
	cartesian_coordinatesY <- spherical2Cartesian(c(spherical_coordinates_Y)) 
	
	#Se generan 2 puntos alrededor del centro de masas de la elipse que se situaran sobre el eje X cuando se rote la elipse al plano XY.
	point_X1 <- ellipse_center_of_mass + cartesian_coordinatesX 
	point_X2 <- ellipse_center_of_mass - cartesian_coordinatesX 

	#Se generan 2 puntos alrededor del centro de masas de la elipse que se situaran sobre el eje Y cuando se rote la elipse al plano XY.
	point_Y1 <- ellipse_center_of_mass + cartesian_coordinatesY
	point_Y2 <- ellipse_center_of_mass - cartesian_coordinatesY

	#Se almacenan los 5 puntos que se van a utilizar en el calculo de la elipse
	ellipse4points <- rbind(pointX1, pointX2, pointY1, pointY2, ellipse_center_of_mass)

	#Se colocan los 5 puntos sobre el plano XY para generar la elipse
	PCA <- prcomp(ellipse4points)
	rotated4_points_ellipse <- rotate3d(ellipse4points, matrix = PCA$rotation)

	#Se genera una elipse con numberOfPoints puntos o en su defecto 1000 sobre XY
	random_points <- runif(number_of_points, min = 0, max = 2 * pi)
	a <- sqrt(sum((rotated4_points_ellipse[5, ] - rotated4_points_ellipse[1, ])^2))
	b <- sqrt(sum((rotated4_points_ellipse[5, ] - rotated4_points_ellipse[3, ])^2))
	x <- a * cos(random_points) + rotated4_points_ellipse[5,1]
	y <- b * sin(random_points) + rotated4_points_ellipse[5,2]

	#Una vez se dispone de la elipse se vuelve a colocar la elipse
	reconstructed_ellipse <- cbind(x, y, rotated4_points_ellipse[1,3])
	final_ellipse<-rotate3d(reconstructed_ellipse, matrix = t(PCA$rotation))
	return(final_ellipse)
}

#Funcion que recibe n puntos en sus coordenadas x e y a partir de los cuales calcula la elipse que mejor se ajusta a esos puntos.
#Input:
#	-x:Valores del punto 2 dimensional para el eje x.
#	-y:Valores del punto 2 dimensional para el eje y.
#Output:
#	-Lista de parametros: Se almacenan los parametros necesarios para representar una elipse como el centroide, la longitud de los ejes, el angulo que forman los ejes y los coeficientes de la regresion.
fit.ellipse <- function (x, y = NULL) {
  
  EPS <- 1.0e-8 
  dat <- xy.coords(x, y) 
  
  D1 <- cbind(dat$x * dat$x, dat$x * dat$y, dat$y * dat$y) 
  D2 <- cbind(dat$x, dat$y, 1) 
  S1 <- t(D1) %*% D1 
  S2 <- t(D1) %*% D2 
  S3 <- t(D2) %*% D2 
  T <- -solve(S3) %*% t(S2) 
  M <- S1 + S2 %*% T 
  M <- rbind(M[3,] / 2, -M[2,], M[1,] / 2) 
  evec <- eigen(M)$vec 
  cond <- 4 * evec[1,] * evec[3,] - evec[2,]^2 
  a1 <- evec[, which(cond > 0)] 
  f <- c(a1, T %*% a1) 
  names(f) <- letters[1:6] 
  

  A <- matrix(c(2*f[1], f[2], f[2], 2*f[3]), nrow=2, ncol=2, byrow=T )
  b <- matrix(c(-f[4], -f[5]), nrow=2, ncol=1, byrow=T)
  soln <- solve(A) %*% b

  b2 <- f[2]^2 / 4
  
  center <- c(soln[1], soln[2]) 
  names(center) <- c("x", "y") 
  
  num  <- 2 * (f[1] * f[5]^2 / 4 + f[3] * f[4]^2 / 4 + f[6] * b2 - f[2]*f[4]*f[5]/4 - f[1]*f[3]*f[6]) 
  den1 <- (b2 - f[1]*f[3]) 
  den2 <- sqrt((f[1] - f[3])^2 + 4*b2) 
  den3 <- f[1] + f[3] 
  
  semi.axes <- sqrt(c( num / (den1 * (den2 - den3)),  num / (den1 * (-den2 - den3)) )) 
  
  # calculate the angle of rotation 
  term <- (f[1] - f[3]) / f[2] 
  angle <- atan(1 / term) / 2 
  
  list(coef = f, center = center, major = max(semi.axes), minor = min(semi.axes), angle = unname(angle)) 
}

#Generador de puntos de una elipse en base a los parametros que se ajustaron en fit.ellipse
#Input:
#	-fit:Lista de valores que se obtuvo en fit.ellipse
#	-n:Numero de puntos que se quieren generar para la elipse.
#Output:
#	-x:Valores sobre el eje x para los puntos que se han generado para la elipse.
#	-y:Valores sobre el eje y para los puntos que se han generado para la elipse.
get.ellipse <- function(fit, n = 360 ) 
{
  
  tt <- seq(0, 2*pi, length=n) 
  sa <- sin(fit$angle) 
  ca <- cos(fit$angle) 
  ct <- cos(tt) 
  st <- sin(tt) 
  
  x <- fit$center[1] + fit$maj * ct * ca - fit$min * st * sa 
  y <- fit$center[2] + fit$maj * ct * sa + fit$min * st * ca 
  
  cbind(x = x, y = y) 
}

 
trans3d <- function(vector)
{
        trans_matrix <- rbind(c(1, 0, 0, 0), c(0, 1, 0, 0), c(0, 0, 1, 0), c(vector[1], vector[2], vector[3], 1))
		return(trans_matrix)
}

rotation3dX<-function(radians)
{
        c <- cos(radians)
        s <- sin(radians)
        return( rbind(c(1,0,0,0),c( 0,c,s,0),c( 0,-s,c,0),c( 0,0,0,1)))
}

rotation3dY<-function(radians)
{
        c <- cos(radians)
        s <- sin(radians)
        return( rbind(c(c,0,-s,0),c( 0,1,0,0),c( s,0,c,0),c( 0,0,0,1)))
}

#u es un vector perpendicular a ambos vectores que se desean alinear y theta el angulo entre los vectores
rotation_matrix<-function(u, theta)
{ 
	rotationMatrix <- t(matrix(c(cos(theta)+u[1]^2*(1-cos(theta)), u[1]*u[2]*(1-cos(theta))-u[3]*sin(theta), u[1]*u[3]*(1-cos(theta))+u[2]*sin(theta),
    u[2]*u[1]*(1-cos(theta))+u[3]*sin(theta), cos(theta)+u[2]^2*(1-cos(theta)), u[2]*u[3]*(1-cos(theta))-u[1]*sin(theta),
    u[3]*u[1]*(1-cos(theta))-u[2]*sin(theta), u[3]*u[2]*(1-cos(theta))+u[1]*sin(theta), cos(theta)+u[3]^2*(1-cos(theta))),nrow=3,ncol=3))
	
	return(rotationMatrix)
}

align3dZ <-function(vector)
 {   
	#Translacion al origen
    trans<-trans3d(-vector)
	itrans <- trans3d(vector);
	
    u <- vector/sqrt(sum(vector^2));
    uz <- u;  uz[1] = 0; uz[2] = 0; 
    uzy <- u; uzy[1] = 0;
    	

    modUz <- uz[3];
    modUzy <- sqrt(sum(uzy^2))
    radsX <- acos(modUz/modUzy);
	if (u[2] < 0) {
        radsX <- -radsX;
    }
    rotX <- rotation3dX(radsX);
    irotX <- rotation3dX(-radsX);
    
    modUzx <- 1; 
    modUzBis <- modUzy; 
    radsY <- acos(modUzBis / modUzx);
    if (u[1] < 0) {
        radsY = -radsY;
    }
	
    rotY <- rotation3dY(-radsY);
    irotY <- rotation3dY(radsY);
	return(list(trans, rotX, rotY))
}

#Convertir coordenadas 3D de del sistema cartesiana al sistema esferico. Se basa en la definicion que se prove en 
#http://es.wikipedia.org/wiki/Coordenadas_esf%C3%A9ricas#Relaci.C3.B3n_con_las_coordenadas_cartesianas 
#Input:
#	-Matriz de coordenadas 3D en el que una fila es un punto tridimensional y las columnas son las dimensiones x,y,z.
#Output:
#	-Matriz de coordenadas esfericas en el que cada fila es un punto en coordenadas esfericas y las columnas hacen referencia a la colatitud, azimuth y el radio en ese orden.
# cartesian2Spherical <- function(cartesian_coord_matrix)
# {
# 	 #El caso en el
# 	 if (cartesian_coord_matrix[3] < 1.0e-8 & cartesian_coord_matrix[3] > -1.0e-8)
#      {
# 		colatitude <- pi / 2;
#      }else if (cartesian_coord_matrix[3] > 0)
#            {
# 				colatitude <- atan(sqrt(cartesian_coord_matrix[1]^2 + cartesian_coord_matrix[2]^2) / cartesian_coord_matrix[3]);
#            }else{
#                 colatitude <- atan(sqrt(cartesian_coord_matrix[1]^2 + cartesian_coord_matrix[2]^2) / cartesian_coord_matrix[3]) + pi;
# 	}
#             
#     if (cartesian_coord_matrix[1] < 1.0e-8 & cartesian_coord_matrix[1] > -1.0e-8){
#                     azimuth <- (pi/2) * sign(cartesian_coord_matrix[2])
#     }else if(cartesian_coord_matrix[1] > 0 & cartesian_coord_matrix[2] > 0){
#                     azimuth <- atan(cartesian_coord_matrix[2] / cartesian_coord_matrix[1])
#     }else if(cartesian_coord_matrix[1] > 0 & cartesian_coord_matrix[2] < 0){
#                     azimuth <- 2 * pi + atan(cartesian_coord_matrix[2] / cartesian_coord_matrix[1])
#     }else{ 
#                     azimuth <- pi + atan(cartesian_coord_matrix[2] / cartesian_coord_matrix[1])
#     }      
# 	  
# 	r <- sqrt(sum(cartesian_coord_matrix^2))
#     return(c(colatitude, azimuth, r))  	
# }

cartesian2Spherical <- function(x,y,z)
{
  azimuth <- atan2(y,x)
  elevation <- atan2(z,sqrt(x^2 + y^2))
  r <- sqrt(x^2 + y^2 + z^2)
  return(cbind(azimuth,elevation,r))
}

#Convertir coordenadas 3D de del sistema esferico al cartesiano. Se basa en la definicion que se prove en 
#http://es.wikipedia.org/wiki/Coordenadas_esf%C3%A9ricas#Relaci.C3.B3n_con_las_coordenadas_cartesianas 
#Input:
#	-Matriz de coordenadas 3D en el que una fila es un punto tridimensional y las columnas son las dimensiones x,y,z.
#Output:
#	-Matriz de coordenadas esfericas en el que cada fila es un punto en coordenadas esfericas y las columnas hacen referencia a la colatitud, azimuth y el radio en ese orden.
spherical2Cartesian<-function(azimuth,elevation,r)
{  
  x <- r * cos(elevation) * cos(azimuth)
  y <- r * cos(elevation) * sin(azimuth)
  z <- r * sin(elevation)
	return(cbind(x, y, z))
}

cart2pol <- function(x, y)
{
  r <- sqrt(x^2 + y^2)
  theta <- atan(y/x)
  
  return(c(theta,r))
}

pol2cart <- function(theta,r)
{
  x<-r*cos(theta)
  y<-r*sin(theta)
  return(c(x,y))
}

vcrossp <- function( a, b ) { 
  result <- matrix( NA, 1, 3 ) 
  result[1] <- a[2] * b[3] - a[3] * b[2] 
  result[2] <- a[3] * b[1] - a[1] * b[3] 
  result[3] <- a[1] * b[2] - a[2] * b[1] 
  return(result) 
} 

normv<- function(a)
{
  a<-as.numeric(a)
  return(sqrt(a%*%a))
}

