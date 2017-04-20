#' Compute RMSE 
#' 
#' Compute RMSE between three different techniques to compare processing methods
#'
#' @param technique_1 path to the folder where meshes processed with the first technique are placed. They have to be PLY files
#' @param technique_2 path to the folder where meshes processed with the second technique are placed. They have to be PLY files
#' @param technique_3 path to the folder where meshes processed with the third technique are placed. They have to be PLY files
#' 
#' @return RMSE a matrix Nx3 where N is the number of meshes and first column is RMSE between technique 1 and technique 2, second column is RMSE between technique 1 and technique 3 and third column is RMSE between technique 2 and technique 3
#' 
#' @examples 
#' ######################################################
#' ###################### Interexpert ###################
#' ######################################################
#' #RMSE between the somas of the algorithm and the somas of the experts before repairing the experts' somas
#' path_somas_algorithm<- "./temp/final_result"
#' path_somas_experts<- system.file("test/pre_repaired",package="SomaMS")
#' experts_paths<-list.dirs(path_somas_experts,recursive=F)
#' pre_repaired_RMSE<-RMSE_mesh_distance(path_somas_algorithm, experts_paths[1], experts_paths[2], TRUE)
#' 
#' #RMSE between the somas of the algorithm and the somas of the experts after repairing the experts' somas
#' path_somas_algorithm<- "./temp/final_result"
#' path_somas_experts<- system.file("test/post_repaired",package="SomaMS") 
#' experts_paths<-list.dirs(path_somas_experts,recursive=F)
#' post_repaired_RMSE<-RMSE_mesh_distance(path_somas_algorithm, experts_paths[1], experts_paths[2], TRUE)
#' 
#' X11(width = 18, height = 10.37)
#' par(mfrow=c(1,2))
#' values_barplot <- t(pre_repaired_RMSE)
#' colors <- c(rainbow(3))
#' mp <- barplot2(values_barplot, main = "RMSE before repairing experts' somas", ylab = "RMSE", beside = TRUE,
#'                col = colors, ylim = c(0, 1.7), cex.names = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#' legend("top", legend = c("Procedure  Vs  Expert 1", "Procedure  Vs  Expert 2", "Expert 1     Vs  Expert 2"), fill = colors, box.col = "transparent", x.intersp = 0.8, cex = 1.5)

#' mtext(1, at = mp[2,], text = c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Neuron 5", "Neuron 6", "Neuron 7", "Neuron 8", "Neuron 9"), line = 0.5, cex = 1)
#' legend("topleft", legend = "", title = "A", box.col = "transparent", cex = 2)

#' values_barplot <- t(post_repaired_RMSE)
#' colors <- c(rainbow(3))
#' mp <- barplot2(values_barplot, main = "RMSE after repairing experts' somas", ylab = "RMSE", beside = TRUE,
#'                col = colors, ylim = c(0, 1.7), cex.names = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#' legend("top", legend = c("Procedure  Vs  Expert 1", "Procedure  Vs  Expert 2", "Expert 1     Vs  Expert 2"), fill = colors, box.col = "transparent", x.intersp = 0.8, cex = 1.5)

#' mtext(1, at = mp[2,], text = c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Neuron 5", "Neuron 6", "Neuron 7", "Neuron 8", "Neuron 9"), line = 0.5, cex = 1)
#' legend("topleft", legend = "", title = "B", box.col = "transparent", cex = 2)
#' 
#' print(paste0("Mean interexpert RMSE: ", mean(post_repaired_RMSE[,3])))
#' ######################################################
#' ###################### Intraexpert ###################
#' ######################################################
#' path_somas_experts<- system.file("test/intraexpert",package="SomaMS")
#' experts_paths<-list.dirs(path_somas_experts,recursive=F)
#' 
#' expert_1_days<-list.dirs(experts_paths[1],recursive=F)
#' expert_2_days<-list.dirs(experts_paths[2],recursive=F)
#' 
#' intraexpert_expert_1 <- RMSE_mesh_distance(expert_1_days[1], expert_1_days[2], expert_1_days[3], TRUE)
#' intraexpert_expert_2 <- RMSE_mesh_distance(expert_2_days[1], expert_2_days[2], expert_2_days[3], TRUE)
#' 
#' X11(width = 18, height = 10.37)
#' par(mfrow = c(1,2))
#' 
#' valuesBarplot <- t(intraexpert_expert_1)
#' colors <- c(rainbow(3))
#' mp <- barplot2(valuesBarplot, main="Intra-expert variability of the first expert ", ylab = "RMSE", beside = TRUE,
#'                col = colors, ylim = c(0,1.6), cex.names = 1.5,cex.lab = 1.5, cex.axis = 1.5)
#' legend("top",legend = c("Day 1 Vs Day 2","Day 1 Vs Day 3","Day 2 Vs Day 3"), fill = colors, box.col = "transparent", x.intersp = 0.8, cex = 1.5)
#' 
#' mtext(1, at = mp[2,], text = c("Neuron 1","Neuron 2","Neuron 3","Neuron 4","Neuron 5","Neuron 6"),line = 0.5, cex = 1.3)
#' legend("topleft",legend="",title="A",box.col = "transparent",cex=2)
#' 
#' 
#' values_barplot<-t(intraexpert_expert_2)
#' colors<-c(rainbow(3))
#' mp <- barplot2(values_barplot, main = "Intra-expert variability of the second expert ", ylab = "RMSE", beside = TRUE,
#'                col = colors, ylim = c(0, 1.6), cex.names = 1.5, cex.lab = 1.5, cex.axis = 1.5)
#' legend("top",legend = c("Day 1 Vs Day 2","Day 1 Vs Day 3","Day 2 Vs Day 3"), fill = colors, box.col = "transparent", x.intersp = 0.8, cex = 1.5)
#' 
#' mtext(1, at = mp[2,], text = c("Neuron 1", "Neuron 2", "Neuron 3", "Neuron 4", "Neuron 5", "Neuron 6"), line = 0.5, cex = 1.3)
#' legend("topleft", legend = "", title = "B", box.col = "transparent", cex = 2)
#' 
#' print(paste0("Mean RMSE for the first expert: ", mean(c(intraexpert_expert_1))))
#' print(paste0("Mean RMSE for the second expert: ", mean(c(intraexpert_expert_2))))
#' 
#' @export
RMSE_computation <- function(path_somas_algorithm, path_somas_expert_1, path_somas_expert_2, parallel = TRUE)
{
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

  RMSE <- RMSE_between_meshes(path_somas_algorithm, path_somas_expert_1, cl)
  RMSE <- cbind(RMSE, RMSE_between_meshes(path_somas_algorithm, path_somas_expert_2, cl))
  RMSE <- cbind(RMSE, RMSE_between_meshes(path_somas_expert_1, path_somas_expert_2, cl))
  
  if(ncores > 1)
  {
    stopCluster(cl)
  }
  
  return(RMSE)
}


RMSE_between_meshes <- function(path_somas_algorithm, path_somas_expert, cl = NULL)
{
  list_interexpert_somas <- list.files(path_somas_expert)

  if(is.null(cl))
  {
    distance_matrix<-foreach(i = 1:length(list_interexpert_somas), .combine=c, .packages=c('Rvcg','Morpho')) %do%
    {
      soma_algorithm <- vcgPlyRead(file.path(path_somas_algorithm, list_interexpert_somas[i]))
      soma_expert <- vcgPlyRead(file.path(path_somas_expert, list_interexpert_somas[i]))
  
      max(sqrt(mean((meshDist(soma_algorithm, soma_expert, plot = FALSE)$dists)^2)), sqrt(mean((meshDist(soma_expert, soma_algorithm, plot = FALSE)$dists)^2)))
    }
    
  }else{
    distance_matrix<-foreach(i = 1:length(list_interexpert_somas), .combine=c, .packages=c('Rvcg','Morpho')) %dopar%
    {
      soma_algorithm <- vcgPlyRead(file.path(path_somas_algorithm, list_interexpert_somas[i]))
      soma_expert <- vcgPlyRead(file.path(path_somas_expert, list_interexpert_somas[i]))
   
      max(sqrt(mean((meshDist(soma_algorithm, soma_expert, plot = FALSE)$dists)^2)), sqrt(mean((meshDist(soma_expert, soma_algorithm, plot = FALSE)$dists)^2)))
    }
  }
  return (distance_matrix)
}


#' Wilcoxon test for RMSE
#' 
#' Compute Wilcoxon between RMSE distances to find significant differences between mesh processing techniques
#'
#' @param RMSE_matrix a matrix Nx3 which is the output of 
#' 
#' @return A vector of two p-values. The first value is the p-value resulting from wilcoxon between first column of the RMSE matrix and second column of the RMSE matrix. The second one is the p-value resulting from wilcoxon between first column of the RMSE matrix and second column of the RMSE matrix 
#' 
#' @examples 
#' path_somas_algorithm<- "./test/final_result"
#' path_somas_experts<- system.file("test/pre_repaired",package="SomaMS")
#' experts_paths<-list.dirs(path_somas_experts,recursive=F)
#' pre_repaired_RMSE<-RMSE_mesh_distance(path_somas_algorithm, experts_paths[1], experts_paths[2], TRUE)
#' 
#' #RMSE between the somas of the algorithm and the somas of the experts after repairing the experts' somas
#' path_somas_algorithm<- "./test/final_result"
#' path_somas_experts<- system.file("test/post_repaired",package="SomaMS")
#' experts_paths<-list.dirs(path_somas_experts,recursive=F)
#' post_repaired_RMSE<-RMSE_mesh_distance(path_somas_algorithm, experts_paths[1], experts_paths[2], TRUE)
#' 
#' pvalue_before_repairing <- wilcoxon_RMSE(pre_repaired_RMSE)
#' print(paste0("p-value between algorithm and first expert: ", pvalue_before_repairing[1]))
#' print(paste0("p-value between algorithm and second expert: ", pvalue_before_repairing[2]))
#' 
#' pvalue_after_repairing <- wilcoxon_RMSE(post_repaired_RMSE)
#' print(paste0("p-value between algorithm and first expert: ", pvalue_after_repairing[1]))
#' print(paste0("p-value between algorithm and second expert: ", pvalue_after_repairing[2]))
#'
#' @export
wilcoxon_RMSE<-function(RMSE_matrix)
{
  return (c(wilcox.test(RMSE_matrix[ ,1], RMSE_matrix[ ,3], paired = TRUE, exact = TRUE)$p.value,wilcox.test(RMSE_matrix[ ,2], RMSE_matrix[ ,3], paired = TRUE, exact = TRUE)$p.value))
}


#' Compute the volume of meshes in a folder
#' 
#' Read one by one the meshes in a folder and compute the volume of each one of the them 
#'
#' @param path_to_somas_folder path to the folder where meshes processed with the first technique are placed. They have to be PLY files
#' 
#' @return array of N volumes where N is the number of meshes in the folder
#' 
#' @examples
#' compute_meshes_volumes("./temp/final_result")
#' 
#' @export
compute_meshes_volumes<-function(path_to_mesh_folder)
{
  list_interexpert_somas <- list.files(path_to_mesh_folder)
  volumes <- foreach(i=1:length(list_interexpert_somas), .combine=c, .packages=c('Rvcg','Morpho')) %do%
  {
    soma_mesh <- vcgPlyRead(file.path(path_to_mesh_folder, list_interexpert_somas[i], sep = ""))
    poly_volume(t(soma_mesh$vb[1:3, ])[soma_mesh$it, ])
  }
  
  return (volumes)
}


#' Compute MAQ_S
#' 
#' Compute MAQ_S between three different techniques to compare processing methods
#'
#' @param path_somas_algorithm path to the folder where meshes processed with the first technique are placed. They have to be PLY files
#' @param path_somas_expert_1 path to the folder where meshes processed with the second technique are placed. They have to be PLY files
#' @param path_somas_expert_2 path to the folder where meshes processed with the third technique are placed. They have to be PLY files
#' 
#' @return A vector of three values where each value is the MAQ_S between two techniques. Concretely they are MAQ_S_12, MAQ_S_13, MAQ_S_23
#' 
#' @examples 
#' path_somas_algorithm<- "./temp/final_result"
#' path_somas_experts<- system.file("test/post_repaired",package="SomaMS")
#' experts_paths<-list.dirs(path_somas_experts,recursive=F)
#' path_somas_expert_1<-experts_paths[1]
#' path_somas_expert_2<-experts_paths[2]

#' MAQ_S_result <- MAQ_S(path_somas_algorithm,path_somas_expert_1, path_somas_expert_2)

#' print(paste("MAQ_S_12 is:",   MAQ_S_result[1] * 100, "%"))
#' print(paste("MAQ_S_13 is:",   MAQ_S_result[2] * 100, "%"))
#' print(paste("MAQ_S_23 is:",   MAQ_S_result[3] * 100, "%"))
#' print(paste("Mean MAQ_S between algorithm and experts is:", mean(c(MAQ_S_result[1],MAQ_S_result[2])) * 100, "%"))
#' print(paste("Difference between experts MAQ_S and mean MAQ_S of algorithm is:", abs(MAQ_S_result[3]-mean(c(MAQ_S_result[1],MAQ_S_result[2])))*100, "%"))
#' 
#' @export
MAQ_S<-function(path_somas_algorithm, path_somas_expert_1, path_somas_expert_2)
{
  list_interexpert_somas <- list.files(path_somas_expert_1)
  
  volumes <- matrix(nrow=length(list_interexpert_somas),ncol=3)
  volumes[ ,1] <- compute_meshes_volumes(path_somas_algorithm)
  volumes[ ,2] <- compute_meshes_volumes(path_somas_expert_1)
  volumes[ ,3] <- compute_meshes_volumes(path_somas_expert_2)
  
  #MEAN ABSOLUTE QUOTIENT VOLUME
  MAQ_S_12 <- max(mean(abs((volumes[ ,1] / volumes[ ,2]) - 1)), mean(abs((volumes[ ,2] / volumes[ ,1]) - 1))) 
  MAQ_S_13 <- max(mean(abs((volumes[ ,1] / volumes[ ,3]) - 1)), mean(abs((volumes[ ,3] / volumes[ ,1]) - 1))) 
  MAQ_S_23 <- max(mean(abs((volumes[ ,2] / volumes[ ,3]) - 1)), mean(abs((volumes[ ,3] / volumes[ ,2]) - 1)))

  return(c(MAQ_S_12, MAQ_S_13, MAQ_S_23))
}

