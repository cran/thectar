###################################################################################
### This file contains functions useful during the modelling process
###################################################################################

###################################################################################
### lspoints (lowering Stress by excluding points)
###################################################################################

#' Lowering Stress by excluding points (lspoints)
#' 
#' \code{lspoints} calculates Stress values for all combinations with n - 1 
#' categories. 
#' 
#'  This function is applicable to co-occurrence data. It shows the resulting 
#'  Stress values when single categories are excluded. The output consists in a 
#'  table showing which categories have been excluded and the resulting Stress 
#'  values. The table is sorted such that the the lowest Stress level occurs 
#'  at the top. The MDS model is computed using the package '\code{smacof}' 
#'  (Mair, De Leeuw, Borg, & Groenen). The first column of the input matrix 
#'  \code{data} should contain the ID of the unit of comparison, and the 
#'  following columns the categories for which the similarity matrices and MDS 
#'  maps are calculated (see also the description of \code{simi}). Note that 
#'  the purpose of this function is to assist modeling by helping to identify 
#'  potential problems. It is not, however, meant to be used for excluding 
#'  categories based solely on measures of fit and without substantial 
#'  justification.
#'  
#' @param data Dataset; the first column must be the ID of the unit of 
#' comparison, the other columns must be categories; see dataset requirements 
#' \code{simi}. 
#' @param method Specifies the similarity index used to compute the similarity 
#'  matrix, choose between "\code{as}" (Association Strength Index), 
#'  "\code{jaccard}" (Jaccard  Index), "\code{cosine}" (Cosine Index), 
#'  and "\code{inclusion}" (Inclusion Index). Default is \code{as}. 
#' @param single If \code{TRUE}, single mentionings (i.e. one respondent 
#' mentioning just one category) are included when calculating the similarity 
#' matrix.
#' @param comments If \code{TRUE}, comments relating to exclusion or possible 
#'  exclusion of categories and respondents are displayed. Default is 
#'  \code{FALSE}. 
#' @param type Specifies the type of MDS model used (for more details see
#' the package '\code{smacof}' by Mair, De Leeuw, Borg, & Groenen). Default is
#' \code{ordinal}.
#' @param weight If \code{TRUE}, the MDS model is calcualted using weights, 
#' i.e., similarities of zero are excluded. Default is \code{TRUE}. 
#' @return Matrix showing the categories excluded and the Stress values of the 
#' respective MDS models.
#' @seealso \code{\link[smacof]{smacofSym}} for details on calculating MDS 
#' representations, \code{\link[thectar]{simi}} for details on calculating 
#' similarity matrices.
#' @export
#' @examples
#' ## Calculate Stress values for all combinations with n - 1 (i.e., 16) SDGs
#' data(SDG_coocurrence)
#' SDG_coocurrence <- SDG_coocurrence[,-2] # Drop second column
#' stress <- lspoints(SDG_coocurrence)
#' stress


lspoints <- function(data, method = "as", single = TRUE, comments = FALSE, type = "ordinal", weight = TRUE){
  
  if(any(is.na(data))){

    stop("The dataset contains missing values.")
    
  }else{
    
    if(!(method %in% c("as",  "jaccard","cosine",  "inclusion"))){
      
      method <- "as"
      warning("The parameter 'method' is not properly defined. Default will be used.")
    } #if
    
    # Prepare basic variables 
    
    i <- 2
    output <- matrix(c(rep(0, (length(data) - 1) * 2)), ncol = 2)
    
    # Exclude the points and save the Stress values
    
    while(i < length(data) + 1){
      configuration <- data
      configuration <- configuration[, -i]
      
      similarity <- simi(configuration, method=method, single = single, comments = comments)
      
      if(weight == TRUE){
        w <- as.numeric(similarity < 1e-8)
        weight_mat = matrix(1 - w, nrow = nrow(similarity), ncol = ncol(similarity))
        dissimilarities <- 1 - similarity
        res<-smacof::smacofSym(dissimilarities,  type = type, weightmat = weight_mat)
      }else{
        dissimilarities <- 1 - similarity
        res<-smacof::smacofSym(dissimilarities,  type = type)
      }#ifelse
      output[i -1, 1] <- (colnames(data[i]))
      output[i - 1, 2] <- res$stress
      i <- i + 1
    } #while
    
    # Order the resulting vector such that the lowest Stress value is at the top 
    
    r <- as.numeric(output[, 2])
    output <- output[order(r), ]
    
    # Add labels to the columns 
    
    labels <- c("Excluded object", "Resulting Stress value")
    colnames(output) <- labels
      
  }#ifelse


  # Output 
  
  return(output)

}



###################################################################################
### lscomb (lowering Stress by finding the 'ideal' combination of points)
###################################################################################

#' Lowering Stress by comparing combinations (lscomb)
#' 
#' \code{lscomb} calculates Stress levels for all combinations of p out of 
#' n categories 
#' 
#' This function is applicable to co-occurrence data. It shows the resulting 
#' Stress values for all combinations of p out of n categories. The output 
#' consists in a table showing which categories have been included and the 
#' resulting Stress values. The table is sorted such that the lowest Stress 
#' level occurs at the top. The MDS model is computed using the package 
#' '\code{smacof}' (Mair, De Leeuw, Borg, & Groenen). The first column of the 
#' input matrix \code{data} should contain the ID of the unit of comparison 
#' and the following columns the categories for which the similarity matrices 
#' and MDS maps are calculated (see also the description of \code{simi}). Note 
#' that the purpose of this function is to assist modeling by helping to 
#' identify potential problems. It is not, however, meant to be used for 
#' excluding categories based solely on measures of fit and without substantial 
#' justification.
#' 
#' @param data Dataset; the first column must be the ID of the unit of comparison, 
#'  the other columns must be categories; see dataset requirements \code{simi}. 
#' @param points Number of categories that should be included in the model 
#' (p < n, wherein n equals the number of categories in the dataset). 
#' @param method Specifies the similarity index used to compute the similarity 
#'  matrix, choose between "\code{as}" (Association Strength Index), 
#'  "\code{jaccard}" (Jaccard  Index), "\code{cosine}" (Cosine Index), 
#'  and "\code{inclusion}" (Inclusion Index). Default is \code{as}.  
#' @param single If \code{TRUE}, single mentionings (i.e. one respondent 
#' mentioning just one category) are included when calculating the similarity 
#' matrix. Default is \code{TRUE}. 
#' @param comments If \code{TRUE}, comments relating to exclusion or possible 
#'  exclusion of categories and respondents are displayed. Default is 
#'  \code{FALSE}. 
#' @param type Specifies the type of MDS model used (for more details see
#'  the package '\code{smacof}' of Mair, De Leeuw, Borg, & Groenen). Default is
#' \code{ordinal}.
#' @param weight If \code{TRUE}, the MDS model is calcualted using weights, 
#' i.e., similarities of zero are excluded. Default is \code{TRUE}. 
#' @return Matrix showing the categories included and the Stress values of the 
#' respective MDS models.  
#' @seealso \code{\link[smacof]{smacofSym}} for details on calculating MDS 
#'  representations, \code{\link[thectar]{simi}} for details on calculating  
#'  similarity matrices.
#' @export
#' @examples
#' ## Calculate Stress values for 5 out of 7 SDGs
#' data(SDG_coocurrence)
#' SDG_coocurrence <- SDG_coocurrence[,-2] # Drop second column
#' input <- SDG_coocurrence[,1:8] # For computational reasons, we will work 
#'                                # with 7 SDGs.  
#' stress <- lscomb(input, points = 5)
#' stress


lscomb <- function(data, points, method = "as", single = TRUE, comments = FALSE, type = "ordinal", weight = TRUE){
  
  if(any(is.na(data))){

    stop("The dataset contains missing values.")
    
  }else{
    
    if(!(c(points %% 1) == 0 && points > 0 && points < nrow(data))){
      
      stop("The number of points chosen is either bigger than or equal to the number of objects in data or not a natural number.")
      
    }else{
      
      if(!(method %in% c("as",  "jaccard", "cosine",  "inclusion"))){
        
        method <- "as"
        warning("The parameter 'method' is not properly defined. Default will be used.")
        
      } #if
      
      # Prepare the basic variables 
      
      list <- c(1 : c(length(data) - 1)) ### list with all possible elements
      
      comb <- utils::combn(list, points) ### all possible combinations
      
      combination <- matrix(c(data[, 1]), ncol = 1)
      colnames(combination) <- c(colnames(data[1]))
      output <- matrix(c(rep(0,points + 1)), nrow = 1)
      labels <- c(colnames(data[1]))
      
      # Give a warning to the user that the process may take some time 
      
      ncomb <- ncol(comb) ### number of different combinations
      message(paste0("For ", points, " out of ", c(length(data) - 1), " variables, there are ", 
                  ncomb, " different combinations. It may take a while to calculate the Stress values for all of them.")) 
      
      # Combute the Stress for all possible combinations and save the combination and the Stress values 
      
      i <- 1
      j <- 1
      
      while(i < ncol(comb) + 1){
        
        ### Compute the new matrix 
        while(j < nrow(comb) + 1){
          combination <- cbind(combination, data[, comb[j,i] + 1])
          labels <- cbind(labels, colnames(data[comb[j, i] + 1]))
          j <- j + 1
        } #while
        
        ### Calculate model and save combination and Stress
        colnames(combination) <- labels
        combination <- as.data.frame(combination)
        similarity <- simi(combination, method = method, single = single, comments = comments)
        
        if(weight == TRUE){
          w <- as.numeric(similarity < 1e-8)
          weight_mat = matrix(1 - w, nrow = nrow(similarity), ncol = ncol(similarity))
          dissimilarities <- 1 - similarity
          res <- smacof::smacofSym(dissimilarities, type = type, weightmat = weight_mat)
        }else{
          dissimilarities <- 1 - similarity
          res <- smacof::smacofSym(dissimilarities, type = type)
        } #ifelse
        
        C1 <- matrix(c(colnames(combination[2:c(nrow(comb) + 1)])), nrow = 1) # what was combined
        C2 <- res$stress #stress of the solution
        C3 <- cbind(C1, C2) #what was combined and the stress of the solution
        output <- rbind(output, C3) #add this to output
        
        i <- i + 1
        j <- 1
        
        combination <- data[, 1]
        
        labels <- c(colnames(data[1]))
      } #while
      
      # Order the matrix such that the lowest Stress level is at the top 
      
      output <- output[-1, ]
      output <- as.data.frame(output)
      output <- output[order(output[, points+1]), ]  
      
      # Add labels to output 
      
      i <- 2
      labels <- c("Included object 1")
      
      while(i < points + 1){
        next_label <- paste("Included object", i)
        labels <- cbind(labels, next_label)
        i <- i + 1
      } #while
      
      labels <- cbind(labels, "Resulting Stress value")
      colnames(output) <- labels
      
    } #ifelse
    
  } #ifelse
  
  # Output 
  
  return(output)
  
} #function


