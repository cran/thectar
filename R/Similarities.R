#####################################################################
### This file contains functions useful to calcualte similarities
#####################################################################

####################################################################
### simi (calculate similarities using indices)
####################################################################

#' Similarity matrix (simi) 
#' 
#' \code{simi} calculates a similarity matrix for co-occurrence data. 
#' 
#' This function applies to co-occurrence data. It calculates a similarity 
#' matrix using one of the following indices: Association Strength, Jaccard, 
#' Cosine, or Inclusion (for a detailed discussion see van Eck & Waltman, 2009, 
#' <doi:10.1002/asi.21075>). Additionally, the function can also generate a 
#' sorted, aggregated, or dichotomized version of the input data table. The 
#' first column of the input matrix should contain the ID of the unit of 
#' comparison, and the following columns the categories for which the 
#' similarity is calculated. Lines belonging to the same unit of comparison 
#' (i.e. same ID) will be combined. \code{simi} is particularly suitable for 
#' not sorted, not aggregated, or not dichotomized datasets. For datasets 
#' already sorted, aggregated, and dichotomized, the package \code{proxy} of 
#' Meyer and Buchta offers an alternative to calculate similarity matrices. 
#' \code{simi} does not work with missing data. 
#' @param data Dataset; the first column must be the ID of the unit of 
#' comparison and all other columns must be categories. 
#' @param method Specifies the output, choose between "\code{sort}" (sorted 
#'  version of the data), "\code{aggregate}" (aggregated version of the data), 
#'  "\code{dichotomize}" (dichotomized version of the data), "\code{as}" 
#'  (similarity matrix using Association Strength Index), "\code{jaccard}" 
#'  (similarity matrix using Jaccard  Index), "\code{cosine}" (similarity 
#'  matrix using Cosine Index), and "\code{inclusion}" (similarity matrix using 
#'  Inclusion Index). Default is \code{sort}. 
#' @param single If \code{TRUE}, single mentionings (i.e. one respondent 
#' mentioning just one category) are included. Default is \code{TRUE}. 
#' @param comments If \code{TRUE}, comments relating to exclusion or possible 
#'  exclusion of categories and respondents are displayed. Default is 
#'  \code{TRUE}. 
#' @return Sorted, aggregated, or dichotomized dataset, or similarity matrix. 
#' @seealso \code{\link[proxy]{dist}} from the package '\code{proxy} for 
#' alternative ways to calculate similarity matrices; van Eck and Waltman 
#' (2009, <doi:10.1002/asi.21075>) for a detailed discussion on 
#' similaritiy measues. 
#' @export 
#' @examples
#' ## Calculate similarities using a dichotomized dataset
#' data(SDG_coocurrence)
#' SDG_coocurrence <- SDG_coocurrence[,-2] # Drop second column
#' similarity <- simi(SDG_coocurrence, method = "as", comments = FALSE)
#' head(similarity)


simi <- function (data, method = c("sort", "aggregate", "dichotomize", 
                                   "as", "jaccard", "cosine", "inclusion"), 
                  single = TRUE, 
                  comments = TRUE) 
{
  if(any(is.na(data))){
    
    stop("The dataset contains missing values.")
    
  }else{
    
    if(is.null(method) == FALSE && method %in% c("sort", 
                                                 "aggregate", 
                                                 "dichotomize", 
                                                 "as", 
                                                 "jaccard", 
                                                 "cosine", 
                                                 "inclusion")){
      
      type <- match.arg(method, c("sort", "aggregate", "dichotomize", 
                                  "as", "jaccard", "cosine", "inclusion"), 
                        several.ok = FALSE)
      
    }else{
      
      type <- "sort"
      warning("No valid method chosen. Default will be used.")
      
    } #ifelse
    
    
    # Sort the data
    data <- data[order(data[, 1]), ]
    if(type == "sort"){output <- data} #if
    
    # Aggregate the data, if needed
    if(type != "sort"){
      if(length((unique(data[,1]))) < nrow(data)){
        
        i <- 1
        l <- ncol(data)
        empty_line <- rep(-999, ncol(data))
        data <- rbind(data, empty_line)
        while (i < length(unique(data[,1]))) {
          if (data[i, 1] == data[i + 1, 1]) {
            data[i, ] <- matrix(c(data[i, 1], c(data[i, 2:l] + data[i + 1, 2:l])), nrow = 1)
            data <- data[-(i + 1), ]
          } #if
          else{
            i <- i + 1
          } #ifelse
        } #while
        data <- data[-i, ]
        
      } #if
      if(type == "aggregate"){output <- data}#if
      
      # Get rid of unused categories
      colsums <- colSums(data)
      if (any(colsums == 0)) {
        notused <- which(colsums %in% c(0))
        notused <- sort(notused, decreasing = TRUE)
        for (i in (1:length(notused))) {
          if (comments == TRUE) {
            message(paste0("Not used category. There are no mentionings in category ", 
                        colnames(data)[notused[i]],". It will be excluded from analysis."))
          } #if
          data <- data[, -notused[i]]
        } #for
      } #if
    }#if
    
    # Dichotomize the data, if needed
    if(type == "dichotomize" | type == "as" | type == "jaccard" | 
       type == "cosine" | type == "inclusion"){
      
      if(any(data[,c(2:ncol(data))] > 1)){
        data[, c(2:ncol(data))][data[,(2:ncol(data))] > 1] <- 1
      } #if
      if(type == "dichotomize"){output <- data} #if
      
    }#if
    
    # Get the similarities, if requested
    if (type == "as" | type == "jaccard" | type == "cosine" | type == "inclusion") {
      # Check if any of the rows does not include at least two elements
      # and handle it according to settings. 
      rowsums <- rowSums(data[, c(2:ncol(data))])
      if (any(rowsums < 2)) {
        notused <- which(rowsums %in% c(0, 1))
        notused <- sort(notused, decreasing = TRUE)
        for (i in (1:length(notused))) {
          if (single == FALSE) {
            if (comments == TRUE) {
              message(paste0("Response with only one or no mentioning. There is only one or no mentioning in the row with ",
                          colnames(data)[1]," ", data[notused[i], 1],
                          ". This response will be excluded from analysis."))
            } #if
            data <- data[-notused[i], ]
          } #if
          if (single == TRUE) {
            if (comments == TRUE) {
              message(paste0("Response with only one or no mentioning. There is only one or no mentioning in the row with ", 
                          colnames(data)[1], " ", data[notused[i], 1],
                          ". "))
            } #if
          } #if
        } #for
          } #if
      
      # Get the requested measure
      output <- matrix(rep(0, (ncol(data) - 1)^2), nrow = ncol(data) - 1)
      m <- 2
      w <- m
      while (m < ncol(data) + 1) {
        while (w < ncol(data) + 1) {
          
          si <- sum(data[,m])
          sj <- sum(data[,w])
          cij <- length(which(data[,m] == 1 & data[,w] == 1))
          
          if(type == "as"){
            output[m - 1, w - 1] <- c(cij/(si * sj))
            output[w - 1, m - 1] <- c(cij/(si *sj))
          } #if
          
          if(type == "jaccard"){
            output[m - 1, w - 1] <- c(cij/(si +  sj - cij))
            output[w - 1, m - 1] <- c(cij/(si +  sj - cij))
          } #if
          
          if(type == "cosine"){
            output[m - 1, w - 1] <- c(cij/(sqrt(si *  sj)))
            output[w - 1, m - 1] <- c(cij/(sqrt(si *  sj)))
          } #if
          
          if(type == "inclusion"){
            output[m - 1, w - 1] <- c(cij/(min(si, sj)))
            output[w - 1, m - 1] <- c(cij/(min(si, sj)))  
          } #if
          
          w <- w + 1
        } #while
        w <- m + 1
        m <- m + 1
      } #while
      rownames(output) <- colnames(data)[-1]
      colnames(output) <- colnames(data)[-1]
      } #if
    
  } #ifelse

  return(output)
  
}#function



####################################################################
### count (calculate similarities by counting co-occurrence)
####################################################################


#' Similarity matrix by counting (simicount) 
#' 
#' \code{simicount} calculates a similarity matrix for sorting data. 
#' 
#' This function is applicable to sorting data. It creates a similarity matrix 
#' showing how often two objects were in the same pile. Each line of the 
#' dataset should refer to one sorting. The first column of the input matrix 
#' should contain the ID of the sorting; the following columns refer to the 
#' objects that have been sorted. The allocation of objects to piles is  
#' indicated with numbers; for each line, the objects that were sorted 
#' into the same pile are given the same number (e.g. all objects with a "1" 
#' are in one pile, all objects with a "2" are in one pile, etc.). This 
#' function does not work with missing values.   
#' @param data Dataset; one row represents one sorting, objects in one pile 
#'  must have the same number. 
#' @return Similarity matrix. 
#' @export
#' @examples
#' ## Calculating similarities using sorted data
#' data(SDG_grouping)
#' similarities <- simicount(SDG_grouping)
#' head(similarities)


simicount <- function(data){
  
  if(any(is.na(data))){
    
    stop("The dataset contains missing values.")
    
  }else{
      
    # Prepare basic variables 
    
    length <- (ncol(data) - 1) ^ 2
    
    output <- matrix(rep(0, length), ncol = (ncol(data) - 1))
    
    labels <- colnames(data)
    labels <- labels[-1]
    
    # Compute a similarity matrix by counting co-occurrences 
    
    i <- 1
    j <- 2
    h <- 2
    
    while(i < nrow(data) + 1){
      while(j < ncol(data) + 1){
        while(h < ncol(data) + 1){
          if(data[i, j] == data[i, h]){
            output[j - 1, h - 1] <- output[j - 1, h - 1] + 1
          } #if
          h <- h + 1
        } #while
        j <- j + 1
        h <- 2
      } #while
      i <- i + 1
      j <- 2
      h <- 2
    } #while
    
    colnames(output) <- labels
    rownames(output) <- labels
    
    } #ifelse
  
  # Output 
  
  return(output)

} #function

