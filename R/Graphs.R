#####################################################################
### This file contains the functions producing graphs
#####################################################################


#####################################################################
### highlall (Highlight Similarities version 1)
#####################################################################

#' Highlight highest similarities overall (highlall) 
#' 
#' \code{highlall} draws the highest similarities (overall) into an MDS graph. 
#' 
#' This function is applicable to an MDS solution computed with the package 
#' '\code{smacof}' (Mair, De Leeuw, Borg, & Groenen) or a set of 
#' coordinates. It adds the \code{quantile} percent of highest similarities, as 
#' indicated by the similarity matrix \code{similarity}, to the plot of the
#'  respective map. The objects must occur in the same order in the results/ 
#'  coodinates and the similarity matrix.
#' @param similarity Similarity matrix. 
#' @param results Results of \code{smacofSym} or a set of coordinates. 
#' @param quantile Percentage of the highest similarities that should be drawn.
#' Default is \code{30}. 
#' @param col Color of the graph. Default is \code{black}. 
#' @param coordinates If \code{TRUE}, the input \code{results} consist of a set 
#' of coordinates. Default is \code{FALSE}. 
#' @param add If \code{TRUE}, the points will be added to the existing plot.
#' Default is \code{FALSE}. 
#' @param xlim Numeric vector giving the x coordinates range.
#' @param cex.labels Numeric value indicating font size of the labels.
#' @param main Title of the plot.
#' @return Nothing.
#' @seealso \code{\link[smacof]{smacofSym}}. 
#' @export 
#' @examples
#' ## Calculating an MDS using the package 'smacof' and showing the 25% 
#' ## highest similarities
#' data(SDG_coocurrence)
#' SDG_coocurrence <- SDG_coocurrence[,-2] # Drop second column
#' similarities <- simi(SDG_coocurrence, method = "as")
#' dissimilarities <- 1 - similarities
#' res <- smacof::smacofSym(dissimilarities, type = "ordinal")
#' highlall(similarities, res, quantile = 25, main = "25% Highest Similarities")

highlall <- function(similarity, results, quantile = 30, col = "black", 
                     coordinates = FALSE, add = FALSE, xlim = NULL, 
                     cex.labels = 0.8, main = NULL){
  
  # Find the cut-off value
  
  similarity_dual <- similarity[lower.tri(similarity, diag = FALSE)]
  quantile <- (100 - quantile) / 100
  quantile <- quantile(similarity_dual, quantile)
  
  # Find the values that are above or equal to the cut-off value 
  
  highest_values <- (similarity >= quantile)
  highest_values[upper.tri(highest_values, diag = TRUE)] = FALSE
  
  # Find the points to which the values belong 
  
  ind <- 1
  i <- 1
  j <- 1
  
  coordinates_high <- matrix(c(rep(0, length(which(highest_values)) * 2)), ncol = 2)
  
  while(i < length(highest_values[, 1]) + 1){
    while(j < length(highest_values[, 1]) + 1){
      if(highest_values[i, j] == TRUE){
        coordinates_high[ind, 1] <- i
        coordinates_high[ind, 2] <- j
        ind <- ind + 1
      } #while
      j <- j + 1
    } #while
    j = 1
    i <- i + 1
  } #while
  
  # Plot the graph and draw the lines 
  
  if(coordinates == FALSE){
    coordinates_points <- results$conf 
  }else{
    coordinates_points <- results
  } #ifelse
  
  if(add == FALSE && coordinates == FALSE){
    graphics::plot(results, pch = 19, asp = 1, col = col, xlim = xlim, label.conf = list(label = TRUE, pos =2, col = col, cex = cex.labels), main = main)
  }else if(add == FALSE && coordinates == TRUE){
    graphics::plot(coordinates_points, asp = 1, pch = 19, col = col, 
                   xlim = xlim, main = main)
    graphics::text(coordinates_points + 0.05, rownames(coordinates_points), 
                   cex = cex.labels)
  }else{
    graphics::points(coordinates_points[,1], coordinates_points[,2], 
                     pch = 19, col = col, xlim = xlim)
    graphics::text(coordinates_points + 0.05, rownames(coordinates_points), 
                   cex = cex.labels)
  } #ifelse
  
  i <- 1
  j <- 1
  
  while(i < length(coordinates_high[, 1]) + 1){
    graphics::segments(coordinates_points[coordinates_high[i, 1], 1], 
                       coordinates_points[coordinates_high[i, 1], 2], 
                       coordinates_points[coordinates_high[i, 2], 1], 
                       coordinates_points[coordinates_high[i, 2], 2], 
                       lwd = c(1.9), col = col)
    i <- i + 1
  } #while
  
} # function


#####################################################################
### highlpoints (Highlight Similarities version 2)
#####################################################################

#' Highlight highest similarities per point (highlpoints) 
#' 
#' \code{highlpoints} draws the highest similarities (per point) into an MDS graph. 
#' 
#' This function is applicable to an MDS solution computed with the package 
#' '\code{smacof}' (Mair, De Leeuw, Borg, & Groenen) or a set of 
#' coordinates. It adds the \code{link} highest similarities per point, as 
#' indicated by the similarity matrix \code{similarity}, to the plot of the 
#' respective map. The links belonging to one point are displayed in 
#' the same color. If there is more than one similarity on the last rank 
#' \code{link}, all will be shown. The objects must occur in the same order in 
#' the results/ coodinates and the similarity matrix.
#' @param similarity Similarity matrix. 
#' @param results Results of \code{smacofSym} or a set of coordinates. 
#' @param links Number of similarities that should be drawn per point. 
#' Default is \code{3}. 
#' @param col Color of the points in the graph. Default is \code{black}. 
#' @param coordinates If \code{TRUE}, the input \code{results}  consist of a set 
#' of coordinates. Default is \code{FALSE}. 
#' @param add If \code{TRUE}, the points will be added to the existing plot. 
#' Default is \code{FALSE}. 
#' @param xlim Numeric vector giving the x coordinates range. 
#' @param cex.labels Numeric value indicating font size of the labels.
#' @param main Title of the plot.
#' @return Nothing.
#' @seealso \code{\link[smacof]{smacofSym}}.
#' @export
#' @examples 
#' ## Calculating an MDS using the package 'smacof' and showing the 3 highest 
#' ## similarities per point
#' data(SDG_coocurrence)
#' SDG_coocurrence <- SDG_coocurrence[,-2] # Drop second column
#' similarities <- simi(SDG_coocurrence, method = "as")
#' dissimilarities <- 1 - similarities
#' res <- smacof::smacofSym(dissimilarities, type = "ordinal")
#' highlpoints(similarities, res, links = 3, 
#' main = "3 Highest Similarities Per Point")

highlpoints <- function(similarity, results, links = 3, col = "black", coordinates = FALSE, add = FALSE, xlim = NULL, cex.labels = 0.8, main = NULL){
  
  
  # Find the values that fit the criteria
  
  similarity[row(similarity) == col(similarity)] <- 0
  highest_values <- matrix(c(rep(FALSE, c(nrow(similarity) * 
                                            ncol(similarity)))), 
                           nrow = nrow(similarity), ncol = ncol(similarity))
  
  i <- 1
  j <- 1
  ind <- 1 
  
  v <- 0
  
  while(i < length(similarity[ ,1]) + 1){
    while(ind < links + 1){
      max <- which.max(similarity[i, ])       
      highest_values[i, max] <- TRUE
      number <- similarity[i, max]
      similarity[i, max] <- c(-1)
      v <- rbind(v, i)
      
      if(ind == links) { ### Is this the last turn of the loop? If so we must check that there is not an other equal number. 
        number_k <- sum(similarity[i, ] == number) ### Is there a same number in the column? How many?
        number_i = 1
        while(number_i < number_k + 1){
          max <- which.max(similarity[i, ]) ### We add these numbers to the list. 
          highest_values[i, max] <- TRUE
          similarity[i,max] <- c(-1)
          number_i <- number_i + 1
          v <- rbind(v, i)
        } #while
      } #if
      
      ind <- ind + 1
      
    } #while
    ind <- 1
    i <- i + 1
  } #while
  
  v <- v[-1]
  
  # Find the points to which the values belong 
  
  ind <- 1
  i <- 1 
  j <- 1
  
  coordinates_high <- matrix(c(rep(0, length(which(highest_values)) * 2)), ncol = 2)
  
  while(i < length(highest_values[, 1]) + 1){
    while(j < length(highest_values[, 1]) + 1){
      if(highest_values[i, j] == TRUE){
        coordinates_high[ind, 1] <- i
        coordinates_high[ind, 2] <- j
        ind <- ind + 1
      } #if
      j <- j + 1
    } #while
    j <- 1
    i <- i + 1
  } #while
  
  # Plot the graph and draw the lines 
  
  if(coordinates == FALSE){
    coordinates_points <- results$conf 
  }else{
    coordinates_points <- results
  } #ifelse
  
  if(add == FALSE && coordinates == FALSE){
    graphics::plot(results, pch = 19, asp = 1, col = col, xlim = xlim, 
                   label.conf = list(label = TRUE, pos =2, col = col, 
                                     cex = cex.labels), main = main) 
  }else if(add == FALSE && coordinates == TRUE){
    graphics::plot(coordinates_points, asp = 1, pch = 19, col = col, 
                   xlim = xlim, main = main)
    graphics::text(coordinates_points + 0.05, rownames(coordinates_points), 
                   cex = cex.labels)
    }else{
    graphics::points(coordinates_points, pch = 19, col = col, xlim = xlim)
    graphics::text(coordinates_points + 0.05, rownames(coordinates_points), 
                     cex = cex.labels)
  } #ifelseelse
  
  i <- 1
  j <- 1
  
  while(i < length(coordinates_high[, 1]) + 1){
    graphics::segments(coordinates_points[coordinates_high[i, 1], 1] + 
               c((max(coordinates_points[, 1]) - 
                    min(coordinates_points[, 1])) / 186), 
             coordinates_points[coordinates_high[i, 1], 2] + 
               c((max(coordinates_points[, 1]) - 
                    min(coordinates_points[, 1])) / 172), 
             coordinates_points[coordinates_high[i, 2], 1], 
             coordinates_points[coordinates_high[i, 2], 2], 
             col = v[i], lty = 1)
    i <- i + 1
  } #while
  
} #function


