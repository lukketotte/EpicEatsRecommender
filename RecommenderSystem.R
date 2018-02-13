#### DATA
data <- read.csv("test.csv",
                 colClasses = c("character", "character",
                                "integer", "character",
                                "logical", "logical",
                                "logical", "logical",
                                "logical", "logical"))

# flavour


# rating
ratingmat <- dcast(data, user_phone ~ place_id, value.var = "score", mean)
ratingmat[is.na(ratingmat)] = 0
ratingmat[ratingmat == 0] <- NA
ratingmat2 <- as.data.frame(ratingmat[, -1])
#### ALICE

fSVD <- function(x, k = 40, gamma = 0.01, lambda = 0.01,
                 min_improvement = 1e-6, min_epochs = 100, max_epochs = 500) {
  
  x <- as(x, "matrix")
  
  if (ncol(x) < k || nrow(x) < k)
    stop("k is too large")
  
  # initilize the user-feature and item-feature matrix
  U <- matrix(0.1, nrow = nrow(x), ncol = k)
  V <- matrix(0.1, nrow = ncol(x), ncol = k)
  
  #list of indices pointing to ratings on each item
  itemIDX <- lapply(1:nrow(x), function(temp) which(!is.na(x[temp, ])))
  #list of indices pointing to ratings on each user
  userIDX <- lapply(1:ncol(x), function(temp) which(!is.na(x[, temp])))
  
  errors <- c()
  global_mean = sum(x, na.rm = T) / sum(!is.na(x))
  
  # go through all features
  for (f in 1:k) {
    
    # convergence check
    last_error <- Inf
    delta_error <- Inf
    epoch <- 0L
    p <- tcrossprod(U, V)
    
    while (epoch < max_epochs && delta_error > min_improvement) {
      
      # update user features
      error <- x - p
      temp_U <- U
      for (j in 1:ncol(x)) {
        delta_Uik <- lambda * (error[userIDX[[j]], j] * V[j, f] -
                                 gamma * U[userIDX[[j]], f])
        U[userIDX[[j]], f] <- U[userIDX[[j]], f] + delta_Uik
      }
      
      # update item features
      for (i in 1:nrow(x)) {
        delta_Vjk <- lambda * (error[i, itemIDX[[i]]] * temp_U[i, f] -
                                 gamma * V[itemIDX[[i]], f])
        V[itemIDX[[i]], f] <- V[itemIDX[[i]], f] + delta_Vjk
      }
      
      ### update error
      p <- tcrossprod(U, V)
      new_error <- sqrt(sum(abs(x - p)^2, na.rm = TRUE)/length(x))
      delta_error <- abs(last_error - new_error)
      
      last_error <- new_error
      epoch <- epoch + 1L
      errors[epoch] <- sqrt(sum(abs(x - p)^2, na.rm = TRUE)/length(x)) # not necessary, same as new_error!
      
    }
    
  }
  
  structure(list(U = U, V = V, errors = errors))
}

fm_rec <- function(data, number){
  
  # data processing and get id iN
  #data[data == 0] <- NA
  ratingmat11 <- dcast(data, user_phone~place_id, value.var = "score", mean)
  ratingmat11[is.na(ratingmat11)] = 0
  ratingmat11[ratingmat11 == 0] <- NA
  ratingmat22 <- as.data.frame(ratingmat11[,-1])
  # get id iN
  # num <- formatC(number, width = 10, format = "d", flag = "0")
  num = as.character(number)
  z <- match(ratingmat11$user_phone, num)
  iN <- which(z == 1)
  
  # predict
  z <- fSVD(ratingmat22)
  p <- tcrossprod(z$U,z$V)
  dimnames(p) <- dimnames(ratingmat22)
  #p[p > 5] <- 5
  #p[p < 1] <- 1
  
  # get recommendations for user iN
  usr <- ifelse(is.na(ratingmat22[iN,]), T, F)
  #usr <- ifelse(ratingmat[iN,]==0, T, F)
  usr_rat  <- p[iN,]
  tmpp <- usr_rat[usr]
  rec_sort <- sort(tmpp, decreasing = T)
  rec_list <- names(rec_sort[1:50]) # return 50 restaurants
  
  # commented line is to return both alice result and user rating matrix with predicted values
  # return(list(alice_res = as.matrix(rec_sort[1:50], rownames = rec_list), user_rat_mat = p))
  return(as.matrix(rec_sort[1:50], rownames = rec_list))
  
}



### From a run of the Alice algorithm asuming the user is 067278181
# Not sure how this part is done in azure
top_rest <- fm_rec(data, "067278181")
user_id <- "067278181"


#### HUGO

library(Rcpp)
library(RcppArmadillo)
library(reshape2)

# source the cpp functions for the heavy lifting
sourceCpp("cosine_sim.cpp")

user_flavour_mat <- function(data, flavour){
  
  # flavour <- deparse(substitute(flavour))
  ff <- acast(data, user_phone~place_id, value.var = flavour, fun.aggregate = mean)
  ff[is.na(ff)] = 0
  return(ff)
}

cppFunction(depends = "RcppArmadillo",
            'NumericMatrix similar_users_rcpp(arma::mat X, arma::vec y, int n)
            {
            NumericMatrix res(n, 2);
            for(int i = 0; i != n; i++)
            {
            res(i, 1) = arma::dot(y, X.row(i))/(arma::norm(X.row(i))*arma::norm(y));
            res(i, 0) = i;
            }
            
            return(res);
            }
            ')


hugo <- function(user, data, flav, n = nrow(mat)-1, top_rest, n_sim_users = floor(0.25*nrow(mat)), 
                 b = 10, a = 3)
{
  # current setting of params
  # user: string, user phone
  # data: dataframe
  # flav: flavour (not written as character)
  # top_rest: vector of strings, place id's
  # n_sim_users: how many similar users we use to predict prob of 1, 
  #              default = 0.25 of total users rounded down
  # b: how many restaurants to keep
  # a: is the rule for when to use counter. How many non-zeroes are minimum for 
  #    the user
  mat <- user_flavour_mat(data, flav)
  # do this once, save time ffs
  user_row_id <- which(rownames(mat) == user_id)
  top_rest_col_id <- which(colnames(mat) %in% top_rest)
  # matrix for cosine sim function
  mat_sim <- mat[-user_row_id, ]
  # user vector
  user_vec <- mat[user_row_id, ]
  
  # need to keep track of the other users.
  # do this by taking the row indicators
  # of the all other users in original mat
  # ugly solution, typical adolescents
  if(user_row_id == 1) {
    user_ids <- 2:(n+1)
  } else{
    user_ids <- c(1:(user_row_id-1), (user_row_id+1):(n+1))
  }
  
  # similar users
  cos_sim <- similar_users_rcpp(mat_sim, user_vec, n)
  # append user_ids
  cos_sim[, 1] <- user_ids
  # order matrix by cosine similarity
  cos_sim <- cos_sim[order(cos_sim[, 2], decreasing = T), ]
  # only keep specified number of top users
  # Also we don't need cosine similarity any longer
  # so the second column can be dropped from this 
  # point moving forward
  cos_sim <- as.vector(cos_sim[1:n_sim_users, 1]) # 1:n_sim_users row ids
  
  # next step is to subset the sim_mat by 
  # top restaurants and top similar users
  mat_final <- mat_sim[cos_sim, top_rest_col_id]
  
  # now we model the probability of
  # each restaurants being a 1, which is
  # simply means of each column
  result <- apply(mat_final, 2, function(x) mean(x))
  
  if(length(which(result == 0)) >= a) {
    return(head(sort(result, decreasing = T), b))
  } else {
    # if the result we want actually has less than a non-zero probabilities
    # move to global counter
    cat("WARNING: # Non-zero probabilies <",a, "moving to global counter \n\n")
    result_2 <- apply(mat_sim[, top_rest_col_id], 2, function(x) sum(x))
    return(head(sort(result_2, decreasing = T), b))
  }
}
# hugo <- function(user, data, flavour, n = nrow(mat)-1, top_rest, n_sim_users = floor(0.25*nrow(mat)), 
# b = 10, a = 3)
hugo(user = "067278181", data = data, flav = "gf_friends", top_rest = top_rest)
  


