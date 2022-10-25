########################################################################################################
##' Clinical and radiological predictors of bacteriologically confirmed TB in young children
##' Author: Jonathan Smith, PhD
##' Last updated: August 2, 2022
########################################################################################################
# Clear everything (if desired)
# rm(list=ls()) 

##' ----------------
##' MODEL TRAINING
##' ----------------
##' Step 0: Prep
##' Step 1: Train Models

########################################################################################################
####' Step 0: Programmatic and Data Prep
####'      0.1 - Load libraries and import data
########################################################################################################

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### 0.1 - Load libraries and import data
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
 # install.packages("pacman") # if not installed
 # Load libraries
pacman::p_load(caTools, ROSE, MASS, glmnet, 
               caret, randomForest, xgboost, 
               data.table, e1071) 

# Run data prep script to import data
source("~jonathansmith/Dropbox/CDC/Kenya/Clinical Signs and Symptoms/R Code/Final Coding August 2022/Data Prep - Predicting Microbial Confirmation in Very Young Children.R")

########################################################################################################
####' Step 1: train_data Models
####'   Models
####'      1.01 - Forward Selection
####'      1.02 - Backwards Selection
####'      1.03 - Bidirectional Selection
####'      1.04 - Ridge
####'      1.05 - Lasso
####'      1.06 - Elastic Net
####'      1.07 - Random Forest
####'      1.08 - Gradient Boosted Tree
####'      1.09 - SVM Linear
####'      1.10 - SVM Polynomial
####'      1.11 - SVM RBF
####'      
####'      1.12 - Save all models
########################################################################################################


#fmla <- as.formula(paste0("~ +", paste(varnames, collapse = " + ")))

#' ________________________________________________
#' 1.01 - Forward Selection
#' ________________________________________________
forward.train <- lapply(train_data, function(x) {
  stepAIC(
    glm(tb_outcome ~ ., family = binomial, data = x),
    # scope = list(lower ~ 1,
    #              upper = fmla),
    direction = "forward",
    trace = 0
  )
})

#' ________________________________________________
#' 1.02 - Backward Selection
#' ________________________________________________
backward.train <- lapply(train_data, function(x) {
  stepAIC(
    glm(tb_outcome ~ ., family = binomial, data = x),
    # scope = list(lower ~ 1,
    #              upper = fmla),
    direction = "backward",
    trace = 0
  )
})

#' ________________________________________________
#' 1.03 - Bidirectional Selection
#' ________________________________________________
both.train <- lapply(train_data, function(x) {
  stepAIC(
    glm(tb_outcome ~ ., family = binomial, data = x),
    # scope = list(lower ~ 1,
    #              upper = fmla),
    direction = "both",
    trace = 0
  )
})

#' ________________________________________________
#' 1.04-1.05 - Ridge and Lasso
#' ________________________________________________
##' Write custom function to run ridge and lasso analysis over list data 
##'    @param modeldata - Analytic data 
##'    @param modeltype - model type, either "ridge" or "lasso"
##'    @param cv - cross validation folds, default to 3

ridge_lasso_func <- function(modeldata, modeltype, cv = 3){ 
  mod <- tolower(modeltype)
  if (mod == "lasso") {
    alpha <- 1
  } else{
    if (mod == "ridge") {
      alpha <- 0
    }
    else{
      print("Please specify model type (modeltype = 'ridge' or 'lasso')")
    }
  }
  cvfolds <- cv
  ## matrix of predictor variables
  x.inv <- model.matrix(tb_outcome ~ .,
                        data = modeldata)[, -which(names(modeldata) %in% "tb_outcome")]
  ## vector of outcome
  y.inv <- modeldata$tb_outcome
  ## Create a vector of lambda values
  lambvec <- 10 ^ seq(0.01, -2, length = 1000)
 
  # Find best value for lambda (cross validation)
  ridgecv <- cv.glmnet(x.inv, y.inv, family = "binomial", 
                       alpha = alpha, nfold = cvfolds)
  best_lam <- ridgecv$lambda.min
  
  final_model <- glmnet(x.inv, y.inv, family = "binomial", 
                        alpha = alpha, lambda = best_lam, standardize = TRUE)
  return(final_model)
}

ridge.train <- lapply(train_data, ridge_lasso_func, modeltype = "ridge")
lasso.train <- lapply(train_data, ridge_lasso_func, modeltype = "lasso")

#' ________________________________________________
#' 1.06 - Elastic Net
#' ________________________________________________
##' Write custom function to run elastic net analysis over list data 
##'    @param modeldata - Analytic data 
##'    @param cv - cross validation folds, default to 3
##'    @param tl - tune length, default to 10

en_func <- function(modeldata, cv = 3, tl = 10){
  x.inv <- model.matrix(tb_outcome ~ ., data = modeldata)[,-which(names(modeldata) %in% "tb_outcome")]
  y.inv <- modeldata$tb_outcome
  tune <- train(tb_outcome ~ ., data = modeldata, method = "glmnet", 
                trControl = trainControl("cv", number = cv),
                tuneLength = tl)
  
  model <- glmnet(x.inv, y.inv, family= "binomial", 
                   alpha = tune$bestTune$alpha,
                   lambda = tune$bestTune$lambda,
                   standardize = TRUE)
  return(model)
}

en.train <- lapply(train_data, en_func)

#' ________________________________________________
#' 1.07 - Random Forest
#' ________________________________________________

rf_func <- function(modeldata, trees = 1000, cvfolds = 3){
  ##' First will create a loop to determine the best number of mtrys to include in the model.
  ##' "mtry" is the number of features/predictor variables that are included in each node.
  
  # create vector for models with mtry values ranging from 1 to to total number of factors
  mtrys <- 1:(ncol(modeldata)-1) # minus one to remove the outcome 
  
  # Divide data into n folds for cross validation
  set.seed(05062020)
  f <- sample(rep(1:cvfolds, length = nrow(modeldata)))
  
  #' The following loop will produce a data frame containing the misclassification 
  #' errors for each fold and for each # of factors. This will determine optimal mtry
  
  # Initiate matrix of misclassification error for the loop
  properr <- data.frame(matrix(NA, ncol = length(unique(f)), nrow = length(mtrys)))
  
  for (i in 1:cvfolds){
    for (j in 1:length(mtrys)){
      #Remove the i fold and get the RF model
      temploop.model <- randomForest(tb_outcome ~ ., data = modeldata[f != i, ], mtry = mtrys[j], ntree = trees, importance = TRUE)
      # Use the remaining fold for prediction
      temploop.pred <- predict(temploop.model, newdata = modeldata[f == i,-which(names(modeldata) == "tb_outcome")])
      # Confusion matrix
      cM <- table(temploop.pred, modeldata[f == i, which(names(modeldata) == "tb_outcome")])
      # Error (misclassification)
      properr[j, i] <- (cM[2, 1] + cM[1, 2]) / (cM[1, 1] + cM[1, 2] + cM[2, 1] + cM[2, 2])
    }
  }
  # Calculate the mean error across all three folds
  properr$meanerr <- rowMeans(properr)
  bestmtry <- which(properr$meanerr == min(properr$meanerr))
  bestmtry <- bestmtry[length(bestmtry)] #largest in case there are >1 with the same min error
  
  rf_model <- randomForest(
    tb_outcome ~ ., data = modeldata,
    mtry = bestmtry,
    ntree = trees, importance = TRUE)
 return(rf_model) 
}

rf.train <- lapply(train_data, rf_func)

#' ________________________________________________
#' 1.08 - Gradient Boosted Tree
#' ________________________________________________

gbt_func <- function(modeldata, cvfolds = 3, nrounds = 1000, stopearly = 20){
  
  ## Prepare dataset for XGBoost package
  xgboost_train <- data.table(modeldata, keep.rownames = FALSE)
  sparse_matrix <- sparse.model.matrix(tb_outcome ~ . - 1, data = xgboost_train)
  output_vector <-  xgboost_train[, "tb_outcome"] == 1

  ### Tune for XGBoost
  # set matrix of parameters to evaluate
  xgb.params.matrix <- expand.grid(eta = c(.01, .05, .1, .3),
                                   max_depth = c(1, 3, 5, 7),
                                   min_child_weight = c(1, 3, 5, 7),
                                   subsample = c(0.65, 0.8, 1), 
                                   colsample_bytree = c(0.8, 0.9, 1),
                                   optimal_trees = 0,
                                   min_RMSE = 0)
  
  for (i in 1:nrow(xgb.params.matrix)){
    tempparams.loop <-  list(eta = xgb.params.matrix$eta[i],
                             max_depth = xgb.params.matrix$max_depth[i],
                             min_child_weight = xgb.params.matrix$min_child_weight[i],
                             subsample = xgb.params.matrix$subsample[i],
                             colsample_bytree = xgb.params.matrix$colsample_bytree[i])
    set.seed(05062020)
    xgbloop.tune <- xgb.cv(params = tempparams.loop,
                           data = sparse_matrix,
                           label = output_vector,
                           nrounds = nrounds,
                           nfold = cvfolds,
                           verbose = 0,
                           early_stopping_rounds = stopearly) # stop if no improvement after 20 consecutive trees
    # add min error and trees to grid
    xgb.params.matrix$optimal_trees[i] <- which.min(xgbloop.tune$evaluation_log$test_rmse_mean)
    xgb.params.matrix$min_RMSE[i] <- min(xgbloop.tune$evaluation_log$test_rmse_mean)
  }

  # Choose model with lowest error
  gb.bestparams <- xgb.params.matrix[order(xgb.params.matrix$min_RMSE)[1],]
    
  # Create list from above for use in XGBoost package
  gb.bestparams.list <-  list(eta = gb.bestparams$eta,
                              max_depth = gb.bestparams$max_depth,
                              min_child_weight = gb.bestparams$min_child_weight,
                              subsample = gb.bestparams$subsample,
                              colsample_bytree = gb.bestparams$colsample_bytree)
  # Train final model
  xgb_model <- xgboost(params = gb.bestparams.list, 
                     data = sparse_matrix,
                     label = output_vector,
                     nrounds = gb.bestparams$optimal_trees,
                     verbose = 0)
  return(xgb_model)
}

xgb.train <- lapply(train_data, gbt_func)

#' ________________________________________________
#' 1.09-1.11 - Support Vector Machines
#'                Linear
#'                Polynomial
#'                RBF
#' ________________________________________________

svm_func <- function(modeldata, kernal, cvfolds = 3){
  krnl <- tolower(kernal)
  svm_modeltune <- tune(svm, tb_outcome ~ ., 
                    data = modeldata, 
                    kernel = krnl, 
                    ranges = list(cost = seq(0.01, 2, 0.05)),
                    tunecontrol = tune.control(cross = cvfolds), 
                    scale = TRUE, 
                    probability = TRUE)
  svm_model <- svm_modeltune$best.model
  return(svm_model)
}

svmlin.train <- lapply(train_data, svm_func, kernal = "linear")
svmply.train <- lapply(train_data, svm_func, kernal = "polynomial")
svmrbf.train <- lapply(train_data, svm_func, kernal = "radial")

#' ________________________________________________
#' 1.12 - Save all models
#' ________________________________________________
modelpath <- "~jonathansmith/Dropbox/CDC/Kenya/Clinical Signs and Symptoms/New Models/Models/"
saveRDS(forward.train, paste0(modelpath,"forward.train.rds"))
saveRDS(backward.train, paste0(modelpath,"backward.train.rds"))
saveRDS(both.train, paste0(modelpath,"both.train.rds"))
saveRDS(ridge.train, paste0(modelpath,"ridge.train.rds"))
saveRDS(lasso.train, paste0(modelpath,"lasso.train.rds"))
saveRDS(en.train, paste0(modelpath,"en.train.rds"))
saveRDS(rf.train, paste0(modelpath,"rf.train.rds"))
saveRDS(xgb.train, paste0(modelpath,"xgb.train.rds"))
saveRDS(svmlin.train, paste0(modelpath,"svmlin.train.rds"))
saveRDS(svmply.train, paste0(modelpath,"svmply.train.rds"))
saveRDS(svmrbf.train, paste0(modelpath,"svmrbf.train.rds"))