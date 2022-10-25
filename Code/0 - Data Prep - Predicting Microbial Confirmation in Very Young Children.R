### Data prep

library(ROSE); library(DMwR)
# rm(list = ls())


datapath <- "~jonathansmith/Dropbox/CDC/Kenya/Clinical Signs and Symptoms/Data/Analytic Data/"

alldata <- read.csv(paste0(datapath, "allModelData_FINAL_03AUG2022.csv"))
datnames <- c("casedef", "anypos", "goldstandard", "noninvasive")
varnames <- setdiff(names(alldata), datnames)

#lapply(alldata[varnames], table)
# Variables with no or too few outcomes
#' any_cxr_raisedhemidia (0)
#' any_cxr_pneumothorax (0)
#' any_cxr_hilarmedcalc (1)
#' any_cxr_calc (2)
#' any_cxr_cavitation (3)
#' any_cxr_largeaircomp (6)

# rmv <- c("any_cxr_raisedhemidia", "any_cxr_pneumothorax", "any_cxr_hilarmedcalc",
#          "any_cxr_calc", "any_cxr_cavitation", "any_cxr_largeaircomp")
rmv <- varnames[grep("any_cxr", varnames)]
varnames <- setdiff(varnames, rmv)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### 0.2 - Format and prepare data
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

listdata <- vector(mode = "list", length = length(datnames))
out_names <- c("Clinical", "Any", "Invasive", "Noninvasive")
names(listdata) <- out_names 
for(i in seq_along(listdata)){
  index <- c(datnames[i], varnames)
  listdata[[i]] <- alldata[complete.cases(alldata[,index]), index]
  names(listdata[[i]])[1] <- "tb_outcome"
}

#' Reformat clinical case definitions
#'    - Remove confirmed cases (casedef == 1)
#'    - Recode 2 to 1

listdata[[1]] <- listdata[[1]][listdata[[1]]$tb_outcome %in% c(0,2),]
listdata[[1]][listdata[[1]][,1] == 2,1] <- 1

# Convert everything to factor
factor_cols <- 1:ncol(listdata[[1]])
listdata <- lapply(listdata, function(x){
  x[factor_cols] = lapply(x[factor_cols], as.factor)
  x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
### 0.3 - Split into testing and training data sets
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
seeeed <- 050620
set.seed(seeeed)
test_data <- train_dataprep <- train_data <- vector(mode = "list", length = length(listdata))
names(test_data) <- names(train_data) <- names(train_dataprep) <- out_names
for(i in seq_along(train_data)){
  sample <- sample.split(listdata[[i]]$tb_outcome, SplitRatio = 0.70)
  # Test Data 
  test_data[[i]] <- subset(listdata[[i]], sample == FALSE)
  
  # Train Data
  train_dataprep[[i]] <- subset(listdata[[i]], sample == TRUE)
  
  #ROSE
  #train_data[[i]] <- ROSE(tb_outcome ~ ., data = train_dataprep[[i]], seed = seeeed)$data
  
  #SMOTE
  tt <- train_dataprep[[i]]
  tttab <- table(tt$tb_outcome)
  ttperc <- tttab[2]/sum(tttab)
  train_data[[i]] <- SMOTE(tb_outcome ~ ., data = train_dataprep[[i]],
                           perc.over = 100 * ttperc ^ -1, perc.under = 100 + ttperc * 100)
}

