########################################################################################################
##' Clinical and radiological predictors of bacteriologically confirmed TB in young children
##' Author: Jonathan Smith, PhD
##' Last updated: August 2, 2022
########################################################################################################
# Clear everything (if desired)
# rm(list=ls()) 

##' -------------------------------------------
##' MODEL TESTING, TABLES, AND FIGURES
##' -------------------------------------------
##' Step 0: Import models and bring in testing data
##' Step 1: Test Models
##' Step 2: Tables and Figures

########################################################################################################
####' Step 0: Programmatic and Data Prep
####'      0.1 - Load libraries and import models
########################################################################################################
pacman::p_load(caTools, ROSE, MASS, glmnet, 
               caret, randomForest, xgboost, 
               data.table, e1071, pROC, PRROC,
               psych, mltools, ggplot2, tidyr, dplyr) 

modelpath <- "~jonathansmith/Dropbox/CDC/Kenya/Clinical Signs and Symptoms/New Models/"

# Read in Models
mdls <- lapply(list.files(paste0(modelpath, "Models/"), pattern = "*.rds", full.names = TRUE), readRDS)
names(mdls) <- substr(list.files(paste0(modelpath, "Models/"), pattern = "*.rds", full.names = FALSE), 
                      1, nchar(list.files(paste0(modelpath, "Models/"), pattern = "*.rds", full.names = FALSE)) - 4)
list2env(mdls, envir = .GlobalEnv)

####
# Bring in Data
source("~jonathansmith/Dropbox/CDC/Kenya/Clinical Signs and Symptoms/R Code/Final Coding August 2022/Data Prep - Predicting Microbial Confirmation in Very Young Children.R")

aa <- lapply(train_dataprep[3:4], function(x) mapply(table, x))
train_tab <- do.call(cbind, lapply(aa, function(a) do.call(rbind, lapply(a, as.data.frame))))
train_tab$Invasive_text <- paste0(train_tab$Invasive.Freq, " (", round((train_tab$Invasive.Freq/184)*100,1), ")")
train_tab$Noninvasive_text <- paste0(train_tab$Noninvasive.Freq, " (", round((train_tab$Noninvasive.Freq/183)*100,1), ")")

bb <- lapply(test_data[3:4], function(x) mapply(table, x))
test_tab <- do.call(cbind, lapply(bb, function(b) do.call(rbind, lapply(b, as.data.frame))))
test_tab$Invasive_text <- paste0(test_tab$Invasive.Freq, " (", round((test_tab$Invasive.Freq/78)*100,1), ")")
test_tab$Noninvasive_text <- paste0(test_tab$Noninvasive.Freq, " (", round((test_tab$Noninvasive.Freq/79)*100,1), ")")

cc <- lapply(listdata[3:4], function(x) mapply(table, x))
all_tab <- do.call(cbind, lapply(cc, function(c) do.call(rbind, lapply(c, as.data.frame))))
all_tab$Invasive_text <- paste0(all_tab$Invasive.Freq, " (", round((all_tab$Invasive.Freq/262)*100,1), ")")
all_tab$Noninvasive_text <- paste0(all_tab$Noninvasive.Freq, " (", round((all_tab$Noninvasive.Freq/262)*100,1), ")")

exporttab <- cbind(var = rownames(train_tab),
                   all_inv = all_tab$Invasive_text,
                   inv_train = train_tab$Invasive_text,
                   inv_test = test_tab$Invasive_text,
                   all_non = all_tab$Nonnvasive_text,
                   non_train = train_tab$Noninvasive_text,
                   non_test = test_tab$Noninvasive_text)
write.csv(exporttab, paste0(modelpath, "Tables/All-train-test_breakdown_", Sys.Date(), ".csv"))

write.csv(train_tab, paste0(modelpath, "Tables/Training_breakdown_", Sys.Date(), ".csv"))
write.csv(test_tab, paste0(modelpath, "Tables/Testing_breakdown_", Sys.Date(), ".csv"))



########################################################################################################
####' Step 1: Programmatic and Data Prep
####'      1.1 - Write functions to run models
####'      1.2 - Run Models
########################################################################################################

##' Notes: 
##' Recall = Sens
##' Precision = PPV

## 0.1 - Write Functions
  # define function to calculate F score
  # beta = 1 is F1, where sensitivity and PPV are weighted the same
  # beta = 2 is F2, weights sensitivity 2x higher than PPV, this is likelt what we should use
  # beta = 0.5 is F05, weights PPV 2x more than sensitivity
fscore <- function(PPV, Sens, beta){
  fscore <- (1 + beta ^ 2) * ((PPV * Sens) / ((beta ^ 2 * PPV) + Sens))
  return(fscore)
}


## Define function to test trained models
testmodelresults <- function(testdata, trainedmodel, outcome, sigdig = 2){
  
  trained <- get(trainedmodel)[[outcome]]
  
  if(trainedmodel %in% c("forward.train", "backward.train", "both.train")){
      modpred <- predict(trained,
                         type = "response",
                         newdata = testdata[, -which(names(testdata) == "tb_outcome")])
  } else {
    if(trainedmodel %in% c("ridge.train", "lasso.train", "en.train")){
      x.test <- model.matrix(tb_outcome ~ ., data = testdata)[ , -which(names(testdata) %in% "tb_outcome")]
      modpred <- as.numeric(predict(trained, newx = x.test, type = "response"))
    } else {
      if(trainedmodel %in% c("rf.train")){
        modpred <- predict(trained, 
                           type = "prob", 
                           newdata = testdata[ , - which(names(testdata) %in% "tb_outcome")])[,2]
      } else {
        if(trainedmodel %in% c("xgb.train")){
          xgboost_test <- data.table(testdata, keep.rownames = FALSE)
          sparse_matrix_test <- sparse.model.matrix(tb_outcome ~ . - 1, data = xgboost_test)
          output_vector_test <- xgboost_test[, "tb_outcome"] == 1
          xgbtest <- xgb.DMatrix(data = sparse_matrix_test, label = output_vector_test)
          modpred <- predict(trained, xgbtest)
        }
        else {
          if(trainedmodel %in% c("svmlin.train", "svmply.train", "svmrbf.train")){
            mp <- predict(trained, 
                               probability = TRUE, 
                               newdata = testdata[ , - which(names(testdata) == "tb_outcome")])
            modpred <- attr(mp, "probabilities")[ , 2]
          } else (print("MODEL NOT FOUND"))
        }
      }
    }
  } 
  AUC <- pROC::auc(testdata$tb_outcome, modpred)
  #AUC_ci <- ci.auc(testdata$tb_outcome, modpred)
  
  fg <- modpred[testdata$tb_outcome == 1]
  bg <- modpred[testdata$tb_outcome == 0]
  PRC <- PRROC::pr.curve(scores.class0 = fg, scores.class1 = bg)[[2]]
  
  modroc <- pROC::roc(testdata$tb_outcome, modpred)
  youdbest <- coords(modroc, "b", ret = "threshold", best.method = "youden")
  
  
  # confusion matrix
  cM <- table(test = as.factor(ifelse(modpred > youdbest[[1]], 1, 0)), true = testdata$tb_outcome)

  # Get missclass, sens, spec, ppv, and npv
  Missclass <- round((cM[2, 1] + cM[1, 2]) / sum(cM), sigdig)
  Sens <- round(cM[2, 2] / sum(cM[, 2]), sigdig)
  Spec <- round(cM[1, 1] / sum(cM[, 1]), sigdig)
  PPV <- round(cM[2, 2] / sum(cM[2, ]), sigdig)
  NPV <- round(cM[1, 1] / sum(cM[1, ]), sigdig)
  
  F1 <- round(fscore(Sens, PPV, 1),2)
  F2 <- round(fscore(Sens, PPV, 2),2)
  F05 <- round(fscore(Sens, PPV, 0.5), 2)
  
  kappa <- round(cohen.kappa(cbind(as.factor(ifelse(modpred > youdbest[[1]], 1, 0)), testdata$tb_outcome))$kappa,2)
  mcc <- round(mcc(preds = as.factor(ifelse(modpred > youdbest[[1]], 1, 0)),
                         actuals =  testdata$tb_outcome),2)
  
  # combine results
  results <- cbind(auc =  round(AUC, 2),
                   PRC = round(PRC, 2),#paste0(round(AUC_ci[2],sigdig)," (",round(AUC_ci[1],2),"-",round(AUC_ci[3],2),")"), 
                          Missclass, Sens, Spec, PPV, NPV, 
                          F1, F2, F05,
                          kappa, mcc)
  #names(results) <- c("AUC", "Misclassification", "Sensitivity", "Specificity", "PPV", "NPV", "F1", "F2", "F05", "Kappa", "MCC")
  return(list(modpred, results, cM))
}

# Reorder results
resorder <- c("forward.train", "backward.train", "both.train", 
              "ridge.train", "lasso.train", "en.train",
              "rf.train", "xgb.train", 
              "svmlin.train", "svmply.train", "svmrbf.train")
mdls <- mdls[resorder]


outcomess <- names(train_data)
trainnames <- names(mdls)
results_table <- results_models <- results_cM <-  list()

## Run function on all the data and models
for(ww in outcomess) {
  results_table[[ww]] <- data.frame(matrix(NA, ncol = 12, nrow = 0)) 
  names(results_table[[ww]]) <- c("AUROC", "AUPRC", "Misclassification", "Sensitivity", "Specificity", "PPV", "NPV", "F1", "F2", "F05", "Kappa", "MCC")
  for(xx in trainnames) {
    results_models[[ww]][[xx]] <- testmodelresults(testdata = test_data[[ww]],
                                                   trainedmodel = xx,
                                                   outcome = ww)[[1]]
    results_table[[ww]][xx, ] <- testmodelresults(testdata = test_data[[ww]],
                                                  trainedmodel = xx,
                                                  outcome = ww)[[2]]
    results_cM[[ww]][[xx]] <- testmodelresults(testdata = test_data[[ww]],
                                                  trainedmodel = xx,
                                                  outcome = ww)[[3]]
  }
}



########################################################################################################
####' Step 2: Tables and Figures
####'      2.1 - Write tables to CSV/Store Results
####'      2.2 - AUC Figures
####'      2.3 - Beta Figures
########################################################################################################

writetables <- TRUE # Set to TRUE if you want to save results, else will print to the environment

###' - - - - - - - - - - - - - -
###' 2.1 - Write tables to CSV
###' - - - - - - - - - - - - - -

if(writetables){
  for(cc in names(results_table)){
    write.csv(results_table[[cc]], paste0(modelpath, "Tables/", cc, "_", Sys.Date(), ".csv"))
  }
}
write.csv(results_cM[[cc]], paste0(modelpath, "Tables/", cc, "_cM_", Sys.Date(), ".csv"))

## Figure 1
fig1 <- do.call(rbind, results_table[3:4])
fig1[1:11, "type"] <- "invasive"
fig1[12:22, "type"] <- "noninvasive"
xmax <- dim(fig1)[1] * dim(fig1)[2] + dim(fig1)[2]

makeTransparent <- function(someColor, alpha = 100){
  if(alpha < 1) {
    alpha <- alpha * 100
  } else{
    alpha <- alpha
  }
  newColor <- col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

cols <- c("slateblue","darkgreen", "lightcoral", "magenta4")
tcols <- makeTransparent(cols, alpha = 0.50)
tgrey <- makeTransparent("black", alpha = 0.50)
sshift <- 0.15

pdf(paste0(modelpath,"Figures/Metrics_figure_",Sys.Date(),".pdf"), height = 6, width = 10)
par(xpd = TRUE,
    mar = c(2,1,1,1))
plot(NA, xlim = c(0, 15), ylim = c(0, 1), xlab = "", ylab = "", axes = FALSE)
  for(i in 1:(ncol(fig1)-1)){
    points(x = rep(i - sshift, 11), 
           y = fig1[1:11, i], pch = 18, cex = 0.70, col = tcols[1])
    points(x = rep(i + sshift, 11), 
           y = fig1[12:22, i], pch = 18, cex = 0.70, col = tcols[3])
    inviqr <- summary(fig1[1:11, i])[c(2, 3, 5)]
    noniqr <- summary(fig1[12:22, i])[c(2, 3, 5)]
    segments(x0 = i - sshift, y0 = inviqr[1], y1 = inviqr[3])
    segments(x0 = i + sshift, y0 = noniqr[1], y1 = noniqr[3])
    points(x = i - sshift, y = inviqr[2], pch = 23, col = "black", bg = cols[1], cex = 0.85)
    points(x = i + sshift, y = noniqr[2], pch = 23, col = "black", bg = cols[3], cex = 0.85)
    segments(x0 = i + 0.5, y0 = 0, y1 = 1, lty = 3, lwd = 0.75)
  }
segments(x0 = 0.5, y0 = 0, y1 = 1, lty = 1, lwd = 1)
segments(x0 = 0.5, x1 = 13, y0 = 0, lty = 1, lwd = 1)

axis(1, at = 0:13,
     labels = c(NA, colnames(fig1[, 1:2]), "Missclass.\nError", "Sens.", "Spec",
                colnames(fig1[, 6:12]), NA), tick = FALSE,
     cex.axis = 0.70)
segments(x0 = 0.5, x1 = 0.35, y0 = seq(0,1,0.1))
text(c("0.0", seq(0.1,0.9,0.1), "1.0"), x = 00, y = seq(0,1,0.1), cex = 0.9)

# segments(x0 = 0.5, x1 = 12.5, y0 = 0)
# text(c(colnames(fig1[, 1:2]), "Missclass.\nError", "Sens.", "Spec", colnames(fig1[, 6:12])),
#      x = 1:12, y = -0.05, cex = 0.8)

legend(x = 12.5, y = 1, 
       c("Median", "95% CI", NA, 
         "Invasive", "Noninvasive"),
       pch = c(23, NA, NA, 23, 23), 
       lty = c(NA, 1, NA, NA, NA),
       bty = "n", pt.bg = c("black","black", NA, cols[c(1,3)]),
       cex = 0.85)
dev.off()

###' - - - - - - - - - - - - - -
###' Global Figure Settings
###' - - - - - - - - - - - - - -

writefigures <- TRUE # Set to TRUE if you want to save figures, else will print to the environment

# Global Figure Items
allcols <- c("#FF8D7B", "#A81600", "#F5563D", "#711275", "#F26FF7", "#733575",
             "#297527", "#56F551", "#003EA8", "#538FF5", "#B8D2FF")

figlabs <- c("Clinical Diagnosis", 
             "Microbial Confirmation (Any Specimen)",
             "Microbial Confirmation (Invasive Techniques)", 
             "Microbial Confirmation (Noninvasive Techniques)")

###' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###' 2.2 - AUC Figure
###' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
mconly <- TRUE
#writefigures <- FALSE
if(writefigures){
  if(mconly){
    pdf(paste0(modelpath,"Figures/ROC_combined_MConly_",Sys.Date(),".pdf"), height = 6, width = 8)
  } else {
  pdf(paste0(modelpath,"Figures/ROC_combined_",Sys.Date(),".pdf"), height = 8, width = 8)
  }
}

if(mconly){
  par(mfrow = c(1, 2),
      pty = "s") # default)  
} else {
  par(mfrow = c(2, 2),
      pty = "s")
}
if(mconly){
  plotoutcomes <- outcomess[3:4]
} else {
  plotoutcomes <- outcomess
}
for(ww in 1:length(plotoutcomes)){
  plot.roc(smooth(roc(test_data[[plotoutcomes[ww]]]$tb_outcome, 
                      results_models[[plotoutcomes[ww]]][[trainnames[1]]])),
    col = allcols[1],
    lwd = 2,
    #xaxs = "i", #yaxs = "i",
    xlab = "", ylab = "",
    #xlab = "1 - Specificity",
    #ylab = "Sensitivity"
    #legacy.axes = TRUE,
    axes = FALSE, identity = FALSE
  )
  for(i in 2:length(trainnames)) {
    plot.roc(smooth(roc(
      test_data[[plotoutcomes[ww]]]$tb_outcome, results_models[[plotoutcomes[ww]]][[trainnames[i]]])),
    col = allcols[i],
    lwd = 2,
    add = TRUE)
  }
  axis(1, at = seq(0, 1, 0.2), label = rev(c("0.0", seq(0.2, 0.8, 0.2), "1.0")),
       pos = 0, las = 1)
  axis(2, at = seq(0, 1, 0.2),
       pos = 1, las = 1)
  mtext("1 - Specificity", side = 1, padj = 4)
  mtext("Sensitivity", side = 2, padj = -4)
  mtext(paste0(LETTERS[ww], ") "), padj = -4, adj = -0.25, side = 3)
  
  #if(ww == 1){
  legend(x = 0.50, y = 0.75, 
         bty = "n", 
         c("Forward", "Backward", "Bidirectional",
             "Ridge", "LASSO","Elastic Net",
             "Random Forest", "GBT",
             "SVM Linear", "SVM Polynomial", "SVM RBF"),
         lty = 1, lwd = 2, col = allcols, cex = 0.70)
  #}
}
if(writefigures){
  dev.off()
}


###' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###' 2.3 - Betas Figure
###' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Find and Isolates all the betas from the various models
betamodels <- names(mdls)[c(1:6)]
betas <- list()
for(qq in outcomess){
  for(rr in betamodels){
    if(rr %in% betamodels[1:3]){
      x <- data.frame(summary(get(rr)[[qq]])[[12]],
                                      Feature = row.names(summary(get(rr)[[qq]])[[12]]))
      names(x)[4] <- "pz"
      x$Model <- rr
      x$Outcome <- qq
      betas[[qq]][[rr]] <- x[x$pz <= 0.05, c("Estimate", "Feature", "Model", "Outcome")]
    } else { 
      if(rr %in% betamodels[4:6]){
        y <- data.frame(as.matrix(coef(get(rr)[[qq]])),
                                           Feature = row.names(as.matrix(coef(get(rr)[[qq]]))))
        names(y)[1] <- "Estimate"
        y$Model <- rr
        y$Outcome <- qq
        betas[[qq]][[rr]] <- y
      }}
  }
}

# Organize data
betasfig <- Reduce(function(...) merge(..., all = TRUE), unlist(betas, recursive = FALSE))
betasfig <- betasfig[betasfig$Feature != "(Intercept)",]
# betasfig[2:4] <- lapply(betasfig[2:4], as.factor)

# levels(betasfig$Model) <- c("Backward", "Both", "Elastic Net", "Forward", "LASSO", "Ridge")
# 
# levels(betasfig$Feature) <- c("Age 1-2 yrs", "Age 2-3 yrs", "Age 3-4 yrs", "Age 4-5 yrs", 
#                               "Cough", "CXR-TB", "Fever", "Hx Exposure", "HIV+", "IE", "Lethargy",
#                               "Female", "WFA: -2<z<-1", "WFA: z < -2")


#adjj <- 2


############################################################
# writefigures <- FALSE

## Create Figure

ccols <- ggsci::pal_lancet(palette = c("lanonc"), alpha = 1)(6)
yyaxis <- 6
y0aux <- -9.50
y1aux <- -9.25


#writefigures <- TRUE
if(writefigures){
  if(mconly){
    pdf(paste0(modelpath, "Figures/betas_combined_new_MConly_", Sys.Date(), ".pdf"),
        height = 6,
        width = 12)   
  } else {
  pdf(paste0(modelpath, "Figures/betas_combined_new_", Sys.Date(), ".pdf"), height = 8, width = 10)
  }
}
# bottom, left, top and right
# par(mar = c(5.1, 4.1, 4.1, 2.1)) # default
if(mconly){
  par(mfrow = c(1, 2),
      mar = c(8.0, 4.1, 4.1, 2.1),
      xpd = TRUE)
} else {
  par(mfrow = c(2,2),
    mar = c(8.0, 4.1, 4.1, 2.1),
    xpd = TRUE)
}
if(mconly){
  plotoutcomes <- outcomess[3:4]
  plotfiglabs <- figlabs[3:4]
} else {
  plotoutcomes <- outcomess
  plotfiglabs <- figlabs
}
for(pp in 1:length(plotoutcomes)){ 
  zz <- betasfig[betasfig$Outcome == plotoutcomes[pp],]
  aa <- as.matrix(tidyr::pivot_wider(zz, names_from = Feature, values_from = Estimate)[, -c(1:2)])
  aa[is.na(aa)] <- 0
  aa <- aa[,c(9, 1:2, 6:5, 12:11, 13, 4, 14, 10, 7, 8, 3)]

  temp <- barplot(height = aa, beside = TRUE, axes = FALSE,
         names.arg = c("1-2 years", "2-3 years", "3-4 years", "4-5 years", 
                       "Cough", "CXR-TB", "Fever", "Hx. of Exposure", "HIV Positive", "TST/IGRA+", "Lethargy",
                       "Female", "Moderate", "Severe"),
          las = 2, ylim = c(-yyaxis, yyaxis), col = ccols)
  axis(1, tick = TRUE, at = c(0, 100), labels = FALSE)
  axis(2, las = 1)
  segments(x0 = 0, x1 = 100, y0 = 0)
  if(pp == 1){
  legend("topleft", c("Forward","Backward","Both","Ridge","Lasso","Elastic Net"),
         pch = 22, pt.bg = ccols, bty = "n", cex = 0.9, pt.cex = 1.2)
  }
  mtext(expression(paste("Coefficient estimates (", hat(beta), ")")), side = 2, padj = -1.5)
  segments(x0 = 2, x1 = 27, y0 = y0aux)
  segments(x0 = 2, y0 = y0aux, y1 = y1aux)
  segments(x0 = 27, y0 = y0aux, y1 = y1aux)
  text(x = 29/2, y = y0aux - 1, "Age\n(Reference: <1 years)")
  
  segments(x0 = 85, x1 = 99, y0 = y0aux)
  segments(x0 = 85, y0 = y0aux, y1 = y1aux)
  segments(x0 = 99, y0 = y0aux, y1 = y1aux)
  text(x = ((85+99)/2), y = y0aux - 1, "Malnutrition\n(Reference: None)")
  mtext(paste0(LETTERS[pp],") "),  adj = -0.16, padj = -2) # plotfiglabs[pp]),
  
}
if(writefigures){
  dev.off()
}

######################
######################
# Forest plot figure
######################
######################

rfplot <- function(x){
  rfplot <- data.frame(x[[9]])
  rfplot$Feature <- rownames(rfplot)
  return(rfplot[,3:5])
}

rfmodel <- names(mdls)[7]
treebetas <- list()
for(jj in outcomess){
  treebetas[jj] <- mdls[[rfmodel]][jj]
}
rf_accgini <- lapply(treebetas, rfplot)

for(hh in outcomess){
  rf_accgini[[hh]][,"Outcome"] <- hh
}
rf_accgini_long <- lapply(rf_accgini, function(x) {
  tidyr::pivot_longer(x, cols = c(1:2))
  })
rf_accgini_long <- do.call(rbind, rf_accgini_long)

tt <- do.call(rbind, rf_accgini)


## Accuracy
acc <- rf_accgini_long[rf_accgini_long$name == "MeanDecreaseAccuracy",]
accwide <- data.frame(pivot_wider(acc, names_from = Outcome, values_from = value))
rownames(accwide) <- accwide$Feature
accwide <- as.matrix(accwide[ , -c(1:2)])
accwideplot <- accwide[order(rowMeans(accwide)),]

# GINI
gini <- rf_accgini_long[rf_accgini_long$name == "MeanDecreaseGini",]
giniwide <- data.frame(pivot_wider(gini, names_from = Outcome, values_from = value))
rownames(giniwide) <- giniwide$Feature
giniwide <- as.matrix(giniwide[ , -c(1:2)])
giniwideplot <- giniwide[order(rowMeans(giniwide)),]
 

# bottom, left, top and right
# par(mar = c(5.1, 4.1, 4.1, 2.1)) # default
# writefigures <- TRUE
if(mconly){
  accplot <- accwideplot[, 3:4]
  giniplot <- giniwideplot[, 3:4]
} else {
  accplot <- accwideplot
  giniplot <- giniwideplot
}

if(writefigures){
  if(mconly){
    pdf(paste0(modelpath, "Figures/rf_combined_new_MConly_", Sys.Date(), ".pdf"),
        height = 6,
        width = 10)  
  } else {
  pdf(paste0(modelpath, "Figures/rf_combined_new_", Sys.Date(), ".pdf"),
      height = 6,
      width = 12)
  }
}

par(mfrow = c(1,2),
    mar = c(4, 7, 1.75, 1))
barplot(
  t(accplot),
  horiz = TRUE,
  beside = TRUE,
  xlim = c(0, 0.25),
  names.arg = c("Cough", "Lethargy","HIV Positive","Fever","Female",
                "Malnutrition","CXR-TB","Age","TST/IGRA\nPositive","History of\nExposure"),
  las = 1,
  col = ccols[1:ncol(accplot)]
)
segments(x0 = 0, y0 = 1, y1 = 30, lwd = 2)
mtext("A)", adj = -0.25, padj = -1)
mtext("Mean Decrease in Acccuracy", side = 1, padj = 4)
if(mconly){
  legend("bottomright",
         c("Invasive Techniques", "Noninvasive Techniques"),
         pch = 22, pt.bg = ccols[1:2], bty = "n", pt.cex = 1.3)  
} else {
legend("bottomright",
       c("Clinical Diagnosis", "Any Specimen (MC)", "Invasive Techniques (MC)", "Noninvasive Techniques (MC)"),
       pch = 22, pt.bg = ccols[1:4], bty = "n", pt.cex = 1.3)
}
# segments(x0 = 0.09, y0 = 3, y1 = 10)
# text("Microbial\nConfirmation", srt=90, x = 0.08,  y = 6, cex = 0.8)

barplot(
  t(giniplot),
  horiz = TRUE,
  beside = TRUE,
  xlim = c(0, 55),
  names.arg = c("Cough", "Lethargy","HIV Positive","Fever","Female",
                "Malnutrition","CXR-TB","Age","TST/IGRA\nPositive","History of\nExposure"),
  las = 1,
  col = ccols[1:ncol(accplot)]
)
mtext("Mean Decrease in Gini Coefficient", side = 1, padj = 4)
segments(x0 = 0, y0 = 1, y1 = 30, lwd = 2)
mtext("B)", adj = -0.25, padj = -1)
if(writefigures){
 dev.off()
}


##############################
## Gradient Boosted Tree
##############################

xgbmat <- list()
for (ff in outcomess){
  xgbmat[[ff]] <- data.frame(xgb.importance(model = mdls[["xgb.train"]][[ff]]))
  xgbmat[[ff]]$Outcome <- ff
  rownames(xgbmat[[ff]]) <- xgbmat[[ff]]$Feature
  xgbmat[[ff]] <- xgbmat[[ff]][,-c(3:4)]
}
xgb <- do.call(rbind, xgbmat) 

xgbplot <- data.frame(tidyr::pivot_wider(xgb, names_from = Outcome, values_from = Gain))
rownames(xgbplot) <- xgbplot$Feature
xgbplot <- as.matrix(xgbplot[-10,-1])
xgbplot <- xgbplot[order(rowMeans(xgbplot)),]

# bottom, left, top and right
# par(mar = c(5.1, 4.1, 4.1, 2.1)) # default

if(mconly){
  xgbplot_final <- xgbplot[,3:4]
} else {
  xgbplot_final <- xgbplot
}

##'Importance is calculated for a single decision tree 
##'by the amount that each attribute split point improves 
##'the performance measure, weighted by the number of observations the node is responsible for. 

if(writefigures){
  if(mconly){
    pdf(paste0(modelpath, "Figures/xgb_combined_new_MConly_", Sys.Date(), ".pdf"),
        height = 8,
        width = 10)  
  } else {
  pdf(paste0(modelpath, "Figures/xgb_combined_new_", Sys.Date(), ".pdf"),
      height = 6,
      width = 8)
  }
}

par(mar = c(5.1, 7.1, 1, 2.1))
barplot(t(xgbplot_final), beside = TRUE, horiz = TRUE, 
        xlim = c(0, 0.5),
        names.arg = c(
          "HIV Positive",
          "Age 2-3 years",
          "Severe\nMalnutrition",
          "Lethargy",
          "Cough",
          "Age 4-5 years",
          "Age 3-4 years",
          "Age 1-2 years",
          "Fever",
          "Moderate\nMalnutrition",
          "Female", 
          "CXR-TB",
          "TST/IGRA\nPositive",
          "History of\nExposure"
        ), axes = FALSE, las = 2,
        col = ccols[1:ncol(xgbplot_final)])
axis(1)
mtext("Feature Importance", side = 1, padj = 4.5, cex = 1.1)
if(mconly){
  legend("bottomright", c("Invasive Techniques","Noninvasive Techniques"),
         pch = 22, pt.bg = ccols[1:2], bty = "n", pt.cex = 1.3)  
} else {
  legend("bottomright", c("Clinical Diagnosis", "Any Specimen (MC)",
                        "Invasive Techniques (MC)","Noninvasive Techniques (MC)"),
       pch = 22, pt.bg = ccols[1:4], bty = "n", pt.cex = 1.3)
}
if(writefigures){
 dev.off()
}

