###############################################################################
# SOCIAL NETWORKS AS A UTILITY MAXIMIZATION DECISION
# Stat. processing by Aldric Labarthe @ ENS Paris-Saclay
# Licensed under GNU GPL v3 - 2O24
###############################################################################

source("utilsGraphs.R")
source("rebuildCompMatrix.R")

###############################################################################
#### Table1 - Similar networks
###############################################################################

morenoTuple <- buildGraphFromStrangeCSV("MorenoTrain/adjMatrixMorenoTrain.csv")
rebuiltM <- rebuildCompMatrix(morenoTuple$adj_mat, 10, 9, 2)
reducedRebuiltM <- reduceCompatibilitiesAtThreshold(morenoTuple$adj_mat, rebuiltM)
pmatrixfile <- "MorenoTrain/clonesFromPDouble/pMatrixSaveMorenoDouble.csv"
pmatrixMoreno <- as.matrix(read.csv(pmatrixfile, header = FALSE))

for (n in 1:100){
  write.table(insertWNinCompMatrix(rebuiltM, 0.1), file=paste("MorenoTrain/sd0.1ClonesNonReduced/CMatrices/", n, ".csv", sep=""),row.names = FALSE,col.names = FALSE, sep = ",")
  write.table(insertWNinCompMatrix(rebuiltM, 0.0001), file=paste("MorenoTrain/sd0.0001ClonesNonReduced/CMatrices/", n, ".csv", sep=""),row.names = FALSE,col.names = FALSE, sep = ",")
}
original <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "",T)
sd0.1ClonesNR <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "MorenoTrain/sd0.1ClonesNonReduced/",T)
sd0.0001ClonesNR <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "MorenoTrain/sd0.0001ClonesNonReduced/",T)
sd1ClonesP <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "MorenoTrain/clonesFromPDouble/sd1/",T)
sd0.2ClonesP <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "MorenoTrain/clonesFromPDouble/sd0.2/",T)
sd0.1ClonesP <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "MorenoTrain/clonesFromPDouble/sd0.1/",T)
sd0.05ClonesP <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "MorenoTrain/clonesFromPDouble/sd0.05/",T)
sd0.0001ClonesP <- importAndCompareSimulatedGraphs(morenoTuple$adj_mat, rebuiltM, "MorenoTrain/clonesFromPDouble/sd0.0001/",T)

meansTable1 <- combine_means(original, sd0.1ClonesNR[id!=0,], sd0.0001ClonesNR[id!=0,], sd1ClonesP[id!=0,], sd0.2ClonesP[id!=0,], sd0.1ClonesP[id!=0,], sd0.05ClonesP[id!=0,], sd0.0001ClonesP[id!=0,])
sdTable1 <- combine_sd(original, sd0.1ClonesNR[id!=0,], sd0.0001ClonesNR[id!=0,], sd1ClonesP[id!=0,], sd0.2ClonesP[id!=0,], sd0.1ClonesP[id!=0,], sd0.05ClonesP[id!=0,], sd0.0001ClonesP[id!=0,])
cat(generate_latex_table(meansTable1, sdTable1, c("id", "mse", "countWeightChanges")))



###############################################################################
#### Invisible hand section
###############################################################################

morenoTuple <- buildGraphFromStrangeCSV("MorenoTrain/adjMatrixMorenoTrain.csv")
for (n in 1:50){
  write.table(buildCompMatrixFromRelationVector((morenoTuple$adj_mat[n,]+morenoTuple$adj_mat[64-n+1,])/2, 0.9)$compMat, file=paste("MorenoTrain/invisibleHand/CMatrices4R0.9/", n, ".csv", sep=""),row.names = FALSE,col.names = FALSE, sep = ",")
}

ihSocialImports <- importSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i4run3/",T, mode="optimal")
ihSocial <- compareSimulatedGraphs(ihSocialImports)
ihSocial <- ihSocial[order(ihSocial$id), ]
ihSocial$utilitySd <- ifelse(ihSocial$utilitySd == 0, 1e-8, ihSocial$utilitySd)
ihSocial$deg_sec_cmom <- ifelse(ihSocial$deg_sec_cmom == 0, 1e-8, ihSocial$deg_sec_cmom)

ihFullScopeI4Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i4run1/",T, mode="local")
ihFullScopeI4Run1 <- ihFullScopeI4Run1[order(ihFullScopeI4Run1$id), ]
ihFullScopeI4Run1Comp <- ihFullScopeI4Run1/ihSocial
ihFullScopeI4Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i4run2/",T, mode="local")
ihFullScopeI4Run2 <- ihFullScopeI4Run2[order(ihFullScopeI4Run2$id), ]
ihFullScopeI4Run2Comp <- ihFullScopeI4Run2/ihSocial
ihFullScopeI4Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i4run3/",T, mode="local")
ihFullScopeI4Run3 <- ihFullScopeI4Run3[order(ihFullScopeI4Run3$id), ]
ihFullScopeI4Run3Comp <- ihFullScopeI4Run3/ihSocial
ihFullScopeI4CompMeans <- (ihFullScopeI4Run1Comp + ihFullScopeI4Run2Comp + ihFullScopeI4Run3Comp) / 3
ihFullScopeI4CompSD <- sqrt(((ihFullScopeI4Run1Comp - ihFullScopeI4CompMeans)^2 + 
                               (ihFullScopeI4Run2Comp - ihFullScopeI4CompMeans)^2 + 
                               (ihFullScopeI4Run3Comp - ihFullScopeI4CompMeans)^2) / 3)
ihFullScopeI4CompMeans$id <- ihFullScopeI4Run1$id
ihFullScopeI4CompSD$id <- ihFullScopeI4Run1$id

ihFullScopeI16Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i16run1/",T, mode="local")
ihFullScopeI16Run1 <- ihFullScopeI16Run1[order(ihFullScopeI16Run1$id), ]
ihFullScopeI16Run1Comp <- ihFullScopeI16Run1/ihSocial
ihFullScopeI16Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i16run2/",T, mode="local")
ihFullScopeI16Run2 <- ihFullScopeI16Run2[order(ihFullScopeI16Run2$id), ]
ihFullScopeI16Run2Comp <- ihFullScopeI16Run2/ihSocial
ihFullScopeI16Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i16run3/",T, mode="local")
ihFullScopeI16Run3 <- ihFullScopeI16Run3[order(ihFullScopeI16Run3$id), ]
ihFullScopeI16Run3Comp <- ihFullScopeI16Run3/ihSocial
ihFullScopeI16CompMeans <- (ihFullScopeI16Run1Comp + ihFullScopeI16Run2Comp + ihFullScopeI16Run3Comp) / 3
ihFullScopeI16CompSD <- sqrt(((ihFullScopeI16Run1Comp - ihFullScopeI16CompMeans)^2 + 
                                (ihFullScopeI16Run2Comp - ihFullScopeI16CompMeans)^2 + 
                                (ihFullScopeI16Run3Comp - ihFullScopeI16CompMeans)^2) / 3)
ihFullScopeI16CompMeans$id <- ihFullScopeI16Run1$id
ihFullScopeI16CompSD$id <- ihFullScopeI16Run1$id

ihFullScopeI32Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i32run1/",T, mode="local")
ihFullScopeI32Run1 <- ihFullScopeI32Run1[order(ihFullScopeI32Run1$id), ]
ihFullScopeI32Run1Comp <- ihFullScopeI32Run1/ihSocial
ihFullScopeI32Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i32run2/",T, mode="local")
ihFullScopeI32Run2 <- ihFullScopeI32Run2[order(ihFullScopeI32Run2$id), ]
ihFullScopeI32Run2Comp <- ihFullScopeI32Run2/ihSocial
ihFullScopeI32Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i32run3/",T, mode="local")
ihFullScopeI32Run3 <- ihFullScopeI32Run3[order(ihFullScopeI32Run3$id), ]
ihFullScopeI32Run3Comp <- ihFullScopeI32Run3/ihSocial
ihFullScopeI32CompMeans <- (ihFullScopeI32Run1Comp + ihFullScopeI32Run2Comp + ihFullScopeI32Run3Comp) / 3
ihFullScopeI32CompSD <- sqrt(((ihFullScopeI32Run1Comp - ihFullScopeI32CompMeans)^2 + 
                                (ihFullScopeI32Run2Comp - ihFullScopeI32CompMeans)^2 + 
                                (ihFullScopeI32Run3Comp - ihFullScopeI32CompMeans)^2) / 3)
ihFullScopeI32CompMeans$id <- ihFullScopeI32Run1$id
ihFullScopeI32CompSD$id <- ihFullScopeI32Run1$id

ihFullScopeI48Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i48run1/",T, mode="local")
ihFullScopeI48Run1 <- ihFullScopeI48Run1[order(ihFullScopeI48Run1$id), ]
ihFullScopeI48Run1Comp <- ihFullScopeI48Run1/ihSocial
ihFullScopeI48Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i48run2/",T, mode="local")
ihFullScopeI48Run2 <- ihFullScopeI48Run2[order(ihFullScopeI48Run2$id), ]
ihFullScopeI48Run2Comp <- ihFullScopeI48Run2/ihSocial
ihFullScopeI48Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i48run3/",T, mode="local")
ihFullScopeI48Run3 <- ihFullScopeI48Run3[order(ihFullScopeI48Run3$id), ]
ihFullScopeI48Run3Comp <- ihFullScopeI48Run3/ihSocial
ihFullScopeI48CompMeans <- (ihFullScopeI48Run1Comp + ihFullScopeI48Run2Comp + ihFullScopeI48Run3Comp) / 3
ihFullScopeI48CompSD <- sqrt(((ihFullScopeI48Run1Comp - ihFullScopeI48CompMeans)^2 + 
                                (ihFullScopeI48Run2Comp - ihFullScopeI48CompMeans)^2 + 
                                (ihFullScopeI48Run3Comp - ihFullScopeI48CompMeans)^2) / 3)
ihFullScopeI48CompMeans$id <- ihFullScopeI48Run1$id
ihFullScopeI48CompSD$id <- ihFullScopeI48Run1$id

ihFullScopeI64Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i64run1/",T, mode="local")
ihFullScopeI64Run1 <- ihFullScopeI64Run1[order(ihFullScopeI64Run1$id), ]
ihFullScopeI64Run1Comp <- ihFullScopeI64Run1/ihSocial
ihFullScopeI64Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i64run2/",T, mode="local")
ihFullScopeI64Run2 <- ihFullScopeI64Run2[order(ihFullScopeI64Run2$id), ]
ihFullScopeI64Run2Comp <- ihFullScopeI64Run2/ihSocial
ihFullScopeI64Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope64i64run3/",T, mode="local")
ihFullScopeI64Run3 <- ihFullScopeI64Run3[order(ihFullScopeI64Run3$id), ]
ihFullScopeI64Run3Comp <- ihFullScopeI64Run3/ihSocial
ihFullScopeI64CompMeans <- (ihFullScopeI64Run1Comp + ihFullScopeI64Run2Comp + ihFullScopeI64Run3Comp) / 3
ihFullScopeI64CompSD <- sqrt(((ihFullScopeI64Run1Comp - ihFullScopeI64CompMeans)^2 + 
                                (ihFullScopeI64Run2Comp - ihFullScopeI64CompMeans)^2 + 
                                (ihFullScopeI64Run3Comp - ihFullScopeI64CompMeans)^2) / 3)
ihFullScopeI64CompMeans$id <- ihFullScopeI64Run1$id
ihFullScopeI64CompSD$id <- ihFullScopeI64Run1$id

fullScopeCompMeans <- combine_means(ihFullScopeI4CompMeans, ihFullScopeI4CompSD, ihFullScopeI16CompMeans, ihFullScopeI16CompSD, ihFullScopeI32CompMeans, ihFullScopeI32CompSD, ihFullScopeI48CompMeans, ihFullScopeI48CompSD, ihFullScopeI64CompMeans, ihFullScopeI64CompSD)
fullScopeCompMeans <- fullScopeCompMeans[c("rels","global_clustering_coefficient", "global_gini_coefficient", "global_density", "mean_distance", "deg_mean", "deg_sec_cmom",
                                           "utilityMean", "utilitySd", "utilityMin", "utilityMax"), ]
fullScopeCompSD <- combine_sd(ihFullScopeI4CompMeans, ihFullScopeI4CompSD, ihFullScopeI16CompMeans, ihFullScopeI16CompSD, ihFullScopeI32CompMeans, ihFullScopeI32CompSD, ihFullScopeI48CompMeans, ihFullScopeI48CompSD, ihFullScopeI64CompMeans, ihFullScopeI64CompSD)
fullScopeCompSD <- fullScopeCompSD[c("rels","global_clustering_coefficient", "global_gini_coefficient", "global_density", "mean_distance", "deg_mean", "deg_sec_cmom",
                                     "utilityMean", "utilitySd", "utilityMin", "utilityMax"), ]

cat(generate_latex_table(fullScopeCompMeans, fullScopeCompSD))

ihScope0Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope0i4run1/",T, mode="local")
ihScope0Run1 <- ihScope0Run1[order(ihScope0Run1$id), ]
ihScope0Run1Comp <- ihScope0Run1/ihSocial
ihScope0Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope0i4run2/",T, mode="local")
ihScope0Run2 <- ihScope0Run2[order(ihScope0Run2$id), ]
ihScope0Run2Comp <- ihScope0Run2/ihSocial
ihScope0Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope0i4run3/",T, mode="local")
ihScope0Run3 <- ihScope0Run3[order(ihScope0Run3$id), ]
ihScope0Run3Comp <- ihScope0Run3/ihSocial
ihScope0CompMeans <- (ihScope0Run1Comp + ihScope0Run2Comp + ihScope0Run3Comp) / 3
ihScope0CompSD <- sqrt(((ihScope0Run1Comp - ihScope0CompMeans)^2 + 
                          (ihScope0Run2Comp - ihScope0CompMeans)^2 + 
                          (ihScope0Run3Comp - ihScope0CompMeans)^2) / 3)
ihScope0CompMeans$id <- ihScope0Run1$id
ihScope0CompSD$id <- ihScope0Run1$id

ihScope1Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope1i4run1/",T, mode="local")
ihScope1Run1 <- ihScope1Run1[order(ihScope1Run1$id), ]
ihScope1Run1Comp <- ihScope1Run1/ihSocial
ihScope1Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope1i4run2/",T, mode="local")
ihScope1Run2 <- ihScope1Run2[order(ihScope1Run2$id), ]
ihScope1Run2Comp <- ihScope1Run2/ihSocial
ihScope1Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope1i4run3/",T, mode="local")
ihScope1Run3 <- ihScope1Run3[order(ihScope1Run3$id), ]
ihScope1Run3Comp <- ihScope1Run3/ihSocial
ihScope1CompMeans <- (ihScope1Run1Comp + ihScope1Run2Comp + ihScope1Run3Comp) / 3
ihScope1CompSD <- sqrt(((ihScope1Run1Comp - ihScope1CompMeans)^2 + 
                          (ihScope1Run2Comp - ihScope1CompMeans)^2 + 
                          (ihScope1Run3Comp - ihScope1CompMeans)^2) / 3)
ihScope1CompMeans$id <- ihScope1Run1$id
ihScope1CompSD$id <- ihScope1Run1$id

ihScope2Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope2i4run1/",T, mode="local")
ihScope2Run1 <- ihScope2Run1[order(ihScope2Run1$id), ]
ihScope2Run1Comp <- ihScope2Run1/ihSocial
ihScope2Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope2i4run2/",T, mode="local")
ihScope2Run2 <- ihScope2Run2[order(ihScope2Run2$id), ]
ihScope2Run2Comp <- ihScope2Run2/ihSocial
ihScope2Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope2i4run3/",T, mode="local")
ihScope2Run3 <- ihScope2Run3[order(ihScope2Run3$id), ]
ihScope2Run3Comp <- ihScope2Run3/ihSocial
ihScope2CompMeans <- (ihScope2Run1Comp + ihScope2Run2Comp + ihScope2Run3Comp) / 3
ihScope2CompSD <- sqrt(((ihScope2Run1Comp - ihScope2CompMeans)^2 + 
                          (ihScope2Run2Comp - ihScope2CompMeans)^2 + 
                          (ihScope2Run3Comp - ihScope2CompMeans)^2) / 3)
ihScope2CompMeans$id <- ihScope2Run1$id
ihScope2CompSD$id <- ihScope2Run1$id

ihScope3Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope3i4run1/",T, mode="local")
ihScope3Run1 <- ihScope3Run1[order(ihScope3Run1$id), ]
ihScope3Run1Comp <- ihScope3Run1/ihSocial
ihScope3Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope3i4run2/",T, mode="local")
ihScope3Run2 <- ihScope3Run2[order(ihScope3Run2$id), ]
ihScope3Run2Comp <- ihScope3Run2/ihSocial
ihScope3Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope3i4run3/",T, mode="local")
ihScope3Run3 <- ihScope3Run3[order(ihScope3Run3$id), ]
ihScope3Run3Comp <- ihScope3Run3/ihSocial
ihScope3CompMeans <- (ihScope3Run1Comp + ihScope3Run2Comp + ihScope3Run3Comp) / 3
ihScope3CompSD <- sqrt(((ihScope3Run1Comp - ihScope3CompMeans)^2 + 
                          (ihScope3Run2Comp - ihScope3CompMeans)^2 + 
                          (ihScope3Run3Comp - ihScope3CompMeans)^2) / 3)
ihScope3CompMeans$id <- ihScope3Run1$id
ihScope3CompSD$id <- ihScope3Run1$id

ihScope4Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope4i4run1/",T, mode="local")
ihScope4Run1 <- ihScope4Run1[order(ihScope4Run1$id), ]
ihScope4Run1Comp <- ihScope4Run1/ihSocial
ihScope4Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope4i4run2/",T, mode="local")
ihScope4Run2 <- ihScope4Run2[order(ihScope4Run2$id), ]
ihScope4Run2Comp <- ihScope4Run2/ihSocial
ihScope4Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.5/scope4i4run3/",T, mode="local")
ihScope4Run3 <- ihScope4Run3[order(ihScope4Run3$id), ]
ihScope4Run3Comp <- ihScope4Run3/ihSocial
ihScope4CompMeans <- (ihScope4Run1Comp + ihScope4Run2Comp + ihScope4Run3Comp) / 3
ihScope4CompSD <- sqrt(((ihScope4Run1Comp - ihScope4CompMeans)^2 + 
                          (ihScope4Run2Comp - ihScope4CompMeans)^2 + 
                          (ihScope4Run3Comp - ihScope4CompMeans)^2) / 3)
ihScope4CompMeans$id <- ihScope4Run1$id
ihScope4CompSD$id <- ihScope4Run1$id

smallScopeCompMeans <- combine_means(ihFullScopeI4CompMeans, ihFullScopeI4CompSD, ihScope0CompMeans, ihScope0CompSD, ihScope1CompMeans, ihScope1CompSD, ihScope2CompMeans, ihScope2CompSD, ihScope3CompMeans, ihScope3CompSD, ihScope4CompMeans, ihScope4CompSD)
smallScopeCompMeans <- smallScopeCompMeans[c("rels","global_clustering_coefficient", "global_gini_coefficient", "global_density", "mean_distance", "deg_mean", "deg_sec_cmom",
                                             "utilityMean", "utilitySd", "utilityMin", "utilityMax"), ]
smallScopeCompSD <- combine_sd(ihFullScopeI4CompMeans, ihFullScopeI4CompSD, ihScope0CompMeans, ihScope0CompSD, ihScope1CompMeans, ihScope1CompSD, ihScope2CompMeans, ihScope2CompSD, ihScope3CompMeans, ihScope3CompSD, ihScope4CompMeans, ihScope4CompSD)
smallScopeCompSD <- smallScopeCompSD[c("rels","global_clustering_coefficient", "global_gini_coefficient", "global_density", "mean_distance", "deg_mean", "deg_sec_cmom",
                                       "utilityMean", "utilitySd", "utilityMin", "utilityMax"), ]

cat(generate_latex_table(smallScopeCompMeans, smallScopeCompSD))

#install.packages("dplyr")
library(dplyr)
library(plm)
ihFullScopeI4CompMeans_tmp <- ihFullScopeI4CompMeans %>% mutate(scope = 0)
ihScope1CompMeans_tmp <- ihScope1CompMeans %>% mutate(scope = 1)
ihScope2CompMeans_tmp <- ihScope2CompMeans %>% mutate(scope = 2)
ihScope3CompMeans_tmp <- ihScope3CompMeans %>% mutate(scope = 3)
ihScope4CompMeans_tmp <- ihScope4CompMeans %>% mutate(scope = 4)

smallScopeDataMerged <- bind_rows(ihFullScopeI4CompMeans_tmp, ihScope1CompMeans_tmp, ihScope2CompMeans_tmp, ihScope3CompMeans_tmp, ihScope4CompMeans_tmp)
panel_data <- pdata.frame(smallScopeDataMerged, index = c("id", "scope"))
clustering_REModel <- plm(global_clustering_coefficient ~ scope, 
                          data = panel_data, 
                          model = "random")
clustering_FEModel <- plm(global_clustering_coefficient ~ scope, 
                          data = panel_data, 
                          model = "within")
phtest(clustering_FEModel,clustering_REModel)

gini_REModel <- plm(global_gini_coefficient ~ scope, 
                    data = panel_data, 
                    model = "random")
distance_REModel <- plm(mean_distance ~ scope, 
                        data = panel_data, 
                        model = "random")
meanDegree_REModel <- plm(deg_mean ~ scope, 
                          data = panel_data, 
                          model = "random")
utility_REModel <- plm(utilityMean ~ scope, 
                       data = panel_data, 
                       model = "random")
sdUtility_REModel <- plm(utilitySd ~ scope, 
                         data = panel_data, 
                         model = "random")
density_REModel <- plm(global_density ~ scope, 
                       data = panel_data, 
                       model = "random")
summary(clustering_REModel)

m1 <- clustering_REModel
m2 <- gini_REModel
m3 <- distance_REModel
m4 <- meanDegree_REModel
m5 <- utility_REModel
m6 <- density_REModel

stargazer(m1, m2, m3, m4, m5, m6,
          type = "latex",
          title = "Random Effects Models for Different Dependent Variables",
          dep.var.labels = c("Global Gini Coefficient", "Mean Distance", "Mean Degree", "Utility Mean", "Global Density"),
          column.labels = c("Gini", "Distance", "Degree", "Utility", "Density"),
          covariate.labels = c("Scope"),
          digits = 3)

ihFullScopeI4R0.9Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope64i4run1/",T, mode="local")
ihFullScopeI4R0.9Run1 <- ihFullScopeI4R0.9Run1[order(ihFullScopeI4R0.9Run1$id), ]
ihFullScopeI4R0.9Run1Comp <- ihFullScopeI4R0.9Run1/ihSocial
ihFullScopeI4R0.9Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope64i4run2/",T, mode="local")
ihFullScopeI4R0.9Run2 <- ihFullScopeI4R0.9Run2[order(ihFullScopeI4R0.9Run2$id), ]
ihFullScopeI4R0.9Run2Comp <- ihFullScopeI4R0.9Run2/ihSocial
ihFullScopeI4R0.9Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope64i4run3/",T, mode="local")
ihFullScopeI4R0.9Run3 <- ihFullScopeI4R0.9Run3[order(ihFullScopeI4R0.9Run3$id), ]
ihFullScopeI4R0.9Run3Comp <- ihFullScopeI4R0.9Run3/ihSocial
ihFullScopeI4R0.9CompMeans <- (ihFullScopeI4R0.9Run1Comp + ihFullScopeI4R0.9Run2Comp + ihFullScopeI4R0.9Run3Comp) / 3
ihFullScopeI4R0.9CompSD <- sqrt(((ihFullScopeI4R0.9Run1Comp - ihFullScopeI4R0.9CompMeans)^2 + 
                                   (ihFullScopeI4R0.9Run2Comp - ihFullScopeI4R0.9CompMeans)^2 + 
                                   (ihFullScopeI4R0.9Run3Comp - ihFullScopeI4R0.9CompMeans)^2) / 3)
ihFullScopeI4R0.9CompMeans$id <- ihFullScopeI4R0.9Run1$id
ihFullScopeI4R0.9CompSD$id <- ihFullScopeI4R0.9Run1$id

ihScope1R0.9Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope1i4run1/",T, mode="local")
ihScope1R0.9Run1 <- ihScope1R0.9Run1[order(ihScope1R0.9Run1$id), ]
ihScope1R0.9Run1Comp <- ihScope1R0.9Run1/ihSocial
ihScope1R0.9Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope1i4run2/",T, mode="local")
ihScope1R0.9Run2 <- ihScope1R0.9Run2[order(ihScope1R0.9Run2$id), ]
ihScope1R0.9Run2Comp <- ihScope1R0.9Run2/ihSocial
ihScope1R0.9Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope1i4run3/",T, mode="local")
ihScope1R0.9Run3 <- ihScope1R0.9Run3[order(ihScope1R0.9Run3$id), ]
ihScope1R0.9Run3Comp <- ihScope1R0.9Run3/ihSocial
ihScope1R0.9CompMeans <- (ihScope1R0.9Run1Comp + ihScope1R0.9Run2Comp + ihScope1R0.9Run3Comp) / 3
ihScope1R0.9CompSD <- sqrt(((ihScope1R0.9Run1Comp - ihScope1R0.9CompMeans)^2 + 
                              (ihScope1R0.9Run2Comp - ihScope1R0.9CompMeans)^2 + 
                              (ihScope1R0.9Run3Comp - ihScope1R0.9CompMeans)^2) / 3)
ihScope1R0.9CompMeans$id <- ihScope1R0.9Run1$id
ihScope1R0.9CompSD$id <- ihScope1R0.9Run1$id

ihScope2R0.9Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope2i4run1/",T, mode="local")
ihScope2R0.9Run1 <- ihScope2R0.9Run1[order(ihScope2R0.9Run1$id), ]
ihScope2R0.9Run1Comp <- ihScope2R0.9Run1/ihSocial
ihScope2R0.9Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope2i4run2/",T, mode="local")
ihScope2R0.9Run2 <- ihScope2R0.9Run2[order(ihScope2R0.9Run2$id), ]
ihScope2R0.9Run2Comp <- ihScope2R0.9Run2/ihSocial
ihScope2R0.9Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope2i4run3/",T, mode="local")
ihScope2R0.9Run3 <- ihScope2R0.9Run3[order(ihScope2R0.9Run3$id), ]
ihScope2R0.9Run3Comp <- ihScope2R0.9Run3/ihSocial
ihScope2R0.9CompMeans <- (ihScope2R0.9Run1Comp + ihScope2R0.9Run2Comp + ihScope2R0.9Run3Comp) / 3
ihScope2R0.9CompSD <- sqrt(((ihScope2R0.9Run1Comp - ihScope2R0.9CompMeans)^2 + 
                              (ihScope2R0.9Run2Comp - ihScope2R0.9CompMeans)^2 + 
                              (ihScope2R0.9Run3Comp - ihScope2R0.9CompMeans)^2) / 3)
ihScope2R0.9CompMeans$id <- ihScope2R0.9Run1$id
ihScope2R0.9CompSD$id <- ihScope2R0.9Run1$id

ihScope3R0.9Run1 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope3i4run1/",T, mode="local")
ihScope3R0.9Run1 <- ihScope3R0.9Run1[order(ihScope3R0.9Run1$id), ]
ihScope3R0.9Run1Comp <- ihScope3R0.9Run1/ihSocial
ihScope3R0.9Run2 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope3i4run2/",T, mode="local")
ihScope3R0.9Run2 <- ihScope3R0.9Run2[order(ihScope3R0.9Run2$id), ]
ihScope3R0.9Run2Comp <- ihScope3R0.9Run2/ihSocial
ihScope3R0.9Run3 <- importAndCompareSimulatedGraphs(NULL, NULL, "MorenoTrain/invisibleHand/FinalData3/R0.9/scope3i4run3/",T, mode="local")
ihScope3R0.9Run3 <- ihScope3R0.9Run3[order(ihScope3R0.9Run3$id), ]
ihScope3R0.9Run3Comp <- ihScope3R0.9Run3/ihSocial
ihScope3R0.9CompMeans <- (ihScope3R0.9Run1Comp + ihScope3R0.9Run2Comp + ihScope3R0.9Run3Comp) / 3
ihScope3R0.9CompSD <- sqrt(((ihScope3R0.9Run1Comp - ihScope3R0.9CompMeans)^2 + 
                              (ihScope3R0.9Run2Comp - ihScope3R0.9CompMeans)^2 + 
                              (ihScope3R0.9Run3Comp - ihScope3R0.9CompMeans)^2) / 3)
ihScope3R0.9CompMeans$id <- ihScope3R0.9Run1$id
ihScope3R0.9CompSD$id <- ihScope3R0.9Run1$id

ihFullScopeI4CompMeans_tmp <- ihFullScopeI4CompMeans %>% mutate(scope = 0, R=0.5)
ihFullScopeI4CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihFullScopeI4CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope0CompMeans_tmp <- ihScope0CompMeans %>% mutate(scope = 1, R=0.5)
ihScope0CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope0CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope1CompMeans_tmp <- ihScope1CompMeans %>% mutate(scope = 2, R=0.5)
ihScope1CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope1CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope2CompMeans_tmp <- ihScope2CompMeans %>% mutate(scope = 3, R=0.5)
ihScope2CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope2CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope3CompMeans_tmp <- ihScope3CompMeans %>% mutate(scope = 4, R=0.5)
ihScope3CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope3CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope4CompMeans_tmp <- ihScope4CompMeans %>% mutate(scope = 5, R=0.5)
ihScope4CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope4CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient

ihFullScopeI4R0.9CompMeans_tmp <- ihFullScopeI4R0.9CompMeans %>% mutate(scope = 0, R=0.9)
ihFullScopeI4R0.9CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihFullScopeI4R0.9CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope1R0.9CompMeans_tmp <- ihScope1R0.9CompMeans %>% mutate(scope = 2, R=0.9)
ihScope1R0.9CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope1R0.9CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope2R0.9CompMeans_tmp <- ihScope2R0.9CompMeans %>% mutate(scope = 3, R=0.9)
ihScope2R0.9CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope2R0.9CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient
ihScope3R0.9CompMeans_tmp <- ihScope3R0.9CompMeans %>% mutate(scope = 4, R=0.9)
ihScope3R0.9CompMeans_tmp[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihScope3R0.9CompMeans_tmp[, "originalGini"] <- ihSocial$global_gini_coefficient

smallScopeWithRDataMerged <- bind_rows(ihFullScopeI4CompMeans_tmp, ihScope0CompMeans_tmp, ihScope1CompMeans_tmp, ihScope2CompMeans_tmp, ihScope3CompMeans_tmp,ihScope4CompMeans_tmp,
                                       ihFullScopeI4R0.9CompMeans_tmp, ihScope1R0.9CompMeans_tmp, ihScope2R0.9CompMeans_tmp, ihScope3R0.9CompMeans_tmp)
panel_data_withR <- pdata.frame(smallScopeWithRDataMerged, index = c("id", "scope", "R"))
clustering_REModel_withR <- plm(global_clustering_coefficient ~ scope + R, 
                                data = panel_data_withR, 
                                model = "random")
clustering_FEModel_withR <- plm(global_clustering_coefficient ~ scope + R, 
                                data = panel_data_withR, 
                                model = "within")
phtest(clustering_FEModel_withR,clustering_REModel_withR)


gini_REModel_withR <- plm(global_gini_coefficient ~ scope + R, 
                          data = panel_data_withR, 
                          model = "random")
gini_FEModel_withR <- plm(global_gini_coefficient ~ scope + R, 
                          data = panel_data_withR, 
                          model = "within")
phtest(gini_REModel_withR,gini_FEModel_withR)


distance_REModel_withR <- plm(mean_distance ~ scope + R, 
                              data = panel_data_withR, 
                              model = "random")
distance_FEModel_withR <- plm(mean_distance ~ scope + R, 
                              data = panel_data_withR, 
                              model = "within")
phtest(distance_REModel_withR,distance_FEModel_withR)

meanDegree_REModel_withR <- plm(deg_mean ~ scope + R, 
                                data = panel_data_withR, 
                                model = "random")
meanDegree_FEModel_withR <- plm(deg_mean ~ scope + R, 
                                data = panel_data_withR, 
                                model = "within")
phtest(meanDegree_REModel_withR,meanDegree_FEModel_withR)

utility_REModel_withR <- plm(utilityMean ~ scope + R, 
                             data = panel_data_withR, 
                             model = "random")
utility_FEModel_withR <- plm(utilityMean ~ scope + R, 
                             data = panel_data_withR, 
                             model = "within")
phtest(utility_REModel_withR,utility_FEModel_withR)

sdUtility_REModel_withR <- plm(utilitySd ~ scope + R, 
                               data = panel_data_withR, 
                               model = "random")
sdUtility_FEModel_withR <- plm(utilitySd ~ scope + R, 
                               data = panel_data_withR, 
                               model = "within")
phtest(sdUtility_REModel_withR,sdUtility_FEModel_withR)

density_REModel_withR <- plm(global_density ~ scope + R, 
                             data = panel_data_withR, 
                             model = "random")
density_FEModel_withR <- plm(global_density ~ scope + R, 
                             data = panel_data_withR, 
                             model = "within")
phtest(density_REModel_withR,density_FEModel_withR)

m1R <- clustering_REModel_withR
m2R <- gini_REModel_withR
m3R <- distance_REModel_withR
m4R <- meanDegree_REModel_withR
m5R <- utility_REModel_withR
m6R <- sdUtility_REModel_withR
m7R <- density_REModel_withR

stargazer(m1R, m2R, m3R, m4R, m5R, m6R, m7R)

clustering_REModel_withRandOld <- plm(global_clustering_coefficient ~ scope + R + scope*originalClustering, 
                                      data = panel_data_withR, 
                                      model = "random")
clustering_FEModel_withRandOld <- plm(global_clustering_coefficient ~ scope + R + scope*originalClustering, 
                                      data = panel_data_withR, 
                                      model = "within")
phtest(clustering_REModel_withRandOld,clustering_FEModel_withRandOld)


distance_REModel_withRandOld <- plm(mean_distance ~ scope + R + scope*originalClustering, 
                                    data = panel_data_withR, 
                                    model = "random")
meanDegree_REModel_withRandOld <- plm(deg_mean ~ scope + R + scope*originalClustering, 
                                      data = panel_data_withR, 
                                      model = "random")
utility_REModel_withRandOld <- plm(utilityMean ~ scope + R + scope*originalClustering, 
                                   data = panel_data_withR, 
                                   model = "random")
sdUtility_REModel_withRandOld <- plm(utilitySd ~ scope + R + scope*originalClustering, 
                                     data = panel_data_withR, 
                                     model = "random")
density_REModel_withRandOld <- plm(global_density ~ scope + R + scope*originalClustering, 
                                   data = panel_data_withR, 
                                   model = "random")
gini_REModel_withRandOld <- plm(global_gini_coefficient ~ scope + R + scope*originalGini, 
                                data = panel_data_withR, 
                                model = "random")

m1Ro <- clustering_REModel_withRandOld
m3Ro <- distance_REModel_withRandOld
m4Ro <- meanDegree_REModel_withRandOld
m5Ro <- utility_REModel_withRandOld
m6Ro <- sdUtility_REModel_withRandOld
m7Ro <- density_REModel_withRandOld
m2Ro <- gini_REModel_withRandOld
stargazer(m1Ro, m3Ro, m4Ro, m5Ro, m6Ro, m7Ro)
stargazer(m2Ro)

ihSocial_tmp3 <- ihSocial
ihSocial_tmp3[,c("global_clustering_coefficient", "global_gini_coefficient", "mean_distance", "deg_mean", "utilityMean", "utilitySd", "global_density")] <- 1
ihSocial_tmp3 <- ihSocial_tmp3 %>% mutate(i=0)
ihSocial_tmp3[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihSocial_tmp3[, "originalGini"] <- ihSocial$global_gini_coefficient
ihFullScopeI4CompMeans_tmp3 <- ihFullScopeI4CompMeans %>% mutate(i=4)
ihFullScopeI4CompMeans_tmp3[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihFullScopeI4CompMeans_tmp3[, "originalGini"] <- ihSocial$global_gini_coefficient
ihFullScopeI16CompMeans_tmp3 <- ihFullScopeI16CompMeans %>% mutate(i=16)
ihFullScopeI16CompMeans_tmp3[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihFullScopeI16CompMeans_tmp3[, "originalGini"] <- ihSocial$global_gini_coefficient
ihFullScopeI32CompMeans_tmp3 <- ihFullScopeI32CompMeans %>% mutate(i=32)
ihFullScopeI32CompMeans_tmp3[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihFullScopeI32CompMeans_tmp3[, "originalGini"] <- ihSocial$global_gini_coefficient
ihFullScopeI48CompMeans_tmp3 <- ihFullScopeI48CompMeans %>% mutate(i=48)
ihFullScopeI48CompMeans_tmp3[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihFullScopeI48CompMeans_tmp3[, "originalGini"] <- ihSocial$global_gini_coefficient
ihFullScopeI64CompMeans_tmp3 <- ihFullScopeI64CompMeans %>% mutate(i=64)
ihFullScopeI64CompMeans_tmp3[, "originalClustering"] <- ihSocial$global_clustering_coefficient
ihFullScopeI64CompMeans_tmp3[, "originalGini"] <- ihSocial$global_gini_coefficient

fullScopeWithRDataMerged <- bind_rows(ihSocial_tmp3, ihFullScopeI4CompMeans_tmp3, ihFullScopeI16CompMeans_tmp3, ihFullScopeI32CompMeans_tmp3, ihFullScopeI48CompMeans_tmp3, ihFullScopeI64CompMeans_tmp3)
panel_data_FS <- pdata.frame(fullScopeWithRDataMerged, index = c("id", "i"))
clustering_REModel_FS <- plm(global_clustering_coefficient ~ i, 
                             data = panel_data_FS, 
                             model = "random")
clustering_FEModel_FS <- plm(global_clustering_coefficient ~ i, 
                             data = panel_data_FS, 
                             model = "within")
phtest(clustering_REModel_FS,clustering_FEModel_FS)

gini_REModel_FS <- plm(global_gini_coefficient ~ i, 
                       data = panel_data_FS, 
                       model = "random")
distance_REModel_FS <- plm(mean_distance ~ i, 
                           data = panel_data_FS, 
                           model = "random")
meanDegree_REModel_FS <- plm(deg_mean ~ i, 
                             data = panel_data_FS, 
                             model = "random")
utility_REModel_FS <- plm(utilityMean ~ i, 
                          data = panel_data_FS, 
                          model = "random")
sdUtility_REModel_FS <- plm(utilitySd ~ i, 
                            data = panel_data_FS, 
                            model = "random")
density_REModel_FS <- plm(global_density ~ i, 
                          data = panel_data_FS, 
                          model = "random")
m1FS <- clustering_REModel_FS
m2FS <- gini_REModel_FS
m3FS <- distance_REModel_FS
m4FS <- meanDegree_REModel_FS
m5FS <- utility_REModel_FS
m6FS <- sdUtility_REModel_FS
m7FS <- density_REModel_FS
stargazer(m1FS, m2FS, m3FS, m4FS, m5FS, m6FS, m7FS)

