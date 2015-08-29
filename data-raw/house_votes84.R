data(HouseVotes84, package = "mlbench")
house_votes84 = (as.matrix(HouseVotes84[, -1]) == "y") * 1.0

rownames(house_votes84) = HouseVotes84$Class
colnames(house_votes84) = c("handicapped-infants", 
                            "water-project-cost-sharing",
                            "adoption-of-the-budget-resolution",
                            "physician-fee-freeze",
                            "el-salvador-aid",
                            "religious-groups-in-schools",
                            "anti-satellite-test-ban",
                            "aid-to-nicaraguan-contras",
                            "mx-missile",
                            "immigration",
                            "synfuels-corporation-cutback",
                            "education-spending",
                            "superfund-right-to-sue",
                            "crime",
                            "duty-free-exports",
                            "export-administration-act-south-africa")

save(house_votes84, file = "data/house_votes84.rdata", compress = "xz")

# library(logisticPCA)
# cv = cv.lpca(votes84_mat, 2, Ms = 1:10)
# plot(cv)
# lpca = logisticPCA(votes84_mat, 2, M = 5)
# 
# plot(lpca, "scores")
# 
# library(ggplot2)
# dfLPCA = data.frame(PC = lpca$PCs, Class = HouseVotes84$Class)
# ggplot(dfLPCA, aes(PC.1, PC.2, colour = Class)) + geom_point()
# 
# 
# lsvd = logisticSVD(votes84_mat, 2)
# 
# plot(lsvd, "scores")
# 
# dfLSVD = data.frame(PC = lsvd$A, Class = HouseVotes84$Class)
# ggplot(dfLSVD, aes(PC.1, PC.2, colour = Class)) + geom_point()
# 
# library(generalizedPCA)
# pca = generalizedPCA(votes84_mat, 2)
# 
# dfPCA = data.frame(PC = pca$PCs, Class = HouseVotes84$Class)
# ggplot(dfPCA, aes(PC.1, PC.2, colour = Class)) + geom_point()
# 
# 
# 
# # LDA ---------------------------------------------------------------------
# library(MASS)
# 
# df_lda = data.frame(Class = HouseVotes84$Class, 
#                     PC = pca$PCs,
#                     LPC = lpca$PCs,
#                     A = lsvd$A
# )
# 
# lda_pca = lda(Class ~ PC.1 + PC.2, data = df_lda)
# lda_lpca = lda(Class ~ LPC.1 + LPC.2, data = df_lda)
# lda_lsvd = lda(Class ~ A.1 + A.2, data = df_lda)
# 
# mean(predict(lda_pca)$class == df_lda$Class)
# mean(predict(lda_lpca)$class == df_lda$Class)
# mean(predict(lda_lsvd)$class == df_lda$Class)
