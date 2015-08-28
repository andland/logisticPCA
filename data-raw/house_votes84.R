data_url = "https://archive.ics.uci.edu/ml/machine-learning-databases/voting-records/house-votes-84.data"
votes = read.csv(data_url)

data(HouseVotes84, package = "mlbench")
house_votes84 = (as.matrix(HouseVotes84[, -1]) == "y") * 1.0
str(house_votes84)
rownames(house_votes84) = HouseVotes84$Class

save(house_votes84, file = "data/house_votes84.rdata", compress = "xz")

library(logisticPCA)
cv = cv.lpca(votes84_mat, 2, Ms = 1:10)
plot(cv)
lpca = logisticPCA(votes84_mat, 2, M = 5)

plot(lpca, "scores")

library(ggplot2)
dfLPCA = data.frame(PC = lpca$PCs, Class = HouseVotes84$Class)
ggplot(dfLPCA, aes(PC.1, PC.2, colour = Class)) + geom_point()


lsvd = logisticSVD(votes84_mat, 2)

plot(lsvd, "scores")

dfLSVD = data.frame(PC = lsvd$A, Class = HouseVotes84$Class)
ggplot(dfLSVD, aes(PC.1, PC.2, colour = Class)) + geom_point()

library(generalizedPCA)
pca = generalizedPCA(votes84_mat, 2)

dfPCA = data.frame(PC = pca$PCs, Class = HouseVotes84$Class)
ggplot(dfPCA, aes(PC.1, PC.2, colour = Class)) + geom_point()



# LDA ---------------------------------------------------------------------
library(MASS)

df_lda = data.frame(Class = HouseVotes84$Class, 
                    PC = pca$PCs,
                    LPC = lpca$PCs,
                    A = lsvd$A
)

lda_pca = lda(Class ~ PC.1 + PC.2, data = df_lda)
lda_lpca = lda(Class ~ LPC.1 + LPC.2, data = df_lda)
lda_lsvd = lda(Class ~ A.1 + A.2, data = df_lda)

mean(predict(lda_pca)$class == df_lda$Class)
mean(predict(lda_lpca)$class == df_lda$Class)
mean(predict(lda_lsvd)$class == df_lda$Class)
