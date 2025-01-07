library(KRLS)
library(Metrics)
library(readxl)
library(mdatools)

# MLDF

fl <-  read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "Fuse",
                 range ='E26:AA46')
ec <-  read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "Fuse",
                  range ='AB26:DW46')
ms <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'av 2-3', range ='B1:P22'))
ms <- ms[-1,]
wpi <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='EL1:EL21'))
cd <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='B1:B21'))
pb <-  data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='C1:C21'))

row.names(wpi) <- row.names(fl)
flpca <-  pca(fl, ncomp = 10, scale = T, center = T)
ecpca <-  pca(ec, ncomp = 10, scale =T, center = T)
mspca <-  pca(ms, ncomp = 10, scale = T, center = T)
ec_scores <- data.frame(ecpca$res$cal$scores[,c(1,2)])
fl_scores <- data.frame(flpca$res$cal$scores[,c(1,2)])
ms_scores <- data.frame(mspca$res$cal$scores[,c(1,2)])

#plot(seq(1,10,1), mspca[["res"]][["cal"]][["cumexpvar"]], type = 'b')

all <- cbind(ec_scores,fl_scores, ms_scores)

all <- data.frame(scale(all))
all <- data.frame(cbind(all, wpi, cd, pb))
colnames(all) <- c('1','2', '3', '4', '5', '6', 'wpi', 'cd', 'pb')
row.names(all) <- seq(1,20,1)
#all <- all[!(row.names(all) %in% c('11','9', '6')),]



#PLS, WPI (3)

m_wpi <- pls(all[-c(11,2), -c(9,8,7)], all[-c(11,2),7], 10, cv = 1, scale = T) 
plot(m_wpi)
show(summary(m_wpi))
# looking for the outliers
plotRegcoeffs(m_wpi, show.labels = T, type = "h")
plotPredictions(m_wpi, show.labels = T)
abline(a = 0, b = 1)


#PLS, Cd
m_wpi <- pls(all[-c(11,5,16,12), -c(1,2,9,8,7)], all[-c(11,5,16,12),8], 10, cv = 1, scale = T)  #delete ec
plot(m_wpi)
show(summary(m_wpi))

# looking for the outliers
plotRegcoeffs(m_wpi, show.labels = T, type = "h")
plotPredictions(m_wpi, show.labels = T)
abline(a = 0, b = 1)


#PLS, Pb
m_pb <- pls(all[-c(2), -c(3,4,9,8,7)], all[-c(2),9], 10, cv = 1, scale = T) 
plot(m_pb)
show(summary(m_pb))

# looking for the outliers
plotRegcoeffs(m_pb, show.labels = T, type = "h")
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)


# PLS, lgaCd
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='AO25:AO45'))
rownames(y) <- seq(1,20, 1)

m_wpi <- pls(all[-c(11,12,5,16), -c(1,2,9,8,7)], y[-c(11,12,5,16),], 10, cv = 1, scale = T) 
plot(m_wpi)
show(summary(m_wpi))

# looking for the outliers
plotRegcoeffs(m_wpi, show.labels = T, type = "h")
plotPredictions(m_wpi, show.labels = T)
abline(a = 0, b = 1)




#KLRS, WPI (3)

fl <-  read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "Fuse",
                  range ='E26:AA46')
ec <-  read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "Fuse",
                  range ='AB26:DW46')
ms <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'av 2-3', range ='B1:P22'))
ms <- ms[-1,]
wpi <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='EL1:EL21'))
cd <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='B1:B21'))
pb <-  data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='C1:C21'))

row.names(wpi) <- row.names(fl)
flpca <-  pca(fl, ncomp = 10, scale = T, center = T)
ecpca <-  pca(ec, ncomp = 10, scale =T, center = T)
mspca <-  pca(ms, ncomp = 10, scale = T, center = T)
ec_scores <- data.frame(ecpca$res$cal$scores[,c(1,2)])
fl_scores <- data.frame(flpca$res$cal$scores[,c(1,2)])
ms_scores <- data.frame(mspca$res$cal$scores[,c(1,2)])

#plot(seq(1,10,1), mspca[["res"]][["cal"]][["cumexpvar"]], type = 'b')

all <- cbind(ec_scores,fl_scores, ms_scores)

all <- data.frame(scale(all))
all <- data.frame(cbind(all, wpi, cd, pb))
colnames(all) <- c('1','2', '3', '4', '5', '6', 'wpi', 'cd', 'pb')
row.names(all) <- seq(1,20,1)

all <- all[!(row.names(all) %in% c('1')),]

# WPI
KRLS <- krls(all[,-c(9,8,7)], all[,7] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y = all[,7], x = all[,-c(9,8,7)], method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(all)-1-3)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ all[,7])
plot(all[,7], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(all[,7], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(all)),
     cex = 0.8, pos = 4, col = "gray")
grid()





# KRLS, Cd


fl <-  read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "Fuse",
                  range ='E26:AA46')
ms <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'av 2-3', range ='B1:P22'))
ms <- ms[-1,]
cd <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='B1:B21'))

flpca <-  pca(fl, ncomp = 10, scale =T, center = T)
mspca <-  pca(ms, ncomp = 10, scale = T, center = T)
fl_scores <- data.frame(flpca$res$cal$scores[,c(1,2)])
ms_scores <- data.frame(mspca$res$cal$scores[,c(1,2)])

#plot(seq(1,10,1), mspca[["res"]][["cal"]][["cumexpvar"]], type = 'b')

all <- cbind(fl_scores, ms_scores)
all <- data.frame(cbind(all, cd))

colnames(all) <- c('3', '4', '5', '6', 'cd')
row.names(all) <- seq(1,20,1)

all <- all[!(row.names(all) %in% c('11', '16', '5', '12')),]

# WPI
KRLS <- krls(all[,-c(-5)], all[,5] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y = all[,5], x = all[,-c(5)], method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(all)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ all[,5])
plot(all[,5], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(all[,5], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(all)),
     cex = 0.8, pos = 4, col = "gray")
grid()







# KRLS, lgaCd


fl <-  read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "Fuse",
                  range ='E26:AA46')
ms <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'av 2-3', range ='B1:P22'))
ms <- ms[-1,]
cd <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='AO25:AO45'))
#lga used

flpca <-  pca(fl, ncomp = 10, scale =T, center = T)
mspca <-  pca(ms, ncomp = 10, scale = T, center = T)
fl_scores <- data.frame(flpca$res$cal$scores[,c(1,2)])
ms_scores <- data.frame(mspca$res$cal$scores[,c(1,2)])

#plot(seq(1,10,1), mspca[["res"]][["cal"]][["cumexpvar"]], type = 'b')

all <- cbind(fl_scores, ms_scores)
all <- data.frame(cbind(all, cd))

colnames(all) <- c('3', '4', '5', '6', 'cd')
row.names(all) <- seq(1,20,1)

all <- all[!(row.names(all) %in% c('16','11','12')),]

KRLS <- krls(all[,-c(-5)], all[,5] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y = all[,5], x = all[,-c(5)], method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(all)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ all[,5])
plot(all[,5], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(all[,5], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(all)),
     cex = 0.8, pos = 4, col = "gray")
grid()




# KRLS, Pb


ec <-  read_excel("~/R/Chemometrics/RawData-China.xlsx",sheet = "Fuse",
                  range ='AB26:DW46')
ms <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'av 2-3', range ='B1:P22'))
ms <- ms[-1,]
pb <-  data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='C1:C21'))

ecpca <-  pca(ec, ncomp = 10, scale =T, center = T)
mspca <-  pca(ms, ncomp = 10, scale = T, center = T)
ec_scores <- data.frame(ecpca$res$cal$scores[,c(1,2)])
ms_scores <- data.frame(mspca$res$cal$scores[,c(1,2)])

#plot(seq(1,10,1), mspca[["res"]][["cal"]][["cumexpvar"]], type = 'b')

all <- cbind(ec_scores, ms_scores)

all <- data.frame(scale(all))
all <- data.frame(cbind(all, pb))
colnames(all) <- c('1','2', '5', '6', 'pb')
row.names(all) <- seq(1,20,1)

all <- all[!(row.names(all) %in% c('2')),]

# WPI
KRLS <- krls(all[,-c(5)], all[,5] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y = all[,5], x = all[,-c(5)], method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(all)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ all[,5])
plot(all[,5], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(all[,5], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(all)),
     cex = 0.8, pos = 4, col = "gray")
grid()



