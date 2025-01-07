# Model waters -- new data, December 22, multisensory system
library(readxl)
library(mdatools)
library(Metrics)

#######PLS, cd
# multisensor system
x<- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                          range ='AA1:AO22'))
rownames(x) <- seq(0,20, 1)
y<- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                          range ='EM1:EM22'))
rownames(y) <- seq(0,20, 1)
x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '15', '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15', '12')),]

m_cd <- pls(x[,-c(1, 8,9,12,11)], y[], 10, cv = 1, scale = T) 
plot(m_cd)
show(summary(m_cd))
plotYCumVariance(m_cd, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_cd, type  = 'h', show.labels = TRUE)

# looking for the outliers
plotRegcoeffs(m_cd, show.labels = T)
plotPredictions(m_cd, show.labels = T)
abline(a = 0, b = 1)

###### KRLS, cd, MS
library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='AA1:AO22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EM1:EM22'))
rownames(y) <- seq(0,20, 1)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('0','2', '5', '11', '15', '12')),]
data <- data[,-c(8,9,11,12, 3, 13, 2, 14, 4)]

KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(0.1, 0.02,0.03,0.04, 0.001, 0.01, KRLS[["lambda"]]), .sigma = (length(data)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ data[,1])
plot(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(data)),
     cex = 0.8, pos = 4, col = "gray")
grid()


################################################################################
########### data fusion PLS, low-level

x  <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='D1:AO22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EM1:EM22'))
rownames(y) <- seq(0,20, 1)

x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '15' , '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15' , '12')),]


m_1 <- pls(x[,-c(24, 31, 32,34)],y, 10, cv = 1, scale = T, center = T)
#plot(m_1)
show(summary(m_1))
plotRMSE(m_1)
plotYCumVariance(m_1)
#look for the outliers
plotPredictions(m_1, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m_1, show.labels = T, type = 'b')

#### LLDF, Cd
library(caret)
library(KRLS)
x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='D1:AO22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EM1:EM22'))
rownames(y) <- seq(0,20, 1)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('0','2', '5', '11', '15', '12')),]

data <- data[,-c(25,30, 29,36,39, 31)]
KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(0.1, 0.02,0.03,0.04, 0.001, 0.01, KRLS[["lambda"]]), .sigma = (length(data)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ data[,1])
plot(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(data)),
     cex = 0.8, pos = 4, col = "gray")
grid()


################################################################################
#### mid-level, PLS

fl <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='D1:Z22'))
rownames(fl) <- seq(0,20, 1)
ms <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='AA1:AO22'))
rownames(ms) <- seq(0,20, 1)
cd <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='EM1:EM22'))
rownames(cd) <- seq(0,20, 1)
fl <-  fl[!(row.names(fl) %in% c('0','2', '5', '11', '15', '12')),]
ms <-  ms[!(row.names(ms) %in% c('0','2', '5', '11', '15', '12')),]
cd <-  cd[!(row.names(cd) %in% c('0','2', '5', '11', '15', '12')),]

flpca <-  pca(fl, ncomp = 1, scale = T, center = T)
plotCumVariance(flpca, show.labels = T, lab.cex = 0.8)

mspca <-  pca(ms, ncomp = 1, scale = T, center = T)
plotCumVariance(mspca, show.labels = T, lab.cex = 0.8)

fl_scores <- data.frame(flpca$res$cal$scores)
ms_scores <- data.frame(mspca$res$cal$scores)

all <- cbind(fl_scores, ms_scores)

m_cd <-  pls(all, cd, 5, cv = 1, scale = T, center = T)
plotRMSE(m_cd)
show(summary(m_cd))
plotPredictions(m_cd, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m_cd, show.labels = T, type = 'h', lab.cex = 0.8)

###############################################################################
################################################################################
#### cd, mid-level KRLS

d <- cbind(cd, fl_scores, ms_scores)

colnames(d) <- c('Cd','1','2')

#d <- d[, -c(3)]

KRLS <- krls(d[,-1], d[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]

test_reg_loo_model <- train(Cd ~ ., data = d, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda =c(KRLS[["lambda"]], 0.001, 0.01,0.1348348, 0.05,0.065, 0.5,0.1,1, 0.2), .sigma = (length(d)-1)))
#0.01,0.1348348, 0.05,0.065, 0.5,0.1,1
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])
rmse_cal <- rmse(test_reg_loo_model[["finalModel"]][["fitted"]], d[,1])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ d[,1])

plot(d[,1], test_reg_loo_model[["pred"]][["pred"]], main = 'KRLS after removing outliers', pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(d[,1], test_reg_loo_model[["pred"]][["pred"]], labels = (row.names(d)), 
     cex = 0.6, pos = 4, col = "gray")
grid()





######################################################PLS, pb
# multisensor system
x<- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                          range ='AA1:AO22'))
rownames(x) <- seq(0,20, 1)
y<- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                          range ='EN1:EN22'))
rownames(y) <- seq(0,20, 1)
x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '15', '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15', '12')),]

m_pb <- pls(x[,], y[], 10, cv = 1, scale = T)
plot(m_pb)
show(summary(m_pb))
plotYCumVariance(m_pb, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_pb, type  = 'h', show.labels = TRUE)

# look for the outliers
plotRegcoeffs(m_pb, show.labels = T)
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)

################################################################################
####### PLS, opt., Cd

fl  <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                             range ='D1:Z22'))
rownames(fl) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EO1:EO22'))
rownames(y) <- seq(0,20, 1)

fl <-  fl[!(row.names(fl) %in% c('0','2', '5', '11', '15' , '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15' , '12')),]


m_1 <- pls(fl,y, 10, cv = 1, scale = T, center = T)
#plot(m_1)
show(summary(m_1))
plotRMSE(m_1)
plotYCumVariance(m_1)
#look for the outliers
plotPredictions(m_1, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m_1, show.labels = T)

###### KRLS, cd
library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='D1:Z22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EO1:EO22'))
rownames(y) <- seq(0,20, 1)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('0','2', '5', '11', '15', '12')),]

KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(0.1, 0.02,0.03,0.04, 0.001, 0.01, KRLS[["lambda"]]), .sigma = (length(data)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ data[,1])
plot(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(data)),
     cex = 0.8, pos = 4, col = "gray")
grid()











######################################################################Pb
library(readxl)
library(mdatools)
library(Metrics)

###############################################################################
# multisensor system, pls
x<- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                          range ='AA1:AO22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EN1:EN22'))
rownames(y) <- seq(0,20, 1)
x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '15', '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15', '12')),]

m_pb <- pls(x[,], y, 10, cv = 1, scale = T, center = T)
plot(m_pb)
show(summary(m_pb))
plotYCumVariance(m_pb, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_pb, type  = 'h', show.labels = TRUE)

# look for the outliers
plotRegcoeffs(m_pb, show.labels = T)
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)



library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='AA1:AO22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EN1:EN22'))
rownames(y) <- seq(0,20, 1)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('0','2', '5', '11', '15', '12')),]

data <- data[,-c(10, 12, 13, 2, 9)]
KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(0.1, 0.02,0.03,0.04, 0.001, 0.01, KRLS[["lambda"]]), .sigma = (length(data)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ data[,1])
plot(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(data)),
     cex = 0.8, pos = 4, col = "gray")
grid()



###############################################################################
# electrochemical sensor, pls 

ec <- x<- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                                range ='AP1:EK22'))
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EP1:EP22'))

row.names(ec) <- seq(0,20,1)
row.names(y) <- seq(0,20,1)

#ec <-  ec[,-c(1,2)] #remove names and Cb
#ec <-  ec[!(row.names(ec) %in% c('0')),]
y <-  y[!(row.names(y) %in% c('0', '2', '5', '11', '12' , '15')),]

ec <-  ec[!(row.names(ec) %in% c('0','2', '5', '11', '12' , '15')),]

#new model
m <- pls(ec, y, 10, cv = 1, scale = T)
#plot(m, ncomp = 5)
plotRMSE(m)
plotYCumVariance(m)
show(summary(m))
#, ncomp = 7
#look for the outliers
plotPredictions(m, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)


################################################################################
#### KRLS, Elect sensor
library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='AP1:EK22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EP1:EP22'))
rownames(y) <- seq(0,20, 1)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('0','2', '5', '11', '15', '12')),]


KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(0.1, 0.02,0.03,0.04, 0.001, 0.01, KRLS[["lambda"]]), .sigma = (length(data)-1)))
show(test_reg_loo_model[["results"]])
barplot(test_reg_loo_model[["finalModel"]][["avgderivatives"]])

st <- lm(test_reg_loo_model[["pred"]][["pred"]] ~ data[,1])
plot(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     main = 'KRLS',
     pch = 16)
abline(0, 1)
abline(st, col = "Blue")
text(data[,1], test_reg_loo_model[["pred"]][["pred"]],
     labels = (row.names(data)),
     cex = 0.8, pos = 4, col = "gray")
grid()


