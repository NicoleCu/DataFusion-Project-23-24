library(readxl)
library(mdatools)
library(Metrics)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='AA1:AO22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EL1:EL22'))
rownames(y) <- seq(0,20, 1)
x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '15', '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15', '12')),]

###############################################################################
########## MS, PLS
#-c(11, 1, 9, 10, 8,12, 2,6)
m <- pls(x[,-c(11,1,8,9,12, 2,3,6, 10)], y, 10, cv = 1, scale = T)
plot(m)
show(summary(m))
plotPredictions(m, show.labels = T)
abline(a =0,b=1)
plotRegcoeffs(m, show.labels = T)

res_cal <- m[["res"]][["cal"]][["y.pred"]][,3,1]
res_cv <- m[["res"]][["cv"]][["y.pred"]][,3,1]
res_ref <-  m[["res"]][["cv"]][["y.ref"]]
res <-data.frame(cbind(res_cv, res_cal, res_ref))



ggplot(res) +
  #geom_smooth((aes(x = c(3,3,4,5,5.5,5.6,5.7,9.8,9,5,6,7,4,4,4) , y =c(3,3,4,5,5.5,5.6,5.7,9.8,9,5,6,7,4,4,4))), method = "lm", se=FALSE, col = 'black',  alpha = 0.7, size = 1)+
  geom_point(aes(x = res_ref , y = res_cal), size = 3, col = 'red', alpha = 0.7) +
  geom_point(aes(x = res_ref , y = res_cv), size = 3, col = 'blue', alpha = 0.7) +
  geom_smooth((aes(x = res_ref , y = res_cal)), method = "lm", se=FALSE, col = 'red',  alpha = 0.7, size = 1) +
  geom_smooth((aes(x = res_ref , y = res_cv)), method = "lm", se=FALSE, col = 'blue',  alpha = 0.7, size = 1)+
  theme_bw() +
  scale_y_continuous(name='Reference value, ppb')+
  scale_x_continuous(name='Predicted value, ppb')+
  theme(text = element_text(size = 15,face="bold", colour ='black'),
        axis.text.x  = element_text(size=13.5,colour = 'black'),
        axis.text.y = element_text(size=13.5,colour = 'black'),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.major = element_line(colour = "grey", linetype= 'dashed'),
        panel.grid.minor = element_line(colour = "grey", linetype= 'dashed'))
abline()




################################################################################
### E-tongue, KRLS
library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='AA1:AO22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EL1:EL22'))
rownames(y) <- seq(0,20, 1)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('0','2', '5', '11', '15', '12')),]

#data <- data[,-c(11+1,1+1,8+1,9+1,12+1, 2+1,3+1,6+1, 10+1, 16)]
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
########### data fusion PLS, low-level (3)

x  <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='D1:EK22'))
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EL1:EL22'))
rownames(y) <- seq(0,20, 1)

x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '15' , '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15' , '12')),]


m_1 <- pls(x[,-c(1:8, 11+23,1+23,8+23,9+23,12+23, 2+23,3+23,6+23, 10+23)],y, 10, cv = 1, scale = T, center = T)
#plot(m_1)
show(summary(m_1))
plotRMSE(m_1)
plotYCumVariance(m_1)
#look for the outliers
plotPredictions(m_1, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m_1, show.labels = T, type = 'b')

### LLDF, KRLS (3)

x  <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='D1:EK22'))
rownames(x) <- seq(0,20, 1)

y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EL1:EL22'))
rownames(y) <- seq(0,20, 1)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('0','2', '5', '11', '15', '12')),]

#data <- data[,-c(2:8, 11+24,1+24,8+24,9+24,12+24, 2+24,3+24,6+24, 10+24)]
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




library(ggplot)

res_cv <- new_m_n$cvres$y.pred[,1, 1]
res_cal <- new_m_n$calres$y.pred[,1, 1]
res_ref <- new_m_n[["calres"]][["y.ref"]]
res <-data.frame(cbind(res_cv, res_cal, res_ref))


ggplot(res) +
  geom_point(aes(x = V3 , y = res_cal), size = 3, col = 'red', alpha = 0.7) +
  geom_point(aes(x = V3 , y = res_cv), size = 3, col = 'blue', alpha = 0.7) +
  geom_smooth((aes(x = V3 , y = res_cal)), method = "lm", se=FALSE, col = 'red',  alpha = 0.7, size = 1) +
  geom_smooth((aes(x = V3 , y = res_cv)), method = "lm", se=FALSE, col = 'blue',  alpha = 0.7, size = 1)+
  theme_bw() +
  scale_color_manual(name='Regression Model',
                     breaks=c('Linear', 'Quadratic', 'Cubic'),
                     values=c('Cubic'='pink', 'Quadratic'='blue', 'Linear'='purple'))+
  scale_y_continuous(name='Predicted value, ppb', breaks = seq(0,10,2))+
  scale_x_continuous(name='Reference value, ppb', breaks = seq(0,10,2))+
  theme(text = element_text(size = 15,face="bold", colour ='black'),
        axis.text.x  = element_text(size=13.5,colour = 'black'),
        axis.text.y = element_text(size=13.5,colour = 'black'),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.major = element_line(colour = "grey", linetype= 'dashed'),
        panel.grid.minor = element_line(colour = "grey", linetype= 'dashed'))



##################################################################################
########### LLDF, optical + ec (2)
x1  <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='D1:Z22'))
x2 <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='AP1:EK22'))
x <- cbind(x1,x2)
rownames(x) <- seq(0,20, 1)

y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                           range ='EL1:EL22'))
rownames(y) <- seq(0,20, 1)

x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '15' , '12')),]
y <-  y[!(row.names(y) %in% c('0','2', '5', '11', '15' , '12')),]


m_1 <- pls(x[,-c(1:8)],y, 10, cv = 1, scale = T, center = T)
#plot(m_1)
show(summary(m_1))
plotRMSE(m_1)
plotYCumVariance(m_1)
#look for the outliers
plotPredictions(m_1, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m_1, show.labels = T, type = 'b')

### KRLS (2)
y <- data.frame(y)
colnames(y) <- 'y'
data <- data.frame(cbind(y,x))

data <- data[,]
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
#### mid-level, PLS (3)

ec <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                                  range ='AP1:EK22'))
rownames(ec) <- seq(0,20, 1)

fl <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='D1:Z22'))
rownames(fl) <- seq(0,20, 1)
ms <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='AA1:AO22'))
rownames(ms) <- seq(0,20, 1)
wpi <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='EL1:EL22'))
rownames(wpi) <- seq(0,20, 1)
ec <-  ec[!(row.names(ec) %in% c('0','2', '5', '11', '15', '12')),]
fl <-  fl[!(row.names(fl) %in% c('0','2', '5', '11', '15', '12')),]
ms <-  ms[!(row.names(ms) %in% c('0','2', '5', '11', '15', '12')),]
wpi <-  wpi[!(row.names(wpi) %in% c('0','2', '5', '11', '15', '12')),]

ecpca <-  pca(ec, ncomp = 2, scale = T, center = T)
plotCumVariance(ecpca, show.labels = T, lab.cex = 0.8)

flpca <-  pca(fl, ncomp = 1, scale = T, center = T)
plotCumVariance(flpca, show.labels = T, lab.cex = 0.8)

mspca <-  pca(ms, ncomp = 1, scale = T, center = T)
plotCumVariance(mspca, show.labels = T, lab.cex = 0.8)

ec_scores <- data.frame(ecpca$res$cal$scores)
fl_scores <- data.frame(flpca$res$cal$scores)
ms_scores <- data.frame(mspca$res$cal$scores)

all <- cbind(ec_scores, fl_scores, ms_scores)

m <-  pls(all[,], wpi, 5, cv = 1, scale = T, center = T)
plotRMSE(m)
show(summary(m))
plotPredictions(m, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T, type = 'h', lab.cex = 0.8)
res_cal <- m[["res"]][["cal"]][["y.pred"]][,3,1]
res_cv <- m[["res"]][["cv"]][["y.pred"]][,3,1]
res_ref <-  m[["res"]][["cv"]][["y.ref"]]
res <-data.frame(cbind(res_cv, res_cal, res_ref))



ggplot(res) +
  geom_point(aes(x = res_ref , y = res_cal), size = 3, col = 'red', alpha = 0.7) +
  geom_point(aes(x = res_ref , y = res_cv), size = 3, col = 'blue', alpha = 0.7) +
  geom_smooth((aes(x = res_ref , y = res_cal)), method = "lm", se=FALSE, col = 'red',  alpha = 0.7, size = 1) +
  geom_smooth((aes(x = res_ref , y = res_cv)), method = "lm", se=FALSE, col = 'blue',  alpha = 0.7, size = 1)+
  theme_bw() +
  scale_y_continuous(name='Reference value, ppb')+
  scale_x_continuous(name='Predicted value, ppb')+
  theme(text = element_text(size = 15,face="bold", colour ='black'),
        axis.text.x  = element_text(size=13.5,colour = 'black'),
        axis.text.y = element_text(size=13.5,colour = 'black'),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.major = element_line(colour = "grey", linetype= 'dashed'),
        panel.grid.minor = element_line(colour = "grey", linetype= 'dashed'))
abline()



###############################################################################
################################################################################
####  mid-level KRLS (3)

d <- cbind(wpi, ec_scores, fl_scores, ms_scores)

colnames(d) <- c('WPI','1','2', '3', '4')

#d <- d[, -c(3)]

KRLS <- krls(d[,-1], d[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]

test_reg_loo_model <- train(WPI ~ ., data = d, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda =c(KRLS[["lambda"]], 0.001, 0.01,0.1348348, 0.05,0.065, 0.5,0.1,1, 0.2, 2), .sigma = (length(d)-1)))

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



################################################################################
##### MID-level DF PLS (2)
ec <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='AP1:EK22'))
rownames(ec) <- seq(0,20, 1)

fl <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='D1:Z22'))
rownames(fl) <- seq(0,20, 1)

wpi <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                             range ='EL1:EL22'))
rownames(wpi) <- seq(0,20, 1)
ec <-  ec[!(row.names(ec) %in% c('0','2', '5', '11', '15', '12')),]
fl <-  fl[!(row.names(fl) %in% c('0','2', '5', '11', '15', '12')),]
wpi <-  wpi[!(row.names(wpi) %in% c('0','2', '5', '11', '15', '12')),]

ecpca <-  pca(ec, ncomp = 2, scale = T, center = T)
plotCumVariance(ecpca, show.labels = T, lab.cex = 0.8)

flpca <-  pca(fl, ncomp = 1, scale = T, center = T)
plotCumVariance(flpca, show.labels = T, lab.cex = 0.8)

ec_scores <- data.frame(ecpca$res$cal$scores)
fl_scores <- data.frame(flpca$res$cal$scores)

all <- cbind(ec_scores, fl_scores)

m <-  pls(all[,], wpi, 5, cv = 1, scale = T, center = T)
plotRMSE(m)
show(summary(m))
plotPredictions(m, show.labels = T, lab.cex = 0.8)
abline(a = 0, b = 1)
plotRegcoeffs(m, show.labels = T, type = 'h', lab.cex = 0.8)

##########################################
library(caret)
library(KRLS)

d <- cbind(wpi, ec_scores, fl_scores)

colnames(d) <- c('WPI','1','2', '3')

#d <- d[, -c(3)]

KRLS <- krls(d[,-1], d[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]

test_reg_loo_model <- train(WPI ~ ., data = d, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda =c(KRLS[["lambda"]], 0.001, 0.01, 0.05,0.065, 0.5,0.1,1, 0.2, 2), .sigma = (length(d)-1)))

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