#LLDF, China, new data for MS
library(readxl)
library(mdatools)
library(Metrics)


# PLS, Cd
x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='C25:AN45'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='B25:B45'))
rownames(y) <- seq(1,20, 1)
x <-  x[!(row.names(x) %in% c('16', '11', '5', '12')),]
y <-  y[!(row.names(y) %in% c('16', '11', '5', '12')),]

m_cd <- pls(x[,-c(1,8,9)], y, 10, cv = 1, scale = T) 
plot(m_cd)
show(summary(m_cd))

# looking for the outliers
plotRegcoeffs(m_cd, show.labels = T, type = "h")
plotPredictions(m_cd, show.labels = T)
abline(a = 0, b = 1)



# KRLS, Cd

library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='C25:AN45'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='B25:B45'))
rownames(y) <- seq(1,20, 1)
colnames(y) <- 'y'

data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('5', '14', '15')),]
#data <- data[,-c(5,11)]

KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(data)-1)))
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




# PLS, lgaCd

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='C25:AN45'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='AO25:AO45'))
rownames(y) <- seq(1,20, 1)
x <-  x[!(row.names(x) %in% c('11', '5', '12')),]
y <-  y[!(row.names(y) %in% c('11', '5', '12')),]

m_cd <- pls(x[,], y, 10, cv = 1, scale = T) 
plot(m_cd)
show(summary(m_cd))

# looking for the outliers
plotRegcoeffs(m_cd, show.labels = T, type = "h")
plotPredictions(m_cd, show.labels = T)
abline(a = 0, b = 1)



#KRLS, lgaCd

library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='C25:AN45'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='AO25:AO45'))
rownames(y) <- seq(1,20, 1)
colnames(y) <- 'y'

data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('1', '5', '14')),]
#data <- data[,-c(5,11)]

KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(data)-1)))
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













# PLS, Pb
x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='AA1:EK21'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='C1:C21'))
rownames(y) <- seq(1,20, 1)
x <-  x[!(row.names(x) %in% c('2')),]
y <-  y[!(row.names(y) %in% c('2')),]

m_pb <- pls(x[,], y, 10, cv = 1, scale = T) 
plot(m_pb)
show(summary(m_pb))

# looking for the outliers
plotRegcoeffs(m_pb, show.labels = T, type = "h")
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)



#KRLS, Pb

library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='AA1:EK21'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='C1:C21'))
rownames(y) <- seq(1,20, 1)
colnames(y) <- 'y'

data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('2')),]
data <- data[,-c(7,10,12)]

KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(data)-1)))
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





# PLS, WPI (3)
x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='D1:EK21'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='EL1:EL21'))
rownames(y) <- seq(1,20, 1)
x <-  x[!(row.names(x) %in% c('2', '17')),]
y <-  y[!(row.names(y) %in% c('2', '17')),]

m_wpi <- pls(x[,], y, 10, cv = 1, scale = T) 
plot(m_wpi)
show(summary(m_wpi))

# looking for the outliers
plotRegcoeffs(m_wpi, show.labels = T, type = "h")
plotPredictions(m_wpi, show.labels = T)
abline(a = 0, b = 1)




# KRLS, WPI (3)
library(KRLS)
library(caret)

x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS', range ='D1:EK21'))    
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='EL1:EL21'))
rownames(y) <- seq(1,20, 1)
colnames(y) <- 'y'

data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('1', '2')),]
#data <- data[,-c(7,10,12)]

KRLS <- krls(data[,-1], data[,1] , vcov = T, derivative = T)
KRLS[["lambda"]]
test_reg_loo_model <- train(y ~ ., data = data, method = "krlsRadial", trControl = trainControl(method = "LOOCV") ,
                            preProc = c("center", "scale"), print.level = 0,
                            tuneGrid = data.frame(.lambda = c(KRLS[["lambda"]]), .sigma = (length(data)-1)))
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

