library(KRLS)
library(caret)
library(readxl)
library(mdatools)

####### WPI, KRLS


x <-  data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                 sheet = 'av 2-3', range ='B1:P22'))
x <- x[-1,]
rownames(x) <- seq(1,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'DF with new MS', range ='EL1:EL21'))
rownames(y) <- seq(1,20, 1)

data <- data.frame(cbind(y, x))
colnames(y) <- 'y'

data <- data.frame(cbind(y,x))
data <-  data[!(row.names(data) %in% c('2', '1')),]
data <- data[,-c(4,7,15,8)]

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