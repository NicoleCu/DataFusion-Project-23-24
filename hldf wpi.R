library(readxl)
library(mdatools)
library(Metrics)

# getting predictions for all points and all datasets (with variable selection)

y_ref <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                               range ='EL1:EL22'))
rownames(y_ref) <- seq(0,20, 1)

MS <-  data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                             range ='AA1:AO22'))
rownames(MS) <- seq(0,20, 1)

x1  <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                             range ='D1:Z22'))
x2 <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'DF with new MS',
                            range ='AP1:EK22'))
x <- cbind(x1,x2)
rownames(x) <- seq(0,20, 1)


MS <-  MS[!(row.names(MS) %in% c('0','2', '5', '11', '12' , '15')),]
x <-  x[!(row.names(x) %in% c('0','2', '5', '11', '12' , '15')),]
y_ref <-  y_ref[!(row.names(y_ref) %in% c('0','2', '5', '11', '15' , '12')),]


m_ms <- pls(MS[,-c(11,1,8,9,12, 2,3,6, 10)], y_ref, 10, cv = 1, scale = T, center = T)   
plotRMSE(m_ms)
show(summary(m_ms))
plotRegcoeffs(m_ms, show.labels = T)

m_x <- pls(x[,-c(1:9)], y_ref, 10, cv = 1, scale = T, center = T)
plotRMSE(m_x)
show(summary(m_x))
plotYCumVariance(m_x)
plotRegcoeffs(m_x, show.labels = T)


y_mscal <-  data.frame(m_ms$res$cal$y.pred[,1,])
y_xcal <-  data.frame(m_x$res$cal$y.pred[,1,])


cd <- data.frame(cbind(y_ref, y_mscal, y_xcal))
colnames(cd) <- c('y_ref', 'y_ms','y_x')
y_ref1 <- y_ref
m <- lm(data = cd, y_ref ~ y_ms+ y_x)
summary(m)
y_ref <-sapply(y_ref, as.numeric)
#plot(y_ref, yp)


################################################################################
##### HLDF calibration and validation

y_msp <-  data.frame(m_ms$res$cv$y.pred[,1,])
y_xp <- data.frame( m_x$res$cv$y.pred[,1,])
daf_p <- data.frame(cbind(y_ref, y_msp, y_xp))
colnames(daf_p) <- c('y_ref', 'y_ms','y_x')

y_p <- predict(m, daf_p)

rmse <- rmse(y_ref, y_p)
plot(y_ref, y_p)
abline(a = 0, b =1, col = 'green')
text(y_ref, y_p,
     labels = (names(y_p)),
     cex = 0.6, pos = 4, col = "gray")
summary(m)

ggplot(res) +
  geom_point(aes(x = res_ref , y = res_cal), size = 3, col = 'red', alpha = 0.7) +
  geom_point(aes(x = res_ref , y = res_pred), size = 3, col = 'blue', alpha = 0.7) +
  geom_smooth((aes(x = res_ref , y = res_cal)), method = "lm", se=FALSE, col = 'red',  alpha = 0.7, size = 1) +
  geom_smooth((aes(x = res_ref , y = res_pred)), method = "lm", se=FALSE, col = 'blue',  alpha = 0.7, size = 1)+
  theme_bw() +
  scale_color_manual(name='Regression Model',
                     breaks=c('Linear', 'Quadratic', 'Cubic'),
                     values=c('Cubic'='pink', 'Quadratic'='blue', 'Linear'='purple'))+
  scale_y_continuous(name='????????????? ????????, ppb')+
  scale_x_continuous(name='???????? ????????, ppb')+
  theme(text = element_text(size = 15,face="bold", colour ='black'),
        axis.text.x  = element_text(size=13.5,colour = 'black'),
        axis.text.y = element_text(size=13.5,colour = 'black'),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.major = element_line(colour = "grey", linetype= 'dashed'),
        panel.grid.minor = element_line(colour = "grey", linetype= 'dashed'))







