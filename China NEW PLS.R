# Model waters -- new data, December 22, multisensory system
library(readxl)
library(mdatools)
library(Metrics)

#PLS, cd
# multisensor system
x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'av 2-3', range ='B1:P22'))    
rownames(x) <- seq(0,20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'av 2-3',
                           range ='Q1:Q22'))
rownames(y) <- seq(0,20, 1)
x <-  x[!(row.names(x) %in% c('0', '16', '11', '17', '5')),]
y <-  y[!(row.names(y) %in% c('0','16', '11', '17', '5')),]

m_pb <- pls(x[,-c(1, 8,9,11,12)], y[], 10, cv = 1, scale = T) 
plot(m_pb)
show(summary(m_pb))
plotYCumVariance(m_pb, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_pb, type  = 'h', show.labels = TRUE)

# looking for the outliers
plotRegcoeffs(m_pb, show.labels = T)
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)



######################################################PLS, pb
# multisensor system
x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'av 2-3',
                           range ='B1:P22'))
rownames(x) <- seq(0, 20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx",  sheet = 'av 2-3',
                           range ='R1:R22'))   
rownames(y) <- seq(0, 20, 1)
x <-  x[!(row.names(x) %in% c('0','2')),]
y <-  y[!(row.names(y) %in% c('0','2')),]

m_pb <- pls(x[], y[], 10, cv = 1, scale = T)
plot(m_pb)
show(summary(m_pb))
plotYCumVariance(m_pb, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_pb, type  = 'h', show.labels = TRUE)

# look for the outliers
plotRegcoeffs(m_pb, show.labels = T)
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)


######################################################PLS, WPI
# multisensor system
x <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                           sheet = 'av 2-3', range ='B1:P22'))
rownames(x) <- seq(0, 20, 1)
y <- data.frame(read_excel("~/R/Chemometrics/China MW December 22.xlsx", sheet = 'av 2-3',
                           range ='U1:U22'))
#y <- y[,-2]
rownames(y) <- seq(0, 20, 1)
x <-  x[!(row.names(x) %in% c('0', '2')),]
y <-  y[!(row.names(y) %in% c('0', '2')),]

m_pb <- pls(x[,-c(14)], y[], 10, cv = 1, scale = T)
plot(m_pb)
show(summary(m_pb))
plotYCumVariance(m_pb, type  = 'h', show.labels = TRUE)
plotXCumVariance(m_pb, type  = 'h', show.labels = TRUE)

# look for the outliers
plotRegcoeffs(m_pb, show.labels = T)
plotPredictions(m_pb, show.labels = T)
abline(a = 0, b = 1)



