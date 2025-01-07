#Spectrums
library(readxl)
library(writexl)


# Fluorescent

x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C20:Y36')))
x <- data.frame(scale(x, center = TRUE, scale = TRUE))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A20:A36'))
row.names(x) <- names

x <- data.frame(t(x))

library(writexl)
write_xlsx(x, "~/R/Chemometrics/spectrums.xlsx", 
           col_names=TRUE)

library(ggplot2)

spectra <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                  sheet = 'Spectrums', range ='A1:C369')))
spectra$sample_name = factor(spectra$name, levels=spectra$name, labels=spectra$name) 
spectrum <- 
  ggplot(spectra) +
  aes(x = Wavelength, y = Intensity_scaled, color = sample_name)+
  geom_line(linewidth = 0.8)+
  theme_bw()+
  facet_wrap(.~sample_name)
spectrum

###### without scaling


x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C20:Y36')))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A20:A36'))
row.names(x) <- names

x <- data.frame(t(x))

library(writexl)
write_xlsx(x, "~/R/Chemometrics/spectrums withot scaling.xlsx", 
           col_names=TRUE)




library(ggplot2)
library(ggrepel)
library(dplyr)

spectra <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                  sheet = 'Spectrums', range ='K1:M369')))
spectra$sample_name = factor(spectra$name, levels=spectra$name, labels=spectra$name) 
spectrum <- 
  ggplot(spectra) +
  aes(x = Wavelength, y = Intensity, color = sample_name, geom_text(~spectra$name))+
  geom_line(linewidth = 0.8)+
  theme_bw()+
  scale_y_continuous(name='Intensity')+
  scale_x_continuous(name='Wavelength, nm')+
  theme(text = element_text(size = 15,face="bold", colour ='black'),
           axis.text.x  = element_text(size=13.5,colour = 'black'),
           axis.text.y = element_text(size=13.5,colour = 'black'),
           panel.border = element_rect(colour = "black", size=1),
           panel.grid.major = element_line(colour = "grey", linetype= 'dashed'),
           panel.grid.minor = element_line(colour = "grey", linetype= 'dashed'))
spectrum





################### Voltammery

x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C39:CX54')))
x <- data.frame(scale(x, center = TRUE, scale = TRUE))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A39:A54'))
row.names(x) <- names
x <- data.frame(t(x))

library(writexl)
write_xlsx(x, "~/R/Chemometrics/voltammetry.xlsx", 
           col_names=TRUE)


library(ggplot2)
spectra <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                       sheet = 'Spectrums', range ='F1:H1501')))

spectra$new = factor(spectra$name, levels=spectra$name, labels=spectra$name) 

spectrum <- 
  ggplot(spectra) +
  aes(x = Potential, y = Current_scaled, color = new)+
  geom_line(linewidth = 0.8)+
  theme_bw()+
  facet_wrap(.~new)
spectrum


##### withot scaling

x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C39:CX54')))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A39:A54'))
row.names(x) <- names
x <- data.frame(t(x))

library(writexl)
write_xlsx(x, "~/R/Chemometrics/voltammetry withot scaling.xlsx", 
           col_names=TRUE)


library(ggplot2)
spectra <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                  sheet = 'Spectrums', range ='P1:R1501')))

spectra$sample_name = factor(spectra$name, levels=spectra$name, labels=spectra$name) 

spectrum <- 
  ggplot(spectra) +
  aes(x = Potential, y = Current, color = sample_name)+
  geom_line(linewidth = 0.8)+
  theme_bw()+
  scale_y_continuous(name='Current, mA')+
  scale_x_continuous(name='Potential, mV')+
  theme(text = element_text(size = 15,face="bold", colour ='black'),
        axis.text.x  = element_text(size=13.5,colour = 'black'),
        axis.text.y = element_text(size=13.5,colour = 'black'),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.major = element_line(colour = "grey", linetype= 'dashed'),
        panel.grid.minor = element_line(colour = "grey", linetype= 'dashed'))
spectrum




##### Multisensor system

x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C1:Q17')))
x <- data.frame(scale(x, center = TRUE, scale = TRUE))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A1:A17'))
row.names(x) <- names
x <- data.frame(t(x))


library(writexl)
write_xlsx(x, "~/R/Chemometrics/ms with scaling.xlsx", 
           col_names=TRUE)

spectra <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                  sheet = 'Spectrums', range ='U1:W241')))

spectra$sample_name = factor(spectra$name, levels=spectra$name, labels=spectra$name) 
spectra$sensor_name = factor(spectra$Sensor, levels=c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15'),
                             labels =c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15')) 

spectrum <- 
  ggplot(spectra) +
  aes(x = sensor_name, y = EMF, fill = sample_name)+
  geom_col()+
  theme_bw()+
  facet_wrap(~sample_name)
spectrum




#new data set
x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C1:Q17')))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A1:A17'))
row.names(x) <- names
x <- data.frame(t(x))

library(writexl)
write_xlsx(x, "~/R/Chemometrics/ms without scaling.xlsx", 
           col_names=TRUE)

spectra <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                  sheet = 'Spectrums', range ='Z1:AB241')))

spectra$sample_name = factor(spectra$name, levels=spectra$name, labels=spectra$name) 
spectra$sensor_name = factor(spectra$Sensor, levels=c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15'),
                             labels =c('S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15')) 


spectrum <- 
  ggplot(spectra) +
  aes(x = sensor_name, y = EMF, fill = sample_name)+
  geom_bar(position = 'dodge', stat = 'identity', alpha = 0.8, color = 'black', width = 0.96)+
  theme_bw()+
  labs(x = 'Sensor', y = 'EMF, mV')+
  theme(text = element_text(size = 15,face="bold", colour ='black'),
        axis.text.x  = element_text(size=13.5,colour = 'black'),
        axis.text.y = element_text(size=13.5,colour = 'black'),
        panel.border = element_rect(colour = "black", size=1),
        panel.grid.major = element_line(colour = "grey", linetype= 'dashed'),
        panel.grid.minor = element_line(colour = "grey", linetype= 'dashed'))
spectrum
