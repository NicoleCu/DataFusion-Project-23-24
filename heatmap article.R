################################################### heatamap by sensor responses
#new dataset

library(pheatmap)
library(mdatools)
library(readxl)

#new data set
x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                      sheet = 'Heatmap', range ='C1:Q17')))
x <- data.frame(scale(x, center = TRUE, scale = TRUE))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                              sheet = 'Heatmap', range ='A1:A17'))
row.names(x) <- names
  
df_ann <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                 sheet = 'Heatmap', range ='B1:B17')))
row.names(df_ann) <- names

library('RColorBrewer')
ms <-  pheatmap(x,  
        treeheight_row = 0, 
        treeheight_col = 0,
        cluster_rows = F,
        cluster_cols = F,
        legend = T,
        fontsize = 14,
        angle_col = 0,
        border_color = 'white',
        annotation_row = df_ann)

# old dataset
x_old <- data.frame(read_excel("~/R/Chemometrics/ModelWaters-Daphnia.xlsx",
                           range ='B2:L22'))
x_old <- data.frame(scale(x_old, center = TRUE, scale = TRUE))
row.names(x_old) <- seq(1, 20, 1)

pheatmap(x_old,  
         treeheight_row = 0, 
         treeheight_col = 0,
         cluster_rows = F,
         cluster_cols = F,
         show_rownames = T, fontsize_row = 20, 
         show_colnames = T, fontsize_col = 20,
         border_color = "white",
         col = rev(brewer.pal(10, 'RdYlBu')),
         border_col = F,
         annotation_row = df_ann)

################################################################################
#################################### spectrum (fluorescence)

x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C20:Y36')))
x <- data.frame(scale(x, center = TRUE, scale = TRUE))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A20:A36'))
row.names(x) <- names

df_ann <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                 sheet = 'Heatmap', range ='B20:B36')))
row.names(df_ann) <- names

library('RColorBrewer')
fl <- pheatmap(x,  
         treeheight_row = 0, 
         treeheight_col = 0,
         cluster_rows = F,
         cluster_cols = F,
         legend = T,
         show_rownames = T, 
         show_colnames = F,
         fontsize = 14,
         angle_col = 0,
         border_color = 'white',
         annotation_row = df_ann)

#library(writexl)
#write_xlsx(x_fl, "~/R/Chemometrics/Fluo.xlsx")
#x_fl <- data.frame((read_excel("~/R/Chemometrics/Fluo.xlsx", range ='A1:C484')))

plot(as.numeric(seq(1,23,1)), as.numeric(x_fl[1,]), type ='l', ylim = c(2,-3))
lines(as.numeric(seq(1,23,1)), x_fl[2,], type ='l', col = 'red')
lines(as.numeric(seq(1,23,1)), x_fl[3,], type ='l', col = 'purple')
lines(as.numeric(seq(1,23,1)), x_fl[4,], type ='l', col = 'green')
lines(as.numeric(seq(1,23,1)), x_fl[5,], type ='l', col = 'orange')



################################################################################
#################################### spectrum (voltammetry)

x <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                            sheet = 'Heatmap', range ='C39:CX54')))
x <- data.frame(scale(x, center = TRUE, scale = TRUE))
names <- t(read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                      sheet = 'Heatmap', range ='A39:A54'))
row.names(x) <- names

df_ann <- data.frame((read_excel("~/R/Chemometrics/China MW December 22.xlsx", 
                                 sheet = 'Heatmap', range ='B39:B54')))
row.names(df_ann) <- names

library('RColorBrewer')
ec <- pheatmap(x,  
         treeheight_row = 0, 
         treeheight_col = 0,
         cluster_rows = F,
         cluster_cols = F,
         legend = T,
         legend_breaks = c(2,1,0,-1,-2),
         show_rownames = T, 
         show_colnames = F,
         fontsize = 14,
         angle_col = 0,
         border_color = 'white',
         annotation_row = df_ann)

#m <- rbind(seq(1,100,1), seq(1,100,1), seq(1,100,1), seq(1,100,1), seq(1,100,1))
#m <- rbind(m, m, m)
par(mar=c(1,1,1,1))

plot(x_el[1,], seq(1,100,1))

lines(x_el[2,], seq(1,100,1))


#matplot(m,x_el[1:15,],xlab = "Wavelength [nm]", ylab = "Absorbance", type = "p", pch = 1)
library(ggpubr)

ggarrange(ggarrange(ms$gtable,fl$gtable, ncol = 2, labels = c("A", "B")), ec$gtable,
          nrow = 2) 