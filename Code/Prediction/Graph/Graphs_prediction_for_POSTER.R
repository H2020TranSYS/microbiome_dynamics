setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Dynamics_in_Microbiome/Results/Prediction/")
library(readxl)
library(RColorBrewer)

Diet <- read_excel("Final_in_the_paper_results_prediction.xlsx", sheet = "Diet")
colnames(Diet) = c("Input", "AUC", "sd")

Input_name = c("Microbial abund. M6 (CLR)", "Microbial abund. M9 (CLR)", "Diff. microbial abund. M9 - M6 (CLR)", 
               "ISN weights M6 (MAGMA)", "ISN weights M9 (MAGMA)", "Diff. ISN weights M9 - M6 (MAGMA)", 
               "Microbial dynamics MNDA (MAGMA)")
Input_name = c("A","B", "C", "D", "E", "F", "G")
Diet$Input = Input_name
Delivery <- read_excel("Final_in_the_paper_results_prediction.xlsx", sheet = "Mode of Delivery")
colnames(Delivery) = c("Input", "AUC", "sd")
Delivery$Input = Input_name

library(ggplot2)

# creating a data frame df
df<-data.frame(Mean=c(0.24,0.25,0.37,0.643,0.54),
               sd=c(0.00362,0.281,0.3068,0.2432,0.322),
               Quality=as.factor(c("good","bad","good",
                                   "very good","very good")), 
               Category=c("A","B","C","D","E"),
               Insert= c(0.0, 0.1, 0.3, 0.5, 1.0))

Diet$Group = c("Abundance","Abundance","Abundance", "Edge", "Edge","Edge", "MNDA") 
Diet$Group = as.factor(Diet$Group)
level_order = c(Diet$Input)


col_brew = rep(brewer.pal(n = 3, name = "Dark2"),each = 3)[1:7]

p<-ggplot(Diet, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=level_order), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  scale_color_manual(name="Cylinders",
                     labels=c("A: Microbial abund. M6 (CLR)", "B: Microbial abund. M9 (CLR)", "C: Diff. microbial abund. M9 - M6 (CLR)", 
               "D: ISN weights M6 (MAGMA)", "E: ISN weights M9 (MAGMA)", "F: Diff. ISN weights M9 - M6 (MAGMA)", 
               "G: Microbial dynamics MNDA (MAGMA)"),
                     values=col_brew)+
  theme_classic(base_size = 25)+
  ylim(0.35,0.85)+
  theme( #axis.text.x = element_text(angle = 25, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
    axis.text=element_text(size = 25,face="bold"),
    axis.title=element_text(face="bold"), 
    legend.text=element_text(face="bold"), legend.title=element_blank(), 
    plot.title = element_text(size=25, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Diet (Persistent vs Non persistent)")#+
# theme_bw(base_size = 17, base_line_size = 1.1)+
# theme_classic()
# theme_classic(base_size = 17, base_line_size = 1.1)+
# theme(axis.text=element_text(face="bold"),
#       axis.title=element_text(face="bold"), title = element_text(face = "bold",),legend.key.size = unit(1, 'cm'),text = element_text(size = 35),legend.position="top")
p




png("Graph_Diet_persistent_vs_NONPERS_1sd_FOR_POSTER_FULL_b.png", width = 600, height = 300, units='mm', res = 300)
p
dev.off()

pdf("Graph_Diet_persistent_vs_NONPERS_1sd_FOR_POSTER_FULL.pdf", width = 15, height = 10)#, units='mm', res = 300)
p
dev.off()


# SAME but only three saving ----------------------------------------------

p<-ggplot(Diet, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=Group), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  scale_color_manual(name="Cylinders",
                     # labels=c("A: Microbial abund. M6 (CLR)", "B: Microbial abund. M9 (CLR)", "C: Diff. microbial abund. M9 - M6 (CLR)", 
                     #          "D: ISN weights M6 (MAGMA)", "E: ISN weights M9 (MAGMA)", "F: Diff. ISN weights M9 - M6 (MAGMA)", 
                     #          "G: Microbial dynamics MNDA (MAGMA)"),
                     values=brewer.pal(n = 3, name = "Dark2"))+
  theme_classic(base_size = 17)+
  ylim(0.35,0.85)+
  theme( #axis.text.x = element_text(angle = 25, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
    axis.text=element_text(size = 25,face="bold"),
    axis.title=element_text(face="bold"), 
    legend.text=element_text(face="bold"), legend.title=element_blank(), 
    plot.title = element_text(size=25, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Diet (Persistent vs Non persistent)")#+


png("Graph_Diet_persistent_vs_NONPERS_1sd_FOR_POSTER_abu_edg_MNDA.png", width = 600, height = 300, units='mm', res = 300)
p
dev.off()

pdf("Graph_Diet_persistent_vs_NONPERS_1sd_FOR_POSTER_abu_edg_MNDA.pdf", width = 15, height = 10)#, units='mm', res = 300)
p
dev.off()




# Mode of delivery  -------------------------------------------------------

Delivery$Group = c("Abundance","Abundance","Abundance", "Edge", "Edge","Edge", "MNDA") 
Delivery$Group = as.factor(Delivery$Group)
level_order = c(Delivery$Input)
p<-ggplot(Delivery, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=level_order), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  scale_color_manual(name="Cylinders",
                     labels=c("A: Microbial abund. M6 (CLR)", "B: Microbial abund. M9 (CLR)", "C: Diff. microbial abund. M9 - M6 (CLR)", 
                              "D: ISN weights M6 (MAGMA)", "E: ISN weights M9 (MAGMA)", "F: Diff. ISN weights M9 - M6 (MAGMA)", 
                              "G: Microbial dynamics MNDA (MAGMA)"),
                     values=col_brew)+
  theme_classic(base_size = 25)+
  ylim(0.35,0.85)+
  theme(# axis.text.x = element_text(angle = 60, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
        axis.text=element_text(size = 25, face="bold"),
        axis.title=element_text(face="bold"), 
        legend.text=element_text(face="bold"), legend.title=element_blank(), 
        plot.title = element_text(size=25, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Delivery type")#+
# theme_bw(base_size = 17, base_line_size = 1.1)+
# theme_classic()
# theme_classic(base_size = 17, base_line_size = 1.1)+
# theme(axis.text=element_text(face="bold"),
#       axis.title=element_text(face="bold"), title = element_text(face = "bold",),legend.key.size = unit(1, 'cm'),text = element_text(size = 35),legend.position="top")
p


png("Graph_Delivery_persistent_vs_NONPERS_1sd_FOR_POSTER_FULL_b.png", width = 600, height = 300, units='mm', res = 300)
p
dev.off()

pdf("Graph_Delivery_persistent_vs_NONPERS_1sd_FOR_POSTER_FULL.pdf", width = 15, height = 10)#, units='mm', res = 300)
p
dev.off()


# Same but only 3 saving --------------------------------------------------



p<-ggplot(Delivery, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=Group), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  scale_color_manual(name="Cylinders",
                     # labels=c("A: Microbial abund. M6 (CLR)", "B: Microbial abund. M9 (CLR)", "C: Diff. microbial abund. M9 - M6 (CLR)", 
                     #          "D: ISN weights M6 (MAGMA)", "E: ISN weights M9 (MAGMA)", "F: Diff. ISN weights M9 - M6 (MAGMA)", 
                     #          "G: Microbial dynamics MNDA (MAGMA)"),
                     values=brewer.pal(n = 3, name = "Dark2"))+
  theme_classic(base_size = 17)+
  ylim(0.35,0.85)+
  theme(# axis.text.x = element_text(angle = 60, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
    axis.text=element_text(size = 25, face="bold"),
    axis.title=element_text(face="bold"), 
    ,legend.text=element_text(face="bold"), legend.title=element_blank(), 
    plot.title = element_text(size=25, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Delivery type")#+



png("Graph_Delivery_persistent_vs_NONPERS_1sd_FOR_POSTER_abu_edg_MNDA.png", width = 600, height = 300, units='mm', res = 300)
p
dev.off()

pdf("Graph_Delivery_persistent_vs_NONPERS_1sd_FOR_POSTER_abu_edg_MNDA.pdf", width = 15, height = 10)#, units='mm', res = 300)
p
dev.off()

