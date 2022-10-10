setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Dynamics_in_Microbiome/Results/Prediction/")
library(readxl)

Diet <- read_excel("Final_in_the_paper_results_prediction.xlsx", sheet = "Diet")
colnames(Diet) = c("Input", "AUC", "sd")

Input_name = c("Microbial abund. m6 (CLR)", "Microbial abund. m9 (CLR)", "Diff. microbial abund. m9 - m6 (CLR)", 
               "ISN weights m6 (MAGMA)", "ISN weights m9 (MAGMA)", "Diff. ISN weights m9 - m6 (MAGMA)", 
               "Microbial dynamics MNDA (MAGMA)")
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
p<-ggplot(Diet, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=Group), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  theme_classic(base_size = 17)+
  theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
        axis.text=element_text(face="bold"),
              axis.title=element_text(face="bold"), 
        ,legend.text=element_text(face="bold"), legend.title=element_blank(), 
        plot.title = element_text(size=16, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Diet (Persistent vs Non persistent)")#+
  # theme_bw(base_size = 17, base_line_size = 1.1)+
  # theme_classic()
  # theme_classic(base_size = 17, base_line_size = 1.1)+
  # theme(axis.text=element_text(face="bold"),
  #       axis.title=element_text(face="bold"), title = element_text(face = "bold",),legend.key.size = unit(1, 'cm'),text = element_text(size = 35),legend.position="top")
p


png("Graph_Diet_persistent_vs_NONPERS_1sd.png", width = 300, height = 300, units='mm', res = 300)
p
dev.off()


p<-ggplot(Diet, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=Group), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-2*sd, ymax=AUC+2*sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  theme_classic(base_size = 17)+
  theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
        axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), 
        ,legend.text=element_text(face="bold"), legend.title=element_blank(), 
        plot.title = element_text(size=16, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Diet (Persistent vs Non persistent)")#+
# theme_bw(base_size = 17, base_line_size = 1.1)+
# theme_classic()
# theme_classic(base_size = 17, base_line_size = 1.1)+
# theme(axis.text=element_text(face="bold"),
#       axis.title=element_text(face="bold"), title = element_text(face = "bold",),legend.key.size = unit(1, 'cm'),text = element_text(size = 35),legend.position="top")
p


png("Graph_Diet_persistent_vs_NONPERS_2sd.png", width = 300, height = 300, units='mm', res = 300)
p
dev.off()



# Mode of delivery  -------------------------------------------------------

Delivery$Group = c("Abundance","Abundance","Abundance", "Edge", "Edge","Edge", "MNDA") 
Delivery$Group = as.factor(Delivery$Group)
level_order = c(Delivery$Input)
p<-ggplot(Delivery, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=Group), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-sd, ymax=AUC+sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  theme_classic(base_size = 17)+
  theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
        axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), 
        ,legend.text=element_text(face="bold"), legend.title=element_blank(), 
        plot.title = element_text(size=16, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Delivery type")#+
# theme_bw(base_size = 17, base_line_size = 1.1)+
# theme_classic()
# theme_classic(base_size = 17, base_line_size = 1.1)+
# theme(axis.text=element_text(face="bold"),
#       axis.title=element_text(face="bold"), title = element_text(face = "bold",),legend.key.size = unit(1, 'cm'),text = element_text(size = 35),legend.position="top")
p


png("Graph_Delivery_1sd.png", width = 300, height = 300, units='mm', res = 300)
p
dev.off()


p<-ggplot(Delivery, aes(x=factor(Input, level =level_order), y=AUC)) + 
  geom_point(aes(color=Group), shape = 18, size = 10)+
  geom_errorbar(aes(ymin=AUC-2*sd, ymax=AUC+2*sd), width=.3,size = 0.8,
                position=position_dodge(0.05))+
  theme_classic(base_size = 17)+
  theme(axis.text.x = element_text(angle = 60, vjust =1, hjust=1), plot.margin = unit(c(1,0,0,1), "cm"),
        axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), 
        ,legend.text=element_text(face="bold"), legend.title=element_blank(), 
        plot.title = element_text(size=16, face="bold",hjust = 0.5))+
  labs(y= "Mean AUC",x = "", fill = "")+
  ggtitle("Delivery type")#+
# theme_bw(base_size = 17, base_line_size = 1.1)+
# theme_classic()
# theme_classic(base_size = 17, base_line_size = 1.1)+
# theme(axis.text=element_text(face="bold"),
#       axis.title=element_text(face="bold"), title = element_text(face = "bold",),legend.key.size = unit(1, 'cm'),text = element_text(size = 35),legend.position="top")
p


png("Graph_Delivery_2sd.png", width = 300, height = 300, units='mm', res = 300)
p
dev.off()



