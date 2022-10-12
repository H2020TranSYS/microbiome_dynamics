library(readxl)
library(ggplot2)

data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Results/Prediction/"
graph_pred_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/Prediction/"
graph_path =  "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/Pie_chart//"


setwd(data_path)


Diet <- read_excel("Final_in_the_paper_results_prediction_plusSUPPL.xlsx", sheet = "Diet")
colnames(Diet) = c("Input", "AUC", "sd")

Input_name = c("Microbial abund. m6 (CLR)", "Microbial abund. m9 (CLR)", "Diff. microbial abund. m9 - m6 (CLR)", 
               "ISN weights m6 (MAGMA)", "ISN weights m9 (MAGMA)", "Diff. ISN weights m9 - m6 (MAGMA)",
               "ISN weights m6 (SparCC)", "ISN weights m9 (SparCC)", "Diff. ISN weights m9 - m6 (SparCC)",
               "Microbial dynamics MNDA (MAGMA)")
Diet$Input = Input_name
Delivery <- read_excel("Final_in_the_paper_results_prediction_plusSUPPL.xlsx", sheet = "Mode of Delivery")
colnames(Delivery) = c("Input", "AUC", "sd")
Delivery$Input = Input_name




Diet$Group = c("Abundance","Abundance","Abundance", "Edge MAGMA","Edge MAGMA","Edge MAGMA", 
               "Edge SparCC","Edge SparCC","Edge SparCC", "MNDA") 
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

p


setwd(graph_pred_path)

png("Suppl_Graph_Diet_persistent_vs_NONPERS.png", width = 300, height = 300, units='mm', res = 300)
p
dev.off()


# Mode of delivery  -------------------------------------------------------

Delivery$Group = c("Abundance","Abundance","Abundance", "Edge MAGMA","Edge MAGMA","Edge MAGMA", 
                   "Edge SparCC","Edge SparCC","Edge SparCC", "MNDA") 
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


p


png("Suppl_Graph_Delivery.png", width = 300, height = 300, units='mm', res = 300)
p
dev.off()
