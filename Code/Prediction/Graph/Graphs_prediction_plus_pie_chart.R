library(RColorBrewer)
library(readxl)
library(ggplot2)
library(dplyr)

data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Results/Prediction/"
graph_pred_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/Prediction/"
graph_path =  "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/Pie_chart//"



setwd(data_path)

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




# Diet --------------------------------------------------------------------

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

p


setwd(graph_pred_path)

png("Graph_Diet_persistent_vs_NONPERS_FOR_POSTER_FULL.png", width = 600, height = 300, units='mm', res = 300)
p
dev.off()

pdf("Graph_Diet_persistent_vs_NONPERS_FOR_POSTER_FULL.pdf", width = 15, height = 10)#, units='mm', res = 300)
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

p


png("Graph_Delivery_persistent_vs_NONPERS_FOR_POSTER_FULL.png", width = 600, height = 300, units='mm', res = 300)
p
dev.off()

pdf("Graph_Delivery_persistent_vs_NONPERS_FOR_POSTER_FULL.pdf", width = 15, height = 10)#, units='mm', res = 300)
p
dev.off()


#  PIE CHART --------------------------------------------------------------


setwd(data_path)


Complete_annotated_list_95microbes <- read_excel("../Selected_95_microbes/Complete_annotated_list_95microbes.xlsx")


table(Complete_annotated_list_95microbes[1,])

t_C = t(Complete_annotated_list_95microbes)
colnames(t_C) = t_C[1,]
t_C = t_C[-1,]
t_C[,2] = gsub("p__","", t_C[,2])
tt = table(t_C[,2])


my_new_df = data.frame(group = names(tt), count = as.numeric(tt))
my_new_df$count = as.numeric(my_new_df$count)/sum(my_new_df$count)
rownames(my_new_df) = paste0("Phylum_", rownames(my_new_df))
right_colours = brewer.pal(n = 6, name = "Dark2")




# WITH OTHERS -------------------------------------------------------------

my_new_df$freq = as.numeric(my_new_df$count)/sum(my_new_df$count)
my_df_b = data.frame(causes =  my_new_df$group, share = 100*my_new_df$freq, freq = my_new_df$count)
mydf = my_df_b


mydf$Label  <- round(((mydf$freq/sum(mydf$freq))*100),0)
mydf <- mydf %>%
  mutate(end = 2 * pi * cumsum(Label)/sum(Label),
         start = lag(end, default = 0),
         middle = 0.5 * (start + end),
         hjust = ifelse(middle > pi, 1, 0),
         vjust = ifelse(middle < pi/2 | middle > 3 * pi/2, 0, 1))
# 

mydf$Label  <- paste0(mydf$Label, "%")

mydf[6,"hjust"] =mydf[6,"hjust"] - 0.85
mydf[6,"vjust"] =mydf[6,"vjust"]+ 0.05
mydf[5,"hjust"] =mydf[5,"hjust"] - 0.1


bp = ggplot(mydf) + 
  ggforce::geom_arc_bar(aes(x0 = 0, y0 = 0, r0 = 0, r = 1,
                            start = start, end = end, fill = causes)) +
  geom_text( aes(x = 0.8 * sin(middle), y = 1.05 * cos(middle), label = Label,
                 hjust = hjust, vjust = vjust), size = 10) +
  coord_fixed() +
  labs(x = NULL, y = NULL, fill = NULL, 
       title = "Pie chart microbes") +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_fill_manual(values =right_colours) +
  theme_classic(base_size = 25, base_line_size = 1.1)+
  # theme(axis.text=element_text(face="bold", size = 20),
  #       axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=25))
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0.5),#, color = "#666666"), 
        legend.text=element_text(size=25) )


setwd(graph_path)


png("Graph_Pie_CHART.png", width = 300, height = 300, units='mm', res = 300)
bp
dev.off()

pdf("Graph_Pie_CHART.pdf", width = 15, height = 10)#, units='mm', res = 300)
bp
dev.off()


