setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Dynamics_in_Microbiome/Results/Prediction/")


Complete_annotated_list_95microbes <- read_excel("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/Dynamics_in_Microbiome/Results/Selected_95_microbes/Complete_annotated_list_95microbes.xlsx")


table(Complete_annotated_list_95microbes[1,])

t_C = t(Complete_annotated_list_95microbes)
colnames(t_C) = t_C[1,]
t_C = t_C[-1,]



t_C[,2] = gsub("p__","", t_C[,2])
tt = table(t_C[,2])

methods(class = class(tt) )
attr(tt, "info")
my_new_df = data.frame(group = names(tt), count = as.numeric(tt))
my_new_df$count = as.numeric(my_new_df$count)/sum(my_new_df$count)
rownames(my_new_df) = paste0("Phylum_", rownames(my_new_df))
right_colours = brewer.pal(n = 6, name = "Dark2")



# BARPLOT -----------------------------------------------------------------



bp = ggplot(my_new_df, aes(x="", y=count, fill=group))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=right_colours, breaks=my_new_df$group)+
  labs(y= "Fractional abundances", x = "",fill='') +
  ggtitle(paste0("Barplot "))+
  theme_classic(base_size = 25, base_line_size = 1.1)+
  # geom_text(aes(label=paste0(sprintf("%1.f", count*100),"%")),
  #           position=position_stack(vjust=0.5), fontface = "bold") +
  theme(axis.text=element_text(face="bold"), #text = element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=25))
bp

geom_text(aes(label=paste0(sprintf("%1.f", count*100),"%")),
          position=position_stack(vjust=0.5), fontface = "bold")
bp2<- bp +   geom_text(aes(label=paste0(sprintf("%1.f",100*count),"%")),
                       position=position_stack(vjust=0.5), fontface = "bold")

bp2



pie <- ggplot(my_new_df, aes(x="", y=count, fill=group))+
  geom_bar(position="fill", stat="identity")+
  scale_fill_manual(values=right_colours, breaks=my_new_df$group)+
  labs(y= "Fractional abundances", x = "",fill='') +
  ggtitle(paste0("Pie Chart" ))+
  theme_minimal(base_size = 25, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold"),legend.text=element_text(size=25))+
  coord_polar("y", start=0)
# geom_text(aes(x=1, y = cumsum(per) - per/2, label=label))



# You can change the border of each area with the classical parameters:
pie(my_new_df$count , labels =paste0(my_new_df$group, " ", sprintf("%1.f", my_new_df$count*100),"%"), border="white", col=right_colours, cex=1.2)

# geom_text(aes(y = count/5 + c(0, cumsum(count)[-length(count)]), 
#               label = scales::percent(count)), size=5)
pie



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
# mydf[6,"vjust"] =mydf[6,"vjust"]+ 0.05
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



png("Graph_Pie_CHART_FOR_POSTER.png", width = 300, height = 300, units='mm', res = 300)
bp
dev.off()

pdf("Graph_Pie_CHART_FOR_POSTER.pdf", width = 15, height = 10)#, units='mm', res = 300)
bp
dev.off()

        