library(vegan)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(xlsx)
library(ggpubr)

data_path   = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Data"
result_path = "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Results/Alpha_beta_diversity"
graph_path =  "C:/Users/fmelo/Documents/GitHub/microbiome_dynamics/Graphs/Alpha_beta_diversity"



# FUNCTIONS ---------------------------------------------------------------
checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}

calculate_violin_plot_test_del = function(dataset, matching_data, ext_pheno = "Delivery", name_g,dataset2 = NULL, low_l = 0.3, high_l = 3){
  if (is.null(dataset2)){ ## ONLY 1 dataset --> 1 timepoint
    vSHAN = vegan::diversity(dataset, index = "shannon", MARGIN = 1, base = exp(1))
    if (ext_pheno == "Delivery"){ ## DIFFERENT per phenotype
      is_variable = ifelse(names(vSHAN) %in% matching_data$Child[matching_data$delivery_type.x == "Vaginal"],"Vaginal", "C-section")
    } else if (ext_pheno == "Diet"){
      is_same = matching_data$Child[matching_data$Diet_TP.x == matching_data$Diet_TP.y] %>% na.omit()
      is_variable = ifelse(names(vSHAN) %in% is_same,"Persist.", "Non persist.")
    }
    
  } else { ## two datasts
    vSHAN1 = vegan::diversity(dataset, index = "shannon", MARGIN = 1, base = exp(1))
    vSHAN2 = vegan::diversity(dataset2, index = "shannon", MARGIN = 1, base = exp(1))
    if (ext_pheno == "Delivery"){
      is_variable1 = ifelse(names(vSHAN1) %in% matching_data$Child[matching_data$delivery_type.x == "Vaginal"],"Vaginal", "C-section")
      is_variable2 = ifelse(names(vSHAN2) %in% matching_data$Child[matching_data$delivery_type.x == "Vaginal"],"Vaginal", "C-section")
    } else if (ext_pheno == "Diet"){
      is_same = matching_data$Child[matching_data$Diet_TP.x == matching_data$Diet_TP.y] %>% na.omit()
      is_variable1 = ifelse(names(vSHAN1) %in% is_same,"Persist.", "Non persist.")
      is_variable2 = ifelse(names(vSHAN2) %in% is_same,"Persist.", "Non persist.")
    }
    vSHAN = c(vSHAN1, vSHAN2) ; is_variable = c(is_variable1, is_variable2)
  }
  # vSHAN = vegan::diversity(dataset, index = "shannon", MARGIN = 1, base = exp(1))
  # is_variable = ifelse(names(vSHAN) %in% matching_data$Child[matching_data$delivery_type.x == "Vaginal"],"Vaginal", "C-section")
  vplotdf = data.frame("Alpha_div" = vSHAN, "Type_of_Delivery" = is_variable) 
  
  graph_name = ifelse(ext_pheno == "Delivery", "Mode of Delivery", "Diet")
  
  g1 = ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Type_of_Delivery)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
    geom_violin() + ggtitle(paste0("Violin plot ",name_g, " Shannon diversity")) +
    stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", size=2, color="gray") +
    # DOT IS THE MEAN ; the extremities are + and - SD
    # geom_boxplot(width = 0.2)+
    theme(plot.title = element_text(hjust = 0.5)) + 
    theme_minimal(base_size = 17, base_line_size = 1.1)+
    labs(y= expression(alpha*" diversity Shannon"), x = graph_name, fill = "")+ #, x = "x axis name")
    coord_cartesian(ylim = c(low_l, high_l)) +
    
    # theme_dark()+
    # theme_classic(base_size = 17, base_line_size = 1.1)+
    theme(axis.text=element_text(face="bold"),
          axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
  
  if (ext_pheno == "Delivery"){
    (w4 = wilcox.test(vSHAN[is_variable == "Vaginal"], vSHAN[is_variable != "Vaginal"], paired = F))
  } else if (ext_pheno == "Diet"){
    (w4 = wilcox.test(vSHAN[is_variable == "Persist."], vSHAN[is_variable != "Persist."], paired = F))
    
  }
  return(list(graph = g1, test = w4, summ =vplotdf))
}
checkStrict(calculate_violin_plot_test_del)

# ALPHa and BETA diversity  -----------------------------------------------
# setwd("C:/Users/fmelo/Desktop/Backup_Federico/Maastricht_data/LucKi_cohort_data/MAGMA_APPROACH/OTUS_node_analysis_difference")
setwd(data_path)

OTU_TABLE_MAGMA6M <- read.delim("MAGMA_data/6M/OTU_TABLE_MAGMA6M_69.tsv")
# 69 x 95
OTU_TABLE_MAGMA9M <- read.delim("MAGMA_data/9M/OTU_TABLE_MAGMA9M_69.tsv")
# OTU_TABLE_MAGMA6M = read.table(file = "OTU_TABLE_69_SUB_6M.tsv")
# OTU_TABLE_MAGMA6M_aa = OTU_TABLE_MAGMA6M[rownames(OTU_TABLE_MAGMA6M_nn),]
# OTU_TABLE_MAGMA9M = read.table(file = "OTU_TABLE_69_SUB_9M.tsv")

head(OTU_TABLE_MAGMA6M)
OTU_TABLE_MAGMA6M[1:5,1:5]


Taxonomy_ASVs <- read.csv("Taxonomy_ASVs.csv")
LuckiMap_anonym_long_format_6m <- read.delim("LuckiMap_anonym_long_format_6m.txt")
LuckiMap_anonym_long_format_9m <- read.delim("LuckiMap_anonym_long_format_9m.txt")
merged_Lucki = merge(LuckiMap_anonym_long_format_6m, LuckiMap_anonym_long_format_9m,by =  "Child")



# USE the matched name CHILD__ instead of LUC -----------------------------

matched = match( rownames(OTU_TABLE_MAGMA6M), merged_Lucki$Sample.name.x)
names_6m = merged_Lucki$Child[matched]
OTU_TABLE_MAGMA6M_nn = OTU_TABLE_MAGMA6M
rownames(OTU_TABLE_MAGMA6M_nn)= names_6m


matched = match( rownames(OTU_TABLE_MAGMA9M), merged_Lucki$Sample.name.y)
names_9m = merged_Lucki$Child[matched]
OTU_TABLE_MAGMA9M_nn = OTU_TABLE_MAGMA9M
rownames(OTU_TABLE_MAGMA9M_nn)= names_9m




# Divided by timepoint ----------------------------------------------------



vSHAN_DIV6M = vegan::diversity(OTU_TABLE_MAGMA6M_nn, index = "shannon", MARGIN = 1, base = exp(1))
vSHAN_DIV9M = vegan::diversity(OTU_TABLE_MAGMA9M_nn, index = "shannon", MARGIN = 1, base = exp(1))

vSHAN_DIV6M
vplotdf = data.frame("Alpha_div" = c(vSHAN_DIV6M, vSHAN_DIV9M), "Timepoint" = c(rep("6M", 69),rep("9M", 69))) 
# Basic violin plot
p <- ggplot(vplotdf, aes(x=Timepoint, y=Alpha_div)) + 
  geom_violin()
p

ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Timepoint)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
  geom_violin() + ggtitle("Violin plot Shannon diversity") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange",shape=3, size=3, color="gray") +
  # geom_boxplot(width = 0.2)+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_gray(base_size = 17, base_line_size = 1.1)+
  labs(y= expression(alpha*" diversity Shannon"), x = expression("Timepoint"), fill = "")+ #, x = "x axis name")
  # theme_dark()+
  # theme_classic(base_size = 17, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")


# 2 graphs --> Mult = 1 and Mult = 2 [ the standard deviations] -----------



# mult = 1 ----------------------------------------------------------------

setwd(result_path)

vplotdf%>%
  group_by(Timepoint)%>% 
  summarise(Mean=mean(Alpha_div), Max=max(Alpha_div), Min=min(Alpha_div), Median=median(Alpha_div), Std=sd(Alpha_div), IQR = IQR(Alpha_div)) %>%
  write.xlsx(.,file = paste0("Statistic_alpha", ".xlsx"),
             sheetName = "Per_timepoint", append = FALSE)


setwd(graph_path)

low_l = 0.3 ; high_l = 3
v_age = ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Timepoint)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
  geom_violin() + ggtitle("Violin plot Shannon diversity") +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult = 1), geom="pointrange", size=2, color="gray") +
  # DOT IS THE MEAN ; the extremities are + and - SD
  # geom_boxplot(width = 0.2)+
  theme(plot.title = element_text(hjust = 0.5)) + 
  theme_minimal(base_size = 17, base_line_size = 1.1)+
  coord_cartesian(ylim = c(low_l, high_l)) +
  labs(y= expression(alpha*" diversity Shannon"), x = "Timepoint", fill = "")+ #, x = "x axis name")
  # theme_dark()+
  # theme_classic(base_size = 17, base_line_size = 1.1)+
  theme(axis.text=element_text(face="bold"),
        axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
# png("Violin_Shannon_alphadiv.png", width = 300, height = 300, units='mm', res = 300)
# print(v_age)
# dev.off()
# png("Violin_Shannon_alphadiv_ds.png", width = 465, height = 225, units='mm', res = 300)
# print(v_age)
# dev.off()



# Mult = 2 ----------------------------------------------------------------



# png("Violin_Shannon_alphadiv_2sd.png", width = 300, height = 300, units='mm', res = 300)
# ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Timepoint)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
#   geom_violin() + ggtitle("Violin plot Shannon diversity") +
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult = 2), geom="pointrange", size=2, color="gray") +
#   # DOT IS THE MEAN ; the extremities are + and - SD
#   # geom_boxplot(width = 0.2)+
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   theme_minimal(base_size = 17, base_line_size = 1.1)+
#   labs(y= expression(alpha*" diversity Shannon"), x = expression("Timepoint"), fill = "")+ #, x = "x axis name")
#   # theme_dark()+
#   # theme_classic(base_size = 17, base_line_size = 1.1)+
#   theme(axis.text=element_text(face="bold"),
#         axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
# dev.off()
# png("Violin_Shannon_alphadiv_ds_2sd.png", width = 465, height = 225, units='mm', res = 300)
# ggplot(data.frame(Trait=vplotdf$Alpha_div, Cluster=as.factor(vplotdf$Timepoint)), aes(x=Cluster, y=Trait, fill = Cluster)) + 
#   geom_violin() + ggtitle("Violin plot Shannon diversity") +
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult = 2), geom="pointrange", size=2, color="gray") +
#   # DOT IS THE MEAN ; the extremities are + and - SD
#   # geom_boxplot(width = 0.2)+
#   theme(plot.title = element_text(hjust = 0.5)) + 
#   theme_minimal(base_size = 17, base_line_size = 1.1)+
#   coord_cartesian(ylim = c(low_l, high_l)) +
#   
#   labs(y= expression(alpha*" diversity Shannon"), x = expression("Timepoint"), fill = "")+ #, x = "x axis name")
#   # theme_dark()+
#   # theme_classic(base_size = 17, base_line_size = 1.1)+
#   theme(axis.text=element_text(face="bold"),
#         axis.title=element_text(face="bold"), title = element_text(face = "bold",),text = element_text(size = 35),legend.position="none")
# dev.off()
# 



# TEST --------------------------------------------------------------------

# w1 = wilcox.test(Alpha_div ~ Timepoint, data=vplotdf) 
# w1_bis = wilcox.test(vSHAN_DIV6M,vSHAN_DIV9M) 
(w3 = wilcox.test(vSHAN_DIV6M,vSHAN_DIV9M, paired = T)  )

## DIFFERENCE in POSITION --> fixed with SORT
sort(names(vSHAN_DIV6M)) == sort(names(vSHAN_DIV9M))
vSHAN_DIV6M[sort(names(vSHAN_DIV6M))]; vSHAN_DIV9M[sort(names(vSHAN_DIV9M))]
(w4 = wilcox.test(vSHAN_DIV6M[sort(names(vSHAN_DIV6M))], vSHAN_DIV9M[sort(names(vSHAN_DIV9M))], paired = T))
## CORRECT one


setwd(result_path)
saveRDS(w4, "TEST_WILCOXON_age.rds")
setwd(graph_path)

# > (w4 = wilcox.test(vSHAN_DIV6M[sort(names(vSHAN_DIV6M))], vSHAN_DIV9M[sort(names(vSHAN_DIV9M))], paired = T))
# 
# Wilcoxon signed rank test with continuity correction
# 
# data:  vSHAN_DIV6M[sort(names(vSHAN_DIV6M))] and vSHAN_DIV9M[sort(names(vSHAN_DIV9M))]
# V = 284, p-value = 3.418e-08
# alternative hypothesis: true location shift is not equal to 0



# Violin plot per delivery ------------------------------------------------



v_6m = calculate_violin_plot_test_del(dataset = OTU_TABLE_MAGMA6M_nn, matching_data = merged_Lucki, name_g = "6M")
v_9m = calculate_violin_plot_test_del(dataset = OTU_TABLE_MAGMA9M_nn, matching_data = merged_Lucki, name_g = "9M")
v_6m$graph
v_9m$graph




colnames(OTU_TABLE_MAGMA6M_nn) == colnames(OTU_TABLE_MAGMA9M_nn)
#OK 
v_gen = calculate_violin_plot_test_del(dataset = OTU_TABLE_MAGMA6M_nn,dataset2 = OTU_TABLE_MAGMA9M_nn, matching_data = merged_Lucki, name_g = "both tp")
v_gen$graph

v_6m$test
v_9m$test
v_gen$test


# png("Alpha_diversity_6M_per_delivery.png", width = 300, height = 300, units='mm', res = 300)
# v_6m$graph
# dev.off()
# png("Alpha_diversity_9M_per_delivery.png", width = 300, height = 300, units='mm', res = 300)
# v_9m$graph
# dev.off()
# png("Alpha_diversity_bothtp_per_delivery.png", width = 300, height = 300, units='mm', res = 300)
# v_gen$graph
# dev.off()




# Same Non persist.s ---------------------------------------------------------


## ELIMINATE MISSING diet in both 
missing_diet = merged_Lucki[is.na(merged_Lucki$Diet_TP.x),]
OTU_TABLE_MAGMA6M_nn_diet = OTU_TABLE_MAGMA6M_nn[!(rownames(OTU_TABLE_MAGMA6M_nn) %in% missing_diet$Child),]
OTU_TABLE_MAGMA9M_nn_diet = OTU_TABLE_MAGMA9M_nn[!(rownames(OTU_TABLE_MAGMA9M_nn) %in% missing_diet$Child),]


vd_6m = calculate_violin_plot_test_del(dataset = OTU_TABLE_MAGMA6M_nn_diet, matching_data = merged_Lucki, ext_pheno = "Diet", name_g = "6M")
vd_9m = calculate_violin_plot_test_del(dataset = OTU_TABLE_MAGMA9M_nn_diet, matching_data = merged_Lucki, ext_pheno = "Diet", name_g = "9M")
vd_gen = calculate_violin_plot_test_del(dataset = OTU_TABLE_MAGMA6M_nn_diet, matching_data = merged_Lucki, ext_pheno = "Diet", name_g = "both tp", dataset2 = OTU_TABLE_MAGMA9M_nn_diet)


# png("Alpha_diversity_6M_per_diet.png", width = 300, height = 300, units='mm', res = 300)
# vd_6m$graph
# dev.off()
# png("Alpha_diversity_9M_per_diet.png", width = 300, height = 300, units='mm', res = 300)
# vd_9m$graph
# dev.off()
# png("Alpha_diversity_bothTP_per_diet.png", width = 300, height = 300, units='mm', res = 300)
# vd_gen$graph
# dev.off()


vd_6m$test
vd_9m$test
vd_gen$test

png("Alpha_diversity_arranged.png", width = 465, height = 225, units='mm', res = 300)
ggarrange(v_age + ggtitle(""), v_gen$graph+ ggtitle(""), vd_gen$graph+ ggtitle(""),  nrow = 1)
dev.off()
ll = list("l6M_delivery" = v_6m,"l9M_delivery" = v_9m, "bothTP_delivery" = v_gen,
          "l6M_diet" = vd_6m,   "l9M_diet" = vd_9m,    "bothTP_diet" = vd_gen)


setwd(result_path)
saveRDS(ll, file = "list_ALPHA_graph_test.rds")

v_gen$summ%>%
  group_by(Type_of_Delivery)%>% 
  summarise(Mean=mean(Alpha_div), Max=max(Alpha_div), Min=min(Alpha_div), Median=median(Alpha_div), Std=sd(Alpha_div), IQR = IQR(Alpha_div)) %>%
  write.xlsx(.,file = paste0("Statistic_alpha", ".xlsx"),
             sheetName = "Per_delivery", append = TRUE)

vd_gen$summ%>%
  group_by(Type_of_Delivery)%>% 
  summarise(Mean=mean(Alpha_div), Max=max(Alpha_div), Min=min(Alpha_div), Median=median(Alpha_div), Std=sd(Alpha_div), IQR = IQR(Alpha_div)) %>%
  write.xlsx(.,file = paste0("Statistic_alpha", ".xlsx"),
             sheetName = "Per_diet", append = TRUE)

# Beta diversity ----------------------------------------------------------



library(compositions)
library(mixOmics)


OTU_TABLE_MAGMA6M_CLR = data.frame( clr( OTU_TABLE_MAGMA6M_nn ))
OTU_TABLE_MAGMA9M_CLR = data.frame( clr( OTU_TABLE_MAGMA9M_nn ))

# OTU_TABLE_MAGMA6M_CLR = read.table(file = "OTU_TABLE_69_SUB_CLR_6M.tsv")
# OTU_TABLE_MAGMA9M_CLR = read.table(file = "OTU_TABLE_69_SUB_CLR_9M.tsv")
## compositional trasformation, if I sum over an individual, the sum = 0 
round(apply(OTU_TABLE_MAGMA9M_CLR,1,sum),2) ; round(apply(OTU_TABLE_MAGMA6M_CLR,1,sum),2)


LuckiMap_6M = LuckiMap_anonym_long_format_6m
LuckiMap_9M = LuckiMap_anonym_long_format_9m
Children =  merge(LuckiMap_6M, LuckiMap_9M, by = "Child")
Children_6M = Children[Children$Child %in% rownames(OTU_TABLE_MAGMA6M_CLR),]
Children_9M = Children[Children$Child %in% rownames(OTU_TABLE_MAGMA9M_CLR),]



# SAME ORDER 6m and 9m ----------------------------------------------------


OTU_TABLE_MAGMA6M_CLR = OTU_TABLE_MAGMA6M_CLR[Children_6M$Child, ]
OTU_TABLE_MAGMA9M_CLR = OTU_TABLE_MAGMA9M_CLR[Children_9M$Child, ]



PCA_6M_CLR = mixOmics::pca(OTU_TABLE_MAGMA6M_CLR)
PCA_6M_CLR
plotIndiv(PCA_6M_CLR, ind.names = F)
plotIndiv(PCA_6M_CLR)

plotVar(PCA_6M_CLR)
plot(PCA_6M_CLR)

PCA_9M_CLR = mixOmics::pca(OTU_TABLE_MAGMA9M_CLR)
PCA_9M_CLR
plotIndiv(PCA_9M_CLR, ind.names = F)
plotIndiv(PCA_9M_CLR)

plotVar(PCA_9M_CLR)
plot(PCA_9M_CLR)



# NOVEL PENDERS asked -----------------------------------------------------



round(apply(OTU_TABLE_MAGMA9M_CLR,1,sum),2) ; round(apply(OTU_TABLE_MAGMA6M_CLR,1,sum),2)

# OTU_TABLE_MAGMA6M_CLR = OTU_TABLE_MAGMA6M_CLR[Children_6M$Child, ]
# OTU_TABLE_MAGMA9M_CLR = OTU_TABLE_MAGMA9M_CLR[Children_9M$Child, ]

OTU_TABLE_MAGMA9M_CLR

OTUS = rbind(OTU_TABLE_MAGMA6M_CLR,OTU_TABLE_MAGMA9M_CLR)
PCA_CLR = mixOmics::pca(OTUS)

setwd(graph_path)
# 
# png("Beta_diversity_per_timepoint.png", width = 300, height = 300, units='mm', res = 300)
# p_age = plotIndiv(PCA_CLR, ind.names = FALSE,
#                   group = c(rep("6M",69), rep("9M",69)),
#                   col.per.group = brewer.pal(n = 4, name = "Set2")[c(2,3)],
#                   pch = 16,#rep(as.numeric(as.factor(Children_6M$delivery_type.x))+15,2),
#                   legend = TRUE, title = "Timepoint",# 'Beta diversity: PCA comp 1 - 2',
#                   legend.title = '', legend.title.pch = 'Diet', legend.position = "top",cex = 3
#                   , size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)
# 
# dev.off()


# per Delivery type -------------------------------------------------------

is_vaginal_6m = rownames(OTU_TABLE_MAGMA6M_CLR) %in% Children$Child[Children$delivery_type.x == "Vaginal"]
is_vaginal_9m = rownames(OTU_TABLE_MAGMA9M_CLR) %in% Children$Child[Children$delivery_type.x == "Vaginal"]

rownames(OTU_TABLE_MAGMA6M_CLR) == rownames(OTU_TABLE_MAGMA9M_CLR)
colnames(OTU_TABLE_MAGMA6M_CLR) == colnames(OTU_TABLE_MAGMA9M_CLR)

is_vag = c(is_vaginal_6m, is_vaginal_9m)
# png("Beta_diversity_per_delivery.png", width = 300, height = 300, units='mm', res = 300)
p_del = plotIndiv(PCA_CLR, ind.names = FALSE,
                  group = ifelse(is_vag, "Vaginal", "C-section"),
                  col.per.group = brewer.pal(n = 4, name = "Set2")[c(2,3)],
                  pch = 16,#rep(as.numeric(as.factor(Children_6M$delivery_type.x))+15,2),
                  legend = TRUE, title = "Mode of delivery", #'Beta diversity: PCA comp 1 - 2',
                  legend.title = '', legend.title.pch = 'Diet', legend.position = "top",cex = 3
                  , size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)
# dev.off()

plotIndiv(PCA_CLR, ind.names = T)


# per Diet ----------------------------------------------------------------

missing_diet = Children[is.na(Children$Diet_TP.x),]
OTU_TABLE_MAGMA6M_CLR_diet = OTU_TABLE_MAGMA6M_CLR[!(rownames(OTU_TABLE_MAGMA6M_CLR) %in% missing_diet$Child),]
OTU_TABLE_MAGMA9M_CLR_diet = OTU_TABLE_MAGMA9M_CLR[!(rownames(OTU_TABLE_MAGMA9M_CLR) %in% missing_diet$Child),]



# SINCE IT IS REPEATED, the second row names = row name +1 ----------------

na_child = Children_6M[is.na(Children_6M$Diet_TP.x), , drop = F]
na_child = rbind(na_child, na_child)
na_child$Child[2] = paste0(na_child$Child[2], "1")
## work only if just one MISSING

PCA_CLR_DIET = PCA_CLR
PCA_CLR_DIET$X = PCA_CLR_DIET$X[!rownames(PCA_CLR_DIET$X) %in% na_child$Child,]

PCA_CLR_DIET$names$sample = PCA_CLR_DIET$names$sample[! PCA_CLR_DIET$names$sample %in% na_child$Child]
PCA_CLR_DIET$variates$X = PCA_CLR_DIET$variates$X[!rownames(PCA_CLR_DIET$variates$X) %in% na_child$Child,]

dim(PCA_CLR_DIET$X)


is_same = Children$Child[Children$Diet_TP.x == Children$Diet_TP.y] %>% na.omit()

is_same_6m = rownames(OTU_TABLE_MAGMA6M_CLR_diet) %in% is_same
is_same_9m = rownames(OTU_TABLE_MAGMA9M_CLR_diet) %in% is_same

rownames(OTU_TABLE_MAGMA6M_CLR) == rownames(OTU_TABLE_MAGMA9M_CLR)
colnames(OTU_TABLE_MAGMA6M_CLR) == colnames(OTU_TABLE_MAGMA9M_CLR)

is_same_f = c(is_same_6m, is_same_9m)
# png("Beta_diversity_per_diet.png", width = 300, height = 300, units='mm', res = 300)
p_diet = plotIndiv(PCA_CLR_DIET, ind.names = FALSE,
                   group = ifelse(is_same_f, "Persist.", "Non persist."),
                   col.per.group = brewer.pal(n = 4, name = "Set2")[c(2,3)],
                   pch = 16,#rep(as.numeric(as.factor(Children_6M$delivery_type.x))+15,2),
                   legend = TRUE, title = "Diet", #'Beta diversity: PCA comp 1 - 2',
                   legend.title = '', legend.title.pch = 'Diet',cex = 3, legend.position = "top",
                   size.title = 30, point.lwd = 2, size.legend = 25, size.legend.title = 30, size.xlabel = 25, size.ylabel = 25, size.axis = 18)

# dev.off()
plotIndiv(PCA_CLR_DIET, ind.names = T)

png("Beta_diversity_arranged.png", width = 465, height = 225, units='mm', res = 300)
ggarrange(p_age$graph, p_del$graph, p_diet$graph, nrow = 1)
dev.off()


