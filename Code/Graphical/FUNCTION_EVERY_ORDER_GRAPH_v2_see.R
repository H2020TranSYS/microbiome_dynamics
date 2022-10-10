
# USED in the ../6m or ../9M folder in the "Importing_data_create_ --------


library(dplyr)

checkStrict <- function(f, silent=FALSE) {
  vars <- codetools::findGlobals(f)
  found <- !vapply(vars, exists, logical(1), envir=as.environment(2))
  if (!silent && any(found)) {
    warning("global variables used: ", paste(names(found)[found], collapse=', '))
    return(invisible(FALSE))
  }
  
  !any(found)
}
# 
# plot_graphical_phylum(
#   OTU_RESULT = binaryresult, 
#   label_name_final = "_graphs_Paddy_Covariates_MAGMA_Phylum_shrinking",
#   otu_table = otu_table,
#   Taxonomy_ASVs = Taxonomy_ASVs, 
#   continuous_or_binary = "BINARY",
#   month = "6M",
#   import_layout = "../6M/LAYOUT_PHYLUM_6M_graph.txt"
# )
# 
# plot_graphical_phylum(
#   OTU_RESULT = continuousresult, 
#   label_name_final = "_graphs_Paddy_Covariates_MAGMA_Phylum_shrinking",
#   otu_table = otu_table,
#   Taxonomy_ASVs = Taxonomy_ASVs, 
#   continuous_or_binary = "CONTINUOUS",
#   month = "6M",
#   import_layout = "../6M/LAYOUT_PHYLUM_6M_graph.txt"
# )
# 
# plot_graphical_phylum(
#   OTU_RESULT = continuousresult, 
#   label_name_final = "_graphs_Paddy_Covariates_MAGMA_Phylum_shrinking",
#   otu_table = otu_table,
#   Taxonomy_ASVs = Taxonomy_ASVs, 
#   continuous_or_binary = "CONTINUOUS",
#   month = "9M",
#   import_layout =  "../6M/LAYOUT_PHYLUM_6M_graph.txt"
# )
# 
# 
# plot_graphical_phylum(
#   OTU_RESULT = binaryresult, 
#   label_name_final = "_graphs_Paddy_Covariates_MAGMA_Phylum_shrinking",
#   otu_table = otu_table,
#   Taxonomy_ASVs = Taxonomy_ASVs, 
#   continuous_or_binary = "BINARY",
#   month = "9M",
#   import_layout =  "../6M/LAYOUT_PHYLUM_6M_graph.txt"
# )


plot_graphical_phylum = function(OTU_RESULT, otu_table, Taxonomy_ASVs,continuous_or_binary,level_taxa_order = "Phylum", label_name_final = "_graphs_Paddy_Covariates_MAGMA_Phylum_shrinking", month = "6M", import_layout = F  )
{
  #   """
  #   OTU_RESULTS: or the binary or the continuous results for the OTU FINAL TABLE
  # Originated from 
  # magma_Stool_AllCov <- magma(data = otu_table9M)
  # index = match(magma_Stool_AllCov$opt.lambda,magma_Stool_AllCov$lambda ) 
  # binaryresult = magma_Stool_AllCov$path[[index]]
  # colnames(binaryresult) = rownames(binaryresult) = colnames(otu_table9M)
  # continuousresult = magma_Stool_AllCov$opt.icov
  # colnames(continuousresult) = rownames(continuousresult) = colnames(otu_table9M)
  # similar code
  #   
  #   otu_table: the otu table with the taxas abundance tused only to filter the taxonomy: format individuals on the rows, taxas on the columns
  #   Taxonomy_ASVs: the taxonomy of every individual in the form 
  #          Kingdom             Phylum                  Class                 Order                 Family               Genus           Species
  # 1   k__Bacteria  p__Proteobacteria c__Gammaproteobacteria  o__Enterobacteriales  f__Enterobacteriaceae      g__Escherichia           s__coli
  # 2   k__Bacteria  p__Actinobacteria      c__Actinobacteria  o__Bifidobacteriales  f__Bifidobacteriaceae  g__Bifidobacterium         s__longum
  # 3   k__Bacteria  p__Actinobacteria      c__Actinobacteria  o__Bifidobacteriales  f__Bifidobacteriaceae  g__Bifidobacterium               s__
  # 4   k__Bacteria  p__Actinobacteria      c__Actinobacteria  o__Bifidobacteriales  f__Bifidobacteriaceae  g__Bifidobacterium         s__longum
  # 
  # continuous_or_binary = a string "continuous" or "binary" depending on which data you put in OTU_RESULT that is used to save 
  
  
  # label_name_final = "" the name that we want to give to the graph 
  # month = "6M" : the month that we need for the graph
  # 
  
  # import_layout = F --> if we want to import layout e.g. "../6M/LAYOUT_PHYLUM_6M_graph.txt"
  # N.b we will use the LOO of the month 6 for the phylum graph
  # 
  # Return: nothing, 
  #  
  # plot the graphs
  #   """
  
  level_dataset = data.frame(name = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), rank = seq(1,7))
  
  rank_current_taxa = level_dataset[level_taxa_order == level_dataset$name, "rank"]
  # Kingdom "k__Bacteria"                                                                                                                                                                                                                                                                                                                                                                                                         
  # Phylum  "p__TM7"                                                                                                                                                                                                                                                                                                                                                                                                              
  # Class   "c__TM7-3"                                                                                                                                                                                                                                                                                                                                                                                                            
  # Order   "o__"                                                                                                                                                                                                                                                                                                                                                                                                                 
  # Family  "f__"                                                                                                                                                                                                                                                                                                                                                                                                                 
  # Genus   "g__"                                                                                                                                                                                                                                                                                                                                                                                                                 
  # Species "s__"       
  
  
  Taxonomy_ASVs[2,]
  ASV_final = Taxonomy_ASVs[Taxonomy_ASVs$ASV %in% colnames(otu_table),]
  df = (t(ASV_final))
  dim(df)
  colnames(df) = df[1,]
  df = df[-1,]
  
  level_taxa = unlist(unique(ASV_final %>% select (level_taxa_order)))
  # CREATE PHYLYMS
  
  Taxas_name = colnames(OTU_RESULT)
  
  phy_sum = matrix(data = 0, nrow = length(level_taxa), ncol = length(level_taxa))
  colnames(phy_sum) = rownames(phy_sum) = level_taxa
  phy_cont = phy_sum
  
  # table(df[5,colnames(continuousresult)])
  for (i in 1:nrow(OTU_RESULT)){
    for (j in 1:ncol(OTU_RESULT)){
      phy_sum[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] =phy_sum[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] + abs(OTU_RESULT[i,j])
      phy_cont[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] = phy_cont[df[rank_current_taxa,Taxas_name[i]],df[rank_current_taxa,Taxas_name[j]] ] +1 
    }
  }
  
  taxas_per_phylum = table(df[rank_current_taxa,])
  ordered_taxa_per_phylum = taxas_per_phylum[colnames(phy_cont)]
  phy_cont_less_diag = phy_cont  
  diag(phy_cont_less_diag) = diag(phy_cont) - ordered_taxa_per_phylum 
  
  
  phy_sum_norm = round((phy_sum / phy_cont_less_diag ) * 100, 1)
  
  
  # NA CREATED because 0 / 0 --> there is no signal, superimpose 0 ----------
  phy_sum_norm[is.na(phy_sum_norm) | is.infinite(phy_sum_norm)] = 0
  
  gg = graph_from_adjacency_matrix(
    phy_sum_norm,
    mode = "lower", ## lower traingolar 
    weighted = TRUE,
    diag = TRUE,
    add.colnames = NULL,
    add.rownames = NA
  )
  
  
  # FROM THE rMAGMA::plot_magma() function ----------------------------------
  
  
  
  # V(gg)$color = (unique(df[rank_current_taxa,] ))#factor(seq(1,6))# as.numeric(factor((unique(df[rank_current_taxa,]))) )
  V(gg)$color = as.numeric(factor((unique(df[rank_current_taxa,]))) )
  V.color.factor = factor(unique(df[rank_current_taxa,]))
  V.color = NA
  if (any(!is.na(V.color.factor))) {
    V.color.factor <- as.factor(V.color.factor)
    ll <- length(levels(V.color.factor))
    if (any(is.na(V.color))) {
      if (ll < 10) {
        V.color <- c("red", "green3", "cyan", "blue", 
                     "yellow", "magenta", "whitesmoke", "gray", 
                     "black")
      }
      else {
        V.color <- grDevices::rainbow(length(levels(V.color.factor)))
      }
      igraph::V(gg)$color <- V.color[V.color.factor]
    }
  }
  # V.color.factor = (unique(df[rank_current_taxa,]))
  
  
  set.seed(123)
  ll <- igraph::layout_with_fr(gg)
  plot(gg,
       vertex.size = (diag(phy_sum_norm) + 5),
       edge.width = (E(gg)$weight / 2),
       V.color.factor = df[2,],
       layout = ll, 
       edge.label = E(gg)$weight)
  graphics::legend(x = +1.2, y = +1.3,ifelse(stri_replace_all_regex(V.color.factor, "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor, "[a-z]__", "")), 
                   pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                   cex = 1.3, bty = "n", ncol = 1)
  
  
  # with edge weight
  
  
  if (is.character(import_layout))
  {
    LO = as.matrix(read.table( import_layout ))
  } else
  {
    set.seed(1234)
    LO = layout_nicely(gg)
    write.table(LO, file = paste0("LAYOUT",level_taxa_order, ".txt"))
  }
  
  png(paste0(month,label_name_final,"_",continuous_or_binary,"_","withEDGE_LABEL.png" ), width = 465, height = 225, units='mm', res = 300)
  if (rank_current_taxa >= 4){
    plot(gg,
         vertex.size = (ordered_taxa_per_phylum^(7/10)+10),
         edge.width = (E(gg)$weight / 10),
         # V.color.factor = df[2,],
         layout = LO, 
         edge.label = paste0("  ",round(E(gg)$weight,1)),
         #edge.color="black", 
         vertex.label = ordered_taxa_per_phylum,
         vertex.label.font=2, edge.label.font = 2, edge.label.cex = 1.2, 
         vertex.label.cex = 1 , edge.label.color = "black",vertex.label.color = "black")
    # graphics::legend(x = +1.2, y = +1.3,ifelse(stri_replace_all_regex(V.color.factor, "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor, "[a-z]__", "")), 
    #                  pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
    #                  cex = 1.3, bty = "n", ncol = 1)
    graphics::legend(x = +1.2, y = +1.3,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                     pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                     cex = 1.3, bty = "n", ncol = 1)
  } else {
    plot(gg,
         vertex.size = (ordered_taxa_per_phylum^(4/5)+10),
         edge.width = (E(gg)$weight / 2),
         # V.color.factor = df[2,],
         layout = LO, 
         edge.label = paste0("  ",round(E(gg)$weight,1)),
         #edge.color="black", 
         vertex.label = ordered_taxa_per_phylum,
         vertex.label.font=2, edge.label.font = 2, edge.label.cex = 1.2, 
         vertex.label.cex = 1 , edge.label.color = "black",vertex.label.color = "black")
        graphics::legend(x = +1.2, y = +1.3,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                     pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                     cex = 1.3, bty = "n", ncol = 1)
  }
  dev.off()
  
  # V.color.factor[order(V.color.factor)]

  
    # graphics::legend(x = +1.2, y = +1.3,ifelse(stri_replace_all_regex(V.color.factor, "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor, "[a-z]__", "")), 
                   # pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                   # cex = 1.3, bty = "n", ncol = 1)
  
  # without edge weight label  
  png(paste0(month,label_name_final,"_",continuous_or_binary,"_","NO_EDGE_LABEL.png" ), width = 465, height = 225, units='mm', res = 300)
  if (rank_current_taxa >= 4){
    ### IF THERE ARE TOO MANY COLORS --> smaller
    plot(gg,
         vertex.size = (ordered_taxa_per_phylum^(7/10)+10),
         edge.width = (E(gg)$weight / 10),
         V.color.factor = df[2,],
         layout = LO, 
         edge.label = NA ,# round(E(gg)$weight,1),
         # edge.color="black", 
         vertex.label = ordered_taxa_per_phylum,
         vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
         vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
    graphics::legend(x = +1.2, y = +1.3,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                     pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                     cex = 1.3, bty = "n", ncol = 1)
  } else {
    plot(gg,
         vertex.size = (ordered_taxa_per_phylum^(4/5)+10),
         edge.width = (E(gg)$weight / 2),
         V.color.factor = df[2,],
         layout = LO, 
         edge.label = NA ,# round(E(gg)$weight,1),
         # edge.color="black", 
         vertex.label = ordered_taxa_per_phylum,
         vertex.label.font=2, edge.label.font =4, edge.label.cex = 1.2, 
         vertex.label.cex = 1 , edge.label.color = "black", vertex.label.color = "black")
    graphics::legend(x = +1.2, y = +1.3,ifelse(stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "") == "", "<Unknown>",  stri_replace_all_regex(V.color.factor[order(V.color.factor)], "[a-z]__", "")), 
                     pch = 21, col = V.color, pt.bg = V.color, pt.cex = 2.5, 
                     cex = 1.3, bty = "n", ncol = 1)
  }
  dev.off()
  
  if (!(is.character(import_layout)))
  {
    write.table(LO,file = paste0("LAYOUT_",level_taxa_order,"_",month,"_graph.txt") ) 
  }
}
checkStrict(plot_graphical_phylum)

