# libraries
library(igraph)
library(scales)

setwd("/home/xavi/Documents/auto-orthology/possom/")
#### Input ####
# input parmeters
dic_fn = "set_family/set_Hsap_names.dict"
man_fn = "set_family/set_manualclass.ml.csv"
aut_fn = "set_family_groups/group0.orthologs.csv"
out_fn = "set_family_groups/group0.orthologs_overlap"

aut = read.table(aut_fn, header = T)
dic = read.table(dic_fn, header = T)
man = read.table(man_fn, header = T)

#### Graph info ####
# is node in reference dict?
aut[aut$node %in% as.character(dic$gene),"ref"]  = "ref"
aut[!aut$gene %in% as.character(dic$gene),"ref"] = "not"
dic = merge(dic,aut, by.x = "gene", by.y = "node")
# aut = merge(aut,dic, by.x = "node", by.y = "gene", all.x=T)

#### Compare manual ####
# compare with manual results
source("../helper_scripts/geneSetAnalysis.R")

for (clustertype in c("cluster","cluster_ref")) {
  pdf(paste(out_fn,".out_",clustertype,"_venns.pdf",sep=""),height=4,width=4)
  
  for (rei in 1:nrow(dic)) {
    
    rei_family  = as.character(dic[rei,"family"])
    rei_cluster = as.character(dic[rei,clustertype])
    rei_gene    = as.character(dic[rei,"gene"])
    
    man_list = as.character(man[man$family == rei_family,"gene"])
    aut_list = as.character(aut[aut[,clustertype] == rei_cluster,"node"])
    
    # plot venn
    # TODO: report lists of intersections, disjoint, etc. (in ven object!)
    ven = venn.two(list1 = aut_list , list2 = man_list, 
                   catname1 = rei_cluster, 
                   catname2 = "manual",
                   main = paste(rei_family, rei_gene))
    
  }
  
  dev.off()
}

# hist(table(net_components$membership))
