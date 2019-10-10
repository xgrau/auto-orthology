# libraries
library(phytools)
library(scales)

#### Input ####
# input parmeters
tau_fn = "orthofinder_Ano14sp/Orthogroups.GeneCount.csv"
nod_fn = "family_SET/set_raxml.out_ete.nodes"
dic_fn = "family_SET/set_Hsap_names.dict"
man_fn = "family_SET/set_manualclass.pb.csv"
# input parmeters
# edg_fn = "family_ADAR/adar.out_ete.sos_0.00.edges"
# nod_fn = "family_ADAR/adar.out_ete.nodes"
# dic_fn = "family_ADAR/adar.dict"
sup_th = 0


# load data
edg = read.table(edg_fn, header = T)
nod = read.table(nod_fn, header = T)
dic = read.table(dic_fn, header = T)
man = read.table(man_fn, header = T)

#### Graph info ####
# is node in reference dict?
nod[nod$gene %in% as.character(dic$gene),"ref"]  = "ref"
nod[!nod$gene %in% as.character(dic$gene),"ref"] = "not"
nod = merge(nod,dic, by.x = "gene", by.y = "gene", all.x=T)

# remove edges with suppor under threshold
edg = edg[edg$branch_support >= sup_th,]

# create network
net = graph_from_data_frame(edg, vertices = nod, directed = F)
net_layout_mds = layout_with_mds(net) # MDS layout
net_layout_nic = layout_components(net) # nice layout

# add gene names as labels, colors and sizes for nodes
V(net)$label = as.character(nod$gene)
net_nod_siz = c(ref = 2, not = 1) # node sizes
net_nod_col = c(ref = "blue", not = "slategray4") # node colors
igraph::V(net)$color = net_nod_col[igraph::V(net)$ref]
igraph::V(net)$size  = net_nod_siz[igraph::V(net)$ref]

# plot igraph
pdf(paste(edg_fn,".network.pdf",sep=""),height=5,width=5)
plot.igraph(net, 
            vertex.label=V(net)$family,
            vertex.label.family="sans", vertex.frame.color=NA, vertex.label.cex=0.7,
            edge.color = alpha("slategray3",0.5),
            layout=net_layout_mds)
plot.igraph(net, 
            vertex.label=V(net)$family,
            vertex.label.family="sans", vertex.frame.color=NA, vertex.label.cex=0.7,
            edge.color = alpha("slategray3",0.5),
            layout=net_layout_nic)
dev.off()



#### Assign families ####
# find components

# assign family (same as ref sequence with which it is sharing component)
# PROBLEMATIC: SAME SEQ CAN SHARE COMPONENT WITH MORE THAN ONE REF
# for (rei in 1:length(nod_ref$component)) {
#   nod[nod$component == nod_ref[rei,"component"], "family_inferred"] = nod_ref[rei,"family"]  
# }

# identify all families each seq is linked to (by orthology)
for (noi in 1:nrow(nod)) {
  noi_bool = nod[noi,"gene"] == edg$in_gene | nod[noi,"gene"] == edg$out_gene
  noi_comp_elements = unique(c(as.character(edg[noi_bool,c("in_gene")]),as.character(edg[noi_bool,c("out_gene")])))
  noi_comp_refvec   = as.character(dic[dic$gene %in% noi_comp_elements,"family"])
  noi_comp_refstr   = paste(noi_comp_refvec, collapse = ',')
  nod[noi,"family_inferred"] = noi_comp_refstr
}

nod$family_inferred_factor= as.factor(nod$family_inferred)
factor_colors = c("slategray4",
                  "red1","red4","purple1","purple4",
                  "blue1","blue4","olivedrab1","olivedrab4",
                  "darkgreen","orange1","orange3",
                  "cyan4","cyan2","gold","limegreen",
                  "violetred1","violetred4",
                  "springgreen1","springgreen4",
                  "slateblue2","slateblue4","sienna2","sienna4",
                  "paleturquoise1","paleturquoise4")
faminf_colors = factor_colors[nod$family_inferred_factor]

# replot, with lots of colors
V(net)$color = faminf_colors
pdf(paste(edg_fn,".network_colorfams.pdf",sep=""),height=5,width=5)
plot.igraph(net, 
            vertex.label=V(net)$family,vertex.size=2,
            vertex.label.family="sans", vertex.frame.color=NA, vertex.label.cex=0.7,
            edge.color = alpha("slategray3",0.5),
            layout=net_layout_mds)
plot.igraph(net, 
            vertex.label=V(net)$family,vertex.size=2,
            vertex.label.family="sans", vertex.frame.color=NA, vertex.label.cex=0.7,
            edge.color = alpha("slategray3",0.5),
            layout=net_layout_nic)
legend("topright", legend = levels(nod$family_inferred_factor), col=factor_colors, pch=20, cex=0.3, bty = "n")
dev.off()


#### Compare manual ####
# compare with manual results
source("helper_scripts/geneSetAnalysis.R")
pdf(paste(edg_fn,".venns.pdf",sep=""),height=4,width=4)
for (rei in 1:nrow(dic)) {
  man_list = as.character(man[man$family          == dic[rei,"family"],"gene"])
  nod_list = as.character(nod[nod$family_inferred == dic[rei,"family"],"gene"])
  nod_list = nod_list[!is.na(nod_list)]
  
  # plot venn
  # TODO: report lists of intersections, disjoint, etc. (in ven object!)
  ven = venn.two(list1 = nod_list , list2 = man_list, catname1 = "inferred", catname2 = "manual", main = as.character(dic[rei,"family"]))
  
}
dev.off()


# hist(table(net_components$membership))
