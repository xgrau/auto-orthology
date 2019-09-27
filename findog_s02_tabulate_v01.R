# libraries
library(igraph)
library(scales)

# input parmeters
edg_fn = "set_raxml.out_ete.sos_0.00.edges"
nod_fn = "set_raxml.out_ete.nodes"
dic_fn = "set_Hsap_names.dict"
man_fn = "set_manualclass.ml.csv"
# input parmeters
# edg_fn = "adar.out_ete.sos_0.00.edges"
# nod_fn = "adar.out_ete.nodes"
# dic_fn = "adar.dict"
sup_th = 10

# load data
edg = read.table(edg_fn, header = T)
nod = read.table(nod_fn, header = T)
dic = read.table(dic_fn, header = T)
man = read.table(man_fn, header = T)

# is node in reference dict?
nod[nod$gene %in% as.character(dic$gene),"ref"]  = "ref"
nod[!nod$gene %in% as.character(dic$gene),"ref"] = "not"
nod = merge(nod,dic, by.x = "gene", by.y = "gene", all.x=T)

# remove edges with suppor under threshold
edg = edg[edg$branch_support >= sup_th,]

# create network
net = graph_from_data_frame(edg, vertices = nod, directed = F)
net_layout = layout_with_mds(net) # MDS layout

# add gene names as labels, colors and sizes for nodes
V(net)$label = as.character(nod$gene)
net_nod_siz = c(ref = 2, not = 1) # node sizes
net_nod_col = c(ref = "blue", not = "slategray4") # node colors
igraph::V(net)$color = net_nod_col[igraph::V(net)$ref]
igraph::V(net)$size  = net_nod_siz[igraph::V(net)$ref]

# plot igraph
pdf(paste(edg_fn,".pdf",sep=""),height=8,width=8)
plot.igraph(net, 
            #vertex.label=ifelse(V(net)$ref == "ref", V(net)$label, NA), 
            vertex.label=V(net)$family,
            vertex.label.family="sans", vertex.frame.color=NA, vertex.label.cex=0.7,
            edge.color = alpha("slategray3",0.5),
            layout=net_layout)
dev.off()

# find components
net_components = components(net)
nod$component  = net_components$membership

nod_ref = nod[!is.na(nod$family),]

for (rei in 1:length(nod_ref$component)) {
  # assign family (same as ref sequence with which it is sharing component)
  nod[nod$component == nod_ref[rei,"component"], "family_inferred"] = nod_ref[rei,"family"]  
}

# compare with manual results
source("../helper_scripts/geneSetAnalysis.R")
pdf(paste(edg_fn,"venns.pdf",sep=""),height=4,width=4)
for (rei in 1:length(nod_ref$component)) {
  man_list = as.character(man[man$family          == nod_ref[rei,"family"],"gene"])
  nod_list = as.character(nod[nod$family_inferred == nod_ref[rei,"family"],"gene"])
  nod_list = nod_list[!is.na(nod_list)]
  
  # plot venn
  ven = venn.two(list1 = nod_list , list2 = man_list, catname1 = "inferred", catname2 = "manual", main = as.character(nod_ref[rei,"family"]))
    
}
dev.off()


# hist(table(net_components$membership))
