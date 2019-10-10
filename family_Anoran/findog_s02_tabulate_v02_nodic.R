# libraries
library(igraph)
library(scales)
library(stringr)

setwd("/home/xavi/Documents/auto-orthology/family_Anoran/")

#### Input ####
# input parmeters
edg_fn = "OG0000000.out_ete.sos_0.00.edges"
nod_fn = "OG0000000.out_ete.nodes"
dic_sp = "Anogam"

# input parmeters
# edg_fn = "family_ADAR/adar.out_ete.sos_0.00.edges"
# nod_fn = "family_ADAR/adar.out_ete.nodes"

sup_th = 0


# load data
edg = read.table(edg_fn, header = T)
nod = read.table(nod_fn, header = T)

#### Graph info ####
# is node in reference dict?
nod[nod$species == dic_sp, "family"] = nod[nod$species == dic_sp, "gene"]
nod[nod$species == dic_sp, "ref"] = "ref"
nod[nod$species != dic_sp, "ref"] = "not"
nod$family = str_remove_all(nod$family, paste(dic_sp,"_",sep=""))

# remove edges with suppor under threshold
edg = edg[edg$branch_support >= sup_th,]

# create network
net = graph_from_data_frame(edg, vertices = nod, directed = F)
net_layout_mds = layout_with_mds(net) # MDS layout
net_layout_nic = layout_components(net) # nice layout

# add gene names as labels, colors and sizes for nodes
V(net)$label  = as.character(nod$gene)
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
  noi_comp_refvec   = na.omit(as.character(nod[nod$gene %in% noi_comp_elements,"family"]))
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
                  "paleturquoise1","paleturquoise4",
                  "turquoise3","cyan")
faminf_colors = factor_colors[nod$family_inferred_factor]

# replot, with lots of colors
V(net)$color = faminf_colors
pdf(paste(edg_fn,".network_colorfams.pdf",sep=""),height=5,width=5)
plot.igraph(net, 
            vertex.label=V(net)$family,vertex.size=2,
            vertex.label.family="sans", vertex.frame.color=NA, vertex.label.cex=0.5,
            edge.color = alpha("slategray3",0.5),
            layout=net_layout_nic)
plot.new()
legend("topright", legend = levels(nod$family_inferred_factor), col=factor_colors, pch=20, cex=0.3, bty = "n")
dev.off()


# hist(table(net_components$membership))

