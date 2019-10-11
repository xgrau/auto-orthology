# libraries
library(ape)
library(phytools)
library(scales)
library(stringr)
library(reshape2)

#### Input ####
# input parmeters
tau_fn = "out_fas.csv"
phy_fn = "tree_dists.iqt.treefile"
ref_sp = "Anogam"


#### Phylogeny ####
phy = ape::read.tree(phy_fn) 
phy_unrooted = ape::unroot(phy)
write.tree(phy = phy_unrooted, file = paste(phy_fn,"unroot", sep="."))

#### Orthogroups ####
# load orthogroups
tau         = read.table(tau_fn, header = T, sep = "\t")
tau         = tau[!is.na(tau$cluster),]
tau$node    = as.character(tau$node)
tau$species = stringr::str_extract(tau$node, pattern = "^[^_]*")

# gene counts per orthogroup
tac = aggregate(tau$cluster, by = list(tau$species,tau$og_cluster), FUN=length)
colnames(tac) = c("species", "og_cluster","num_genes")
# wide table
taw = dcast(tac, og_cluster ~ species, value.var="num_genes")
rownames(taw) = taw$og_cluster
taw = taw[,phy$tip.label]
# binarise table (presence/absence)
taw[is.na(taw)] = 0
taw[taw>1]      = 1

# identify reference species
lis_ogs = rownames(taw)
lis_ref = tau[tau$og_cluster %in% lis_ogs ,"node"] # TODO: this should identify one ref sequence per orthogroup, from the ref sps




#### Loop ogs ####

n = sample(1:nrow(taw), 1)


# phy$edge.length = rep(0.01,length(phy$edge)/2)
# n = sample(1:nrow(taw), 1)
n = 13662
print(n)
phy_states = factor(as.numeric(taw[n,]))
# levels(phy_states) = c(1,0)
names(phy_states) = colnames(taw)

# ancestral reconstruction
ans = ape::ace(phy_states, phy, type = "d", method = "ML")

# plot
plot(phy, type = "phy", use.edge.length = FALSE, label.offset = 1, main=rownames(taw)[n])
co = c("lightblue", "blue")
tiplabels(pch = 21, bg = co[as.numeric(phy_states)], cex = 2, adj = 1)
nodelabels(pie = ans$lik.anc, piecol = co, cex = 0.75)
nodelabels(signif(ans$lik.anc[,"1"],3), cex = 0.75, col="red", frame="none")
# http://www.phytools.org/eqg2015/asr.html



stop("ara")

#### HOW TO OBTAIN GENOME-WIDE TRANSITION RATES?
Q = matrix(c(-1,1,1,-1),2,2) # transition matrix
colnames(Q) = c(0,1)
rownames(Q) = c(0,1)

phy_sim     = sim.history(phy,Q)
phy_sim_map = make.simmap(phy, phy_states, nsim = 100)
phy_sim_den = densityMap(phy_sim_map,lwd=3,outline=TRUE)








stop("AREA")






# phy$node.states = phy_states
dotTree(phy, phy_states)
# plotTree.barplot(phy, phy_states, args.barplot = list(col="slategray",border=NA))


phy$edge.length =0.1
ace(taw[1:1000,], phy, type = "discrete", method = "ML")
fastAnc(phy, phy_states, vars=F)


# phy_states = t(taw[c(1:10),])
# dotTree(phy, phy_states)
# plotTree.barplot(phy, phy_states)
# phylo.heatmap(phy, phy_states)
# getStates(phy)

phy = make.simmap(phy,phy_states,nsim=100)

plotTree(phy, ftype="i")


Q<-matrix(c(-1,1,1,-1),2,2)
rownames(Q)<-colnames(Q)<-c(0,1)
phy<-sim.history(phy,Q)

phy$tip.label


# tree
anole.tree<-read.tree("http://www.phytools.org/eqg2015/data/anole.tre")
## plot tree
plotTree(anole.tree,type="fan",ftype="i")

svl<-read.csv("http://www.phytools.org/eqg2015/data/svl.csv",
              row.names=1)
svl<-as.matrix(svl)[,1]

fit<-fastAnc(anole.tree,svl,vars=TRUE,CI=TRUE)


tree<-read.tree("tree.tre")
x<-as.matrix(read.csv("x.csv",row.names=1))[,1]

dotTree(tree,x,length=10,ftype="i")
plotTree.barplot(tree,x)



#### ACE
library(ape)

phy



plot(phy)
x <- factor(c(rep(0, 5), rep(1, 9)))
ans = ace(phy_states, phy, type = "d")
plot(phy, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(phy_states)], cex = 2, adj = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.75)




data(bird.orders)
bird.orders$tip.label
x <- factor(c(rep(0, 5), rep(1, 18)))
ans <- ace(x, bird.orders, type = "d")

#### Showing the likelihoods on each node:
plot(bird.orders, type = "c", FALSE, label.offset = 1)
co <- c("blue", "yellow")
tiplabels(pch = 22, bg = co[as.numeric(x)], cex = 2, adj = 1)
nodelabels(thermo = ans$lik.anc, piecol = co, cex = 0.75)

# continuous
x <- rnorm(23)
ace(x, bird.orders)



