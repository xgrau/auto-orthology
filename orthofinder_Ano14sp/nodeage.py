import ete3
import numpy as np
import pandas as pd
import markov_clustering
import networkx

# # define gene phylogeny
# nw = """
# ((Dme_001,Dme_002),(((Cfa_001,Mms_001),((((Hsa_001,Hsa_003),Ptr_001)
# ,Mmu_001),((Hsa_004,Ptr_004),Mmu_004))),(Ptr_002,(Hsa_002,Mmu_002))));
# """
# phy = ete3.PhyloTree(nw)
# phy.set_species_naming_function(lambda node: node.name.split("_")[0] )

# # define species phylogeny
# sp = "((((Hsa,Ptr),(Mmu,Mms)),Cfa),Dme);"
# #sp = "((((Hsa,Ptr),(Mmu,Mms)),(Cfa,(Lup,(Can,Cam)))),(Dme,Ano));"
# phs = ete3.PhyloTree(sp)
# phs.set_species_naming_function(lambda node: node.name )
# print(phs)


# gene tree
phy_fn = "/home/xavi/dades/Anotacions/orthofinder_Ano14sps_noclu_9oct19/output/Trees/OG0000000.iqt.treefile"
#phy    = ete3.PhyloTree("%s" % (phy_fn))
#phy.set_species_naming_function(lambda node: node.name.split("_")[0] )

# species tree
phs_fn = "/home/xavi/Documents/auto-orthology/orthofinder_Ano14sp/tree.newick"
phs    = ete3.PhyloTree("%s" % (phs_fn))
phs.set_species_naming_function(lambda node: node.name.split("_")[0] )

# define a dictionary of species-two-species relative ages
# for each species in the phyolgeny
def species_age_dict(phs):
		
	# init dict of dicts
	sps_age_dict = dict()
	sps_list     = phs.get_leaf_names()
	age_root     = 0
	for n,i in enumerate(sps_list):

		# init dict for species i
		sps_age_dict[i] = dict()
		for m,j in enumerate(sps_list):
			if n != m:
				e=phs.get_common_ancestor(i,j)
				d=e.get_farthest_leaf()[1]
			if n == m:
				d=0

			# new entry j in dict for species i
			sps_age_dict[i][j] = int(d)

			# is this the most ancient node?
			if sps_age_dict[i][j] > age_root:
				age_root = sps_age_dict[i][j]

	return sps_age_dict, sps_list, age_root # return dict of dicts

def parse_phylo(phy_fn):

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	phy.set_species_naming_function(lambda node: node.name.split("_")[0] )
	# resolve polytomies in a random fashion
	phy.resolve_polytomy(recursive=True)
	# check if tree is rooted, apply midpoint root if unrooted
	phy_root = phy.get_tree_root()
	phy_outg = phy_root.get_children()
	is_root  = len(phy_outg) == 2
	if is_root:
		pass
	else: 
		phy_outgroup = phy.get_midpoint_outgroup()
		phy.set_outgroup(phy_outgroup)

	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=0)

	# create empty array for network edges
	evo    = np.empty((len(evev)*1000, 5), dtype="object")
	evo[:] = np.nan
	# loop through in and out seqs, create edge table with orthologous events
	n = 0
	for ev in evev:
		if ev.etype == "S":
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					evo[n,0] = ii
					evo[n,1] = oi
					evo[n,2] = ev.branch_supports[0]
					evo[n,3] = ev.etype
					evo[n,4] = ev.sos
					n = n + 1

	evd = pd.DataFrame(evo).dropna()
	evd.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]

	return evd, phy


def clusters_mcl(evd, inf=1.1):

	# MCL clustering: create network
	evou_e = evd[["in_gene","out_gene","branch_support"]]
	evou_n = networkx.convert_matrix.from_pandas_edgelist(evou_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evou_n_nodelist = [ node for i, node in enumerate(evou_n.node()) ]
	evou_m = networkx.to_scipy_sparse_matrix(evou_n, nodelist=evou_n_nodelist)
	# MCL clustering: run clustering
	mcl_m  = markov_clustering.run_mcl(evou_m, inflation=inf)
	mcl_c  = markov_clustering.get_clusters(mcl_m)
	# MCL clustering: save output
	mcl_c_clu = [ i for i, cluster in enumerate(mcl_c) for node in cluster]
	mcl_c_noi = [ node for i, cluster in enumerate(mcl_c) for node in cluster]
	clu = pd.DataFrame( { 
		"node"    : [evou_n_nodelist[i] for i in mcl_c_noi],
		"cluster" : mcl_c_clu,
	}, columns=["node","cluster"])
	clu["cluster"] = clu["cluster"].astype(str)

	return clu


# obtain species-to-species dictionary of relative ages
sps_age_dict, sps_list, age_root = species_age_dict(phs=phs)

# load gene tree
evd,phy = parse_phylo(phy_fn=phy_fn)

# calculate ortholog clusters from gene tree
evc = clusters_mcl(evd=evd)



# age of divergence between all pairs of sequences
# not really useful without knowing relationships (evolutionary events) between sequences
# for n,i in enumerate(phy.get_leaves()):
# 	for m,j in enumerate(phy.get_leaves()):
# 		if i.name != j.name:
# 		# if i.species == "Hsa":
# 			e=phy.get_common_ancestor(i, j)
# 			a=e.get_age(sps_age_dict[i.species])
# 			print("Date of %s - %s divergence age: %d" % (i.name,j.name,a))
# 			#print(e)


# other stuff
clu_lis    = np.unique(evc.cluster)
clu_age = dict()
for c in clu_lis:
	clu_age[c] = 0

# ref_sp = "Anogam"

# clu_age_dic = dict()
# for s in sps_list:
# 	clu_age_dic[s] = dict()
# 	for c in clu_lis:
# 		clu_age_dic[s][c] = 0
# 	for n,i in enumerate(phy.get_leaves()):
# 		if i.species == s:
# 			for m,j in enumerate(phy.get_leaves()):
# 				if evc[evc["node"] == i.name]["cluster"].values == evc[evc["node"] == j.name]["cluster"].values:
# 					c=evc[evc["node"] == i.name]["cluster"].values[0]
# 					e=phy.get_common_ancestor(i, j)
# 					a=e.get_age(sps_age_dict[i.species])
					
# 					if a > clu_age_dic[s][c]:
# 						clu_age_dic[s][c] = a
				

clu_lis    = np.unique(evc.cluster)
clu_age_dic = dict()
for s in ["Anogam"]:
	
	# defin empty dict for ages relative to sps s
	clu_age_dic[s] = dict()
	for c in clu_lis:
		clu_age_dic[s][c] = 0

		# loop through nodes in cluster
		nod_lis = evc[evc["cluster"] == c]["node"].values
		for ni in nod_lis:
			i=phy.get_leaves_by_name(ni)[0]
			if i.species == s:
				for nj in nod_lis:
					j=phy.get_leaves_by_name(nj)[0]
					e=phy.get_common_ancestor(i, j)
					a=e.get_age(sps_age_dict[i.species])

					if a > clu_age_dic[s][c]:
						clu_age_dic[s][c] = a



phy_out = phy
for i in phy_out.get_leaves():
	c=evc[evc["node"] == i.name]["cluster"].values
	if c.size == 0: 
		c="NA"
		a="NA"
	else:
		c=c[0]
		a=clu_age_dic["Anogam"][c]
	i.name = str(i.name) + "_cluster_" + str(c) + "_age_" + str(a)



phy_out.write(outfile="ara.newick")
phy_out.show()
