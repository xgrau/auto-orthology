# libraries
import sys
import numpy as np
import ete3
import pandas as pd
import logging
import os
import markov_clustering
import networkx
import random
import scipy as sp


# input variables
phy_fo = sys.argv[1]
out_fn = sys.argv[2]
#phy_fo = "output/Fasttrees/"
#out_fn = "ara"
phy_li = os.listdir(phy_fo)
phy_li = np.sort(phy_li)

# logging
logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s",
	#handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w"), logging.StreamHandler() ]
	handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w") ]
	)


# loop through phylogenies and recompute orthogropus with
# ETE + species overlap algorithm + MCL clustering
phi_n = 0
phi_l = len(phy_li)
out_d = pd.DataFrame()
print(phi_n,"/",phi_l)
for phi in phy_li: 

	phi_n = phi_n + 1
	if phi_l > 20:
		if phi_n % int(phi_l/20) == 0 : print(phi_n,"/",phi_l)

	# input name
	phy_fn = "%s/%s" % (phy_fo,phi)
	phy_id = phi.split(sep="/")[-1].split(sep=".")[0]
	logging.info("Phylogeny = %s" % phy_id)

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("Nodes = %i" % len(phy))

	# assign species names to tree
	phy.set_species_naming_function(lambda node: node.name.split("_")[0] )
	phy_sps     = [n.species for n in phy.get_leaves()]
	phy_sps_set = set(phy_sps)
	phy_seq     = [n.name for n in phy.get_leaves()]

	# resolve polytomies in a random fashion
	phy.resolve_polytomy(recursive=True)

	# check if tree is rooted, apply midpoint root if unrooted
	phy_root = phy.get_tree_root()
	phy_outg = phy_root.get_children()
	is_root  = len(phy_outg) == 2
	if is_root:
		pass
		logging.info("Tree is rooted, pass")
	else: 
		logging.info("Tree is unrooted, apply midpoint root")
		phy_outgroup = phy.get_midpoint_outgroup()
		phy.set_outgroup(phy_outgroup)

	# find evolutionary events (duplications and speciations)
	evev = phy.get_descendant_evol_events(sos_thr=0)

	# create empty array for network edges
	evou    = np.empty((len(evev)*1000, 5), dtype="object")
	evou[:] = np.nan
	# loop through in and out seqs, create edge table with orthologous events
	n = 0
	for ev in evev:
		if ev.etype == "S":
			for ii in ev.in_seqs:
				for oi in ev.out_seqs:
					evou[n,0] = ii
					evou[n,1] = oi
					evou[n,2] = ev.branch_supports[0]
					evou[n,3] = ev.etype
					evou[n,4] = ev.sos
					n = n + 1

	evou_d = pd.DataFrame(evou).dropna()
	evou_d.columns = ["in_gene","out_gene","branch_support","ev_type","sos"]

	# MCL clustering: create network
	evou_e = evou_d[["in_gene","out_gene","branch_support"]]
	evou_n = networkx.convert_matrix.from_pandas_edgelist(evou_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evou_n_nodelist = [ node for i, node in enumerate(evou_n.node()) ]
	evou_n_noderang = list(range(0,len(evou_n_nodelist)))
	evou_m = networkx.to_scipy_sparse_matrix(evou_n, nodelist=evou_n_nodelist)
	# MCL clustering: run clustering
	mcl_m  = markov_clustering.run_mcl(evou_m, inflation=2)
	mcl_c  = markov_clustering.get_clusters(mcl_m)
	#markov_clustering.draw_graph(mcl_m, mcl_c, node_size=50, with_labels=True, edge_color="k", cmap="Accent")
	# MCL clustering: save output
	mcl_c_clu = [ i for i, cluster in enumerate(mcl_c) for node in cluster]
	mcl_c_noi = [ node for i, cluster in enumerate(mcl_c) for node in cluster]
	mcl_c_out = pd.DataFrame( { 
		"node"    : [evou_n_nodelist[i] for i in mcl_c_noi],
		"id"      : mcl_c_noi,
		"cluster" : mcl_c_clu,
		"og"      : phy_id,
	}, columns=["node","id","cluster","og"])
	mcl_c_out["og_cluster"] = mcl_c_out["og"] +"_"+ mcl_c_out["cluster"].astype(str)

	# save clusters
	if phi_n == 1 : 
		mcl_c_out.to_csv("%s.a.csv" % out_fn, sep="\t", index=None, mode="w")
	if phi_n > 1  : 
		mcl_c_out.to_csv("%s.a.csv" % out_fn, sep="\t", index=None, mode="a", header=False)
	# # # save edges
	# evou_d.to_csv("tmp.edges", sep="\t", index=False, header=False)
	# logging.info("Species overlap score = 0, num speciation events = %i" % (evou_d.shape[0]))

	# # # run mcl
	# os.system("mcl tmp.edges --abc -o tmp.edges.mcl -I 1.5 2> /dev/null")
	# num_clu = sum(1 for line in open('tmp.edges.mcl'))
	# logging.info("Run MCL, num clusters = %i" % (num_clu))
	# os.system("awk '{ for (i=1; i < NF; i++) { print \"%s_\"NR\"\t\"$i  }}' tmp.edges.mcl >> tmp.edges.mcl.out 2> /dev/null" % (phy_id))

print(phi_n,"/",phi_l)

# save total output
out_d.to_csv("%s.a.csv" % out_fn, sep="\t", index=None)

# final housekeeping
# os.system("rm tmp.edges.mcl tmp.edges")
# os.system("mv tmp.edges.mcl.out %s.csv" % out_fn)

# end triumphantly
logging.info("All done in %s" % (phy_fo))


# evou_w = evou_e.pivot(index='in_gene', columns='out_gene', values='branch_support')
# evou_w = evou_w.fillna(0)
# evou_w = evou_w[0:]
# from sklearn.cluster import AgglomerativeClustering
# agc_c = AgglomerativeClustering(affinity="precomputed", linkage="single").fit(evou_w)
