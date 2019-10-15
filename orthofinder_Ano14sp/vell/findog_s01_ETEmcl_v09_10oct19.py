# libraries
import sys
import os
import fnmatch
import numpy as np
import scipy as sp
import pandas as pd
import ete3
import markov_clustering
import logging
import networkx
import argparse

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-p", "--phy", required=True, help="folder with phylogenies")
arp.add_argument("-s", "--suf", required=True, help="suffix of phylogenies in folder")
arp.add_argument("-o", "--out", required=True, help="prefix for output")
arp.add_argument("-t", "--ogt", required=True, help="Orthology.txt table-like file produced by orthofinder")
arl = vars(arp.parse_args())

# input variables
# phy_fo = sys.argv[1] # folder with phylogenies
# phy_su = sys.argv[2] # suffix of newick trees (e.g. "newick" if "XXX.newick")
# out_fn = sys.argv[3] # output prefix
# ort_fn = sys.argv[4] # path to Orthology.txt (to retrieve orthogroups that couldn't be analysed with phylogenies)
phy_fo = arl["phy"]
phy_su = arl["suf"]
out_fn = arl["out"]
ort_fn = arl["ogt"]

# logging
logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s",
	#handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w"), logging.StreamHandler() ]
	handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w") ]
	)

# log input variables
logging.info("Input args: %r", arl)

# list of trees
phy_li = fnmatch.filter(os.listdir(phy_fo), '*.%s' % phy_su)
phy_li = np.sort(phy_li)

# inflation for MCL clustering
inf = 1.5


# find optimal inflation for MCL
def optimise_inflation(matrix, start=1.1, end=2.5, step=0.1):
	I_lis = np.arange(start, end, step).tolist()
	Q_lis = np.zeros(shape=(len(I_lis),1))
	for n,I in enumerate(I_lis):
		result   = markov_clustering.run_mcl(matrix, inflation=I)
		clusters = markov_clustering.get_clusters(result)
		Q        = markov_clustering.modularity(matrix=result, clusters=clusters)
		Q_lis[n] = Q
	max_Q_index = np.argmax(Q_lis)
	# return inflation with maximum modularity and modularities array
	return I_lis[max_Q_index], Q_lis[max_Q_index]

# loop through phylogenies and recompute orthogropus with
# ETE + species overlap algorithm + MCL clustering
phi_n = 0
phi_l = len(phy_li)
print(phi_n,"/",phi_l)
for phi in phy_li: 

	phi_n = phi_n + 1
	if phi_l > 20:
		if phi_n % int(phi_l/20) == 0 : print(phi_n,"/",phi_l)

	# input name
	phy_fn = "%s/%s" % (phy_fo,phi)
	phy_id = phi.split(sep="/")[-1].split(sep=".")[0]

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("%s num nodes = %i" % (phy_id,len(phy)))

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
		logging.info("%s Tree is rooted, pass" % phy_id)
	else: 
		logging.info("%s Tree is unrooted, apply midpoint root" % phy_id)
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

	if len(evou_d) > 1:

		logging.info("%s Create network" % phy_id)
		# MCL clustering: create network
		evou_e = evou_d[["in_gene","out_gene","branch_support"]]
		evou_n = networkx.convert_matrix.from_pandas_edgelist(evou_e, source="in_gene", target="out_gene", edge_attr="branch_support")
		evou_n_nodelist = [ node for i, node in enumerate(evou_n.node()) ]
		evou_n_noderang = list(range(0,len(evou_n_nodelist)))
		evou_m = networkx.to_scipy_sparse_matrix(evou_n, nodelist=evou_n_nodelist)
		# MCL clustering: run clustering
		# inf,_ = optimise_inflation(matrix=evou_m)
		logging.info("%s MCL clustering, optimal inflation = %f (max Q modularity)" % (phy_id, inf))
		mcl_m  = markov_clustering.run_mcl(evou_m, inflation=9)
		mcl_c  = markov_clustering.get_clusters(mcl_m)
		logging.info("%s MCL clustering, num clusters = %i" % (phy_id, len(mcl_c)))
		# markov_clustering.draw_graph(mcl_m, mcl_c, node_size=50, with_labels=True, edge_color="k", cmap="Accent")
		# MCL clustering: save output
		mcl_c_clu = [ i for i, cluster in enumerate(mcl_c) for node in cluster]
		mcl_c_noi = [ node for i, cluster in enumerate(mcl_c) for node in cluster]
		mcl_c_out = pd.DataFrame( { 
			"node"    : [evou_n_nodelist[i] for i in mcl_c_noi],
			"cluster" : mcl_c_clu,
		}, columns=["node","cluster"])
		mcl_c_out["cluster"] = mcl_c_out["cluster"].astype(str)
		logging.info("%s MCL clustering, num clustered genes = %i" % (phy_id, len(mcl_c_out)))

	else:

		logging.info("%s Can't create network (there are no speciation events? all genes from the same sps?), output original clusters instead" % phy_id)
		mcl_c_out = pd.DataFrame( { 
			"node"    : phy_seq,
			"cluster" : np.nan,
		}, columns=["node","cluster"])
		logging.info("%s No MCL clustering, num clusters = %i" % (phy_id, 1))
		logging.info("%s No MCL clustering, num clustered genes = %i" % (phy_id, len(mcl_c_out)))
		
	# create dataframe with old and new clusters, and all genes
	out_d = pd.DataFrame( { 
		"node"    : phy_seq,
		"og"      : phy_id
	}, columns=["node","og"])
	out_d = pd.merge(out_d, mcl_c_out, how="outer", on="node")
	out_d["og_cluster"] = out_d["og"] +"_"+ out_d["cluster"].astype(str)
	# save clusters
	if phi_n == 1 : 
		out_d.to_csv("%s.csv" % out_fn, sep="\t", index=None, mode="w")
	if phi_n > 1  : 
		out_d.to_csv("%s.csv" % out_fn, sep="\t", index=None, mode="a", header=False)

print(phi_n,"/",phi_l)


# end triumphantly
logging.info("All done in %s" % (phy_fo))


