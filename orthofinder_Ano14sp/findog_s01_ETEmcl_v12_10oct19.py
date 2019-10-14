# libraries
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
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# argument parser
arp = argparse.ArgumentParser()

# Add the arguments to the parser
arp.add_argument("-p", "--phy",  required=True,  help="folder with phylogenies")
arp.add_argument("-s", "--suf",  required=True,  help="suffix of phylogenies in folder")
arp.add_argument("-o", "--out",  required=True,  help="prefix for output")
arp.add_argument("-t", "--ort",  required=True,  help="Orthology.txt table-like file produced by orthofinder")
arp.add_argument("-a", "--ani",  required=True,  help="which analysis to perform: \"main\" to analyse all genes, \"opti\" to find optimal inflation value")
arp.add_argument("-i", "--inf",  required=False, default=1.1, help="if analysis is \"main\", which inflation value to use? Default is 1.1")
arp.add_argument("-n", "--nopt", required=False, default=500, help="if analysis is \"opti\", how many phylogenies should we examine for optimisation? Default is 500")
arl = vars(arp.parse_args())

# input variables
phy_fo = arl["phy"]
phy_su = arl["suf"]
out_fn = arl["out"]
ort_fn = arl["ort"]
ani    = arl["ani"]
nopt   = int(arl["nopt"])
inf    = int(arl["inf"])

# logging
logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s",
	#handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w"), logging.StreamHandler() ]
	handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w") ]
	)

# prepare orthology from orthofinder
os.system("awk '{ for (i=2; i <= NF; i++) { print $1\"\t\"$i  }}' %s | sed \"s/://\" > %s.orthology_preete.txt" % (ort_fn,out_fn))
ort = pd.read_csv("%s.orthology_preete.txt" % out_fn, sep="\t", header=None)
ort.columns = ["og","node"]

# list of orthogroups
ort_lis = np.unique(ort["og"])



### FUNCTIONS ###

# loop through some orthogroups to find optimal inflation
def optimisation_loop(nopt=nopt):
	
	# loop through phylogenies
	phi_n = 0
	print(phi_n,"/",nopt)
	inf_lis = np.zeros(nopt)
	mod_lis = np.zeros(nopt)

	ort_ran = np.random.choice(ort_lis, size=nopt)
	for phi in ort_ran:

		# input name
		phy_fn = "%s/%s.%s" % (phy_fo,phi,phy_su)
		phy_id = phi.split(sep="/")[-1].split(sep=".")[0]

		# cluster phylogeny if you can, retrieve original clusters if you can't
		if os.path.exists(phy_fn):
			evou_d = parse_phylo(phy_fn=phy_fn, phy_id=phy_id)
			inf_lis[phi_n], mod_lis[phi_n] = clusters_opt(phy_fn=phy_fn, phy_id=phy_id, evou_d=evou_d)
		else:
			inf_lis[phi_n] = np.nan 
			mod_lis[phi_n] = np.nan

		# add counter
		phi_n = phi_n + 1
		if phi_n % int(nopt/20) == 0 : print(phi_n,"/",nopt)

	# calculate means
	print(phi_n,"/",nopt)
	print("Mean inflation = %f | Median inflation = %f " % (np.nanmean(inf_lis), np.nanmedian(inf_lis)))
	print("Mean modularity = %f | Median modularity = %f " % (np.nanmean(mod_lis), np.nanmedian(mod_lis)))
	
	return inf_lis, mod_lis

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


# loop through all orthogroups and recluster them using 
# ETE + species overlap algorithm + MCL clustering
def orthogroup_loop():

	# loop through phylogenies
	phi_n = 0
	phi_l = len(ort_lis)
	print(phi_n,"/",phi_l)

	for phi in ort_lis:

		# input name
		phy_fn = "%s/%s.%s" % (phy_fo,phi,phy_su)
		phy_id = phi.split(sep="/")[-1].split(sep=".")[0]

		# cluster phylogeny if you can, retrieve original clusters if you can't
		if os.path.exists(phy_fn):
			evou_d = parse_phylo(phy_fn=phy_fn, phy_id=phy_id)
			if len(evou_d) > 1:
				mcl_c_out = clusters_mcl(phy_fn=phy_fn, phy_id=phy_id, evou_d=evou_d)
			else:
				mcl_c_out = clusters_nophylo(phy_fn=phy_fn, phy_id=phy_id)
		else:
			mcl_c_out = clusters_nophylo(phy_fn=phy_fn, phy_id=phy_id)

		# create dataframe with old and new clusters, and all genes
		out_d = ort[ort["og"] == phy_id]
		out_d = pd.merge(out_d, mcl_c_out, how="outer", on="node")
		out_d["og_cluster"] = out_d["og"] +"_"+ out_d["cluster"].astype(str)
		
		# save clusters
		if phi_n == 0 : 
			out_d.to_csv("%s.orthology.csv" % out_fn, sep="\t", index=None, mode="w")
		if phi_n > 0  : 
			out_d.to_csv("%s.orthology.csv" % out_fn, sep="\t", index=None, mode="a", header=False)
	
		# add counter
		phi_n = phi_n + 1
		if phi_l > 20:
			if phi_n % int(phi_l/20) == 0 : print(phi_n,"/",phi_l)

	# end triumphantly
	print(phi_n,"/",phi_l)


# parse phylogenies with ETE to obtain a network-like table defining 
# orthologous relationships, using the species overlap algorithm
def parse_phylo(phy_fn, phy_id):

	# load input
	phy = ete3.PhyloTree("%s" % (phy_fn))
	logging.info("%s num nodes = %i" % (phy_id,len(phy)))
	# assign species names to tree
	phy.set_species_naming_function(lambda node: node.name.split("_")[0] )
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

	return evou_d

# function to calculate optimal inflation for MCL for each phylogeny;
# i.e., that with the highest network modularity (Q)
def clusters_opt(phy_fn, phy_id, evou_d):
	
	# MCL clustering: create network
	logging.info("%s Create network" % phy_id)
	evou_e = evou_d[["in_gene","out_gene","branch_support"]]
	evou_n = networkx.convert_matrix.from_pandas_edgelist(evou_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evou_n_nodelist = [ node for i, node in enumerate(evou_n.node()) ]
	evou_m = networkx.to_scipy_sparse_matrix(evou_n, nodelist=evou_n_nodelist)
	# MCL clustering: find optimal inflation value
	inf,mod = optimise_inflation(matrix=evou_m)
	
	return inf, mod

# function to cluster a network-like table of orthologs (from ETE) using MCL
def clusters_mcl(phy_fn, phy_id, evou_d):

	# MCL clustering: create network
	logging.info("%s Create network" % phy_id)
	evou_e = evou_d[["in_gene","out_gene","branch_support"]]
	evou_n = networkx.convert_matrix.from_pandas_edgelist(evou_e, source="in_gene", target="out_gene", edge_attr="branch_support")
	evou_n_nodelist = [ node for i, node in enumerate(evou_n.node()) ]
	evou_m = networkx.to_scipy_sparse_matrix(evou_n, nodelist=evou_n_nodelist)
	# MCL clustering: run clustering
	# inf,_ = optimise_inflation(matrix=evou_m)
	logging.info("%s MCL clustering, inflation = %f" % (phy_id, inf))
	mcl_m  = markov_clustering.run_mcl(evou_m, inflation=inf)
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

	return mcl_c_out

# if the phylogeny can't be analysed with ETE (no phylogeny, not enough speciation events...), use
# the original orthogroups from orthofinder instead
def clusters_nophylo(phy_fn, phy_id):

	logging.info("%s Can't find phylogeny (small orthogroup?) OR can't perform MCL clusters (no speciations?), output original clusters instead" % phy_id)
	mcl_c_out = pd.DataFrame( { 
		"node"    : ort[ort["og"] == phy_id]["node"].values,
		"cluster" : np.nan
	}, columns=["node","cluster"])
	logging.info("%s No phylo/MCL, num clusters = %i" % (phy_id, 1))
	logging.info("%s No phylo/MCL, num clustered genes = %i" % (phy_id, len(mcl_c_out)))

	return mcl_c_out


### MAIN ####

# log input variables
logging.info("Input args: %r", arl)

# main analysis
if ani == "main":

	# run orthologs loop
	logging.info("Analyse orthogroups in %s" % (phy_fo))
	orthogroup_loop()

elif ani == "opti":

	# run optimisation loop
	logging.info("Find optimal inflation")
	inf_lis, mod_lis = optimisation_loop(nopt=nopt)

	# plot optimisation histograms
	with PdfPages('%s.optimise_inflation.pdf' % out_fn) as pdf:
		# inflation
		plt.figure(figsize=(4,3))
		plt.title("Hist inflation I\nmedian = %f" % np.nanmedian(inf_lis))
		plt.hist(inf_lis[~np.isnan(inf_lis)])
		pdf.savefig(bbox_inches='tight')
		# modularity
		plt.figure(figsize=(4,3))
		plt.title("Hist modularity Q\nmedian = %f" % np.nanmedian(mod_lis))
		plt.hist(mod_lis[~np.isnan(mod_lis)])
		pdf.savefig(bbox_inches='tight')
		# close pdf
		plt.close()

else: 
	print("ERROR: specify analysis with -a/--ani")

# end
logging.info("All done in %s" % (phy_fo))

