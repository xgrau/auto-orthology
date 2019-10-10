# libraries
import sys
import numpy as np
import ete3
import pandas as pd
import logging
import os

# input variables
phy_fn = "output/Fasttrees/OG0000000.tree"
#phy_fn = sys.argv[1]
out_fn = "%s.out_ete" % phy_fn.split(sep=".")[0]

# logging
logging.basicConfig(
	level=logging.DEBUG, 
	format="%(asctime)s [%(levelname)-5.5s]\t%(message)s",
	handlers=[ logging.FileHandler("%s.log" % out_fn, mode="w"), logging.StreamHandler() ]
	)


# load input
phy = ete3.PhyloTree("%s" % (phy_fn))
logging.info("Phylogeny = %s" % phy_fn)
logging.info("Nodes = %i" % len(phy))


# assign species names to tree
phy.set_species_naming_function(lambda node: node.name.split("_")[0] )
phy_sps     = [n.species for n in phy.get_leaves()]
phy_sps_set = set(phy_sps)
phy_seq     = [n.name for n in phy.get_leaves()]

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

# save edges
evou_d.to_csv("%s.edges" % (out_fn), sep="\t", index=False, header=False)
logging.info("Species overlap score = 0, num speciation events = %i" % (evou_d.shape[0]))

# run mcl
os.system("mcl %s.edges --abc -o %s.edges.mcl -I 1.5 2> /dev/null" % (out_fn,out_fn))

# save nodes
noou_d = pd.DataFrame({
	"gene" : phy_seq,
	"species" : phy_sps
})

noou_d.to_csv("%s.nodes" % (out_fn), sep="\t", index=False)
logging.info("Num nodes = %i" % noou_d.shape[0])
logging.info("Num species = %i" % len(np.unique(noou_d["species"])))

# acaba
logging.info("Done!")

