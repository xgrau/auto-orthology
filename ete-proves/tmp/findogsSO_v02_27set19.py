# libraries
import sys
import numpy as np
import ete3
import pandas as pd


# input variables
phy_fn = "set_raxml.newick"
ref_sp = "Hsap"
ref_sl = ref_sp.split(sep=",")
# phy_fn = sys.argv[1]
# ref_sp = sys.argv[2]
man_fn = "set_manualclass.ml.csv"
dic_fn = "set_Hsap_names.dict"


# load input
phy = ete3.PhyloTree(phy_fn)
man = pd.read_table(man_fn)
dic = pd.read_table(dic_fn)
print("Phylogeny = %s\nNodes = %i" % (phy_fn, len(phy)))



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
  print("Tree is rooted, pass")
else: 
  print("Tree is unrooted, apply midpoint root")
  phy_outgroup = phy.get_midpoint_outgroup()
  phy.set_outgroup(phy_outgroup)


# find evolutionary events (duplications and speciations)
evev = phy.get_descendant_evol_events(sos_thr=0.0)

# retrieve lists of orthologous
evev_inseqs = [None] * len(evev)
evev_ouseqs = [None] * len(evev)
evev_genes  = [None] * len(evev)
evev_sps    = [None] * len(evev)
evev_brsup  = [None] * len(evev)
evev_annog  = [None] * len(evev)
evev_annot  = [None] * len(evev)
evev_type   = [None] * len(evev)
for n,ev in enumerate(evev):
	evev_inseqs[n] = ev.in_seqs
	evev_ouseqs[n] = ev.out_seqs
	evev_genes[n]  = ev.orthologs.union(ev.inparalogs)
	evev_sps[n]    = set([i.split("_")[0] for i in evev_genes[n]])
	evev_type[n]   = ev.etype
	evev_brsup[n]  = ev.branch_supports[0]
	evev_annog[n]  = list(filter(None.__ne__, [i if i.split("_")[0] in ref_sp else None for i in evev_genes[n]]))
	evev_annot[n]  = evev_annog[n]


for n,ev in enumerate(evev):
	if evev_brsup[n] >= 25 and evev_type[n] == "S":
		phy_seq[n]

evev_genes[2]
phy_seq[2]

for i in range(100):
	evev[i].branch_supports
# loop through seqs finding human orthologs
for n,si in enumerate(phy_seq):
	print(si)






# print orthologs of each sequence in the list of evolutionary events
for ev in evev:
	if ev.etype == "S":
		for oi in ev.out_seqs:
			for ii in ev.in_seqs:
				if ii.startswith(ref_sp):
					ri = ii
					fi = dic[dic["gene"] == ri]["family"].values
					if len(fi) > 0 :
						fi = fi[0]
					else:
						fi = np.nan
					print("%s\t%s\t%s\t%s" % (ii,oi,ev.etype,fi))

				elif oi.startswith(ref_sp):
					ri = oi
					fi = dic[dic["gene"] == ri]["family"].values
					if len(fi) > 0 :
						fi = fi[0]
					else:
						fi = np.nan
					print("%s\t%s\t%s\t%s" % (oi,ii,ev.etype,fi))









# for ev in evev:
#   if ev.etype == "S":
#     print(ev.orthologs)

