# libraries
import sys
import numpy as np
import ete3


# input variables
phy_fn = "set.newick"
ref_sp = "Hsap"
ref_sl = ref_sp.split(sep=",")
# phy_fn = sys.argv[1]
# ref_sp = sys.argv[2]


# read tree from newick
phy = ete3.PhyloTree(phy_fn)


# assign species names to tree
phy.set_species_naming_function(lambda node: node.name.split("_")[0] )
out_sps = [n.species for n in phy.get_leaves()]
out_seq = [n.name for n in phy.get_leaves()]


# check if tree is rooted, apply midpoint root if unrooted
phy_root = phy.get_tree_root()
phy_outg = phy_root.get_children()
is_root  = len(phy_outg) == 2
if is_root:
  print("Tree is rooted, pass")
else: 
  print("Tree is unrooted, apply midpoint root")
  phy_outgroup = phy.get_midpoint_outgroup()
  phy.set_outgroup(phy_outgroup)


# find evolutionary events (duplications and speciations)
evev = phy.get_descendant_evol_events(sos_thr=0.0)


fseqs = lambda slist: [s for s in slist if s.startswith(ref_sl)]
# print orthologs of each sequence in the list of evolutionary events
for ev in evev:
  if ev.etype == "S":
    for ii in ev.in_seqs:
      for oi in ev.out_seqs:
        print("%s\t%s\t%s" % (ii,oi,ev.etype))

