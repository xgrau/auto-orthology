# libraries
from ete3 import PhyloTree
import sys

# input variables
# phy_fn = sys.argv[1]
phy_fn = "adar_hol.01.iqt.contree.newick"
ref_sp = "Hsap,Dromel"

# read tree from newick
phy = PhyloTree(phy_fn)

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

# print orthologs of each sequence in the list of evolutionary events
for ev in evev:
  for ii in ev.in_seqs:
    for oi in ev.out_seqs:
      print("%s\t%s\t%s" % (ii,oi,ev.etype))

