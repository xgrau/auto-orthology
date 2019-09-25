# In[]:
import ete3
from ete3 import PhyloTree

# In[]:
# read tree from file
tphy = PhyloTree("/home/xavi/Documents/scripts/ete-proves/cyps.newick")
tsps = PhyloTree("/home/xavi/Documents/scripts/ete-proves/cyps_sps_22mosquits.newick")

# In[]:
# assign species names to tree
def get_species_name(node_name_string):
  # Species code is the first part of leaf name (separated by an
  #  underscore character)
  spcode = node_name_string.split("_")[0]
  return spcode

tphy.set_species_naming_function(get_species_name)
for n in tphy.get_leaves():
  print("node:", n.name, "Species name:", n.species)


# In[]:
# find evolutionary events using tree reconciliation
tree_rec, evev_rec = tphy.reconcile(tsps)

# In[]:
print(tree_rec)
tree_rec.show()

# In[]:
# find evolutionary events using species overlap
evev = tphy.get_descendant_evol_events()

for ev in evev:
  if ev.etype == "S":
    for s in ev.in_seqs:
      if s.startswith("Lepsal"):
        print(ev.orthologs)

fseqs = lambda slist: [s for s in slist if s.startswith("Drer") or s.startswith("Hsap")]
print("\nOrthology relationships among Anogam and Anosin")
for ev in evev:
  if ev.etype == "D":
    print('Paralog: ', ','.join(fseqs(ev.in_seqs)), "<====>", ','.join(fseqs(ev.out_seqs)))

for ev in evev:
  if ev.etype == "S":
    print('Ortholog:', ','.join(fseqs(ev.in_seqs)), "<====>", ','.join(fseqs(ev.out_seqs)))






