from ete3 import PhyloTree
from ete3 import Tree

# read tree from file
phy = PhyloTree("adar.newick")


# assign species names to tree
phy.set_species_naming_function(lambda node: node.name.split("_")[0] )
for n in phy.get_leaves():
  print("node:", n.name, "Species name:", n.species)


# find evolutionary events
evev = phy.get_descendant_evol_events(sos_thr=0.9)

for ev in evev:
  if ev.etype == "S":
    print(ev.orthologs)

# all events
for ev in evev:
  print(ev.etype, ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))



# all events involving either Hsap or Drer
fseqs = lambda slist: [s for s in slist if s.startswith("Drer") or s.startswith("Hsap")]
for ev in evev:
  if ev.etype == "D":
    print('Paralog: ', ','.join(fseqs(ev.in_seqs)), "<====>", ','.join(fseqs(ev.out_seqs)))

for ev in evev:
  if ev.etype == "S":
    print('Ortholog:', ','.join(fseqs(ev.in_seqs)), "<====>", ','.join(fseqs(ev.out_seqs)))



# obtain duplication events
ntrees, ndups, sptrees =  phy.get_speciation_trees()
print("Found %d species trees and %d duplication nodes" %(ntrees, ndups))
for spt in sptrees:
   print(spt)

