from ete3 import PhyloTree
from ete3 import Tree

# read tree from file
tree = PhyloTree("adar.newick")

# find evolutionary events
evev = tree.get_descendant_evol_events()

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






