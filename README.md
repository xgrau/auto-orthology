# Automatic orthology

Scripts for automatic orthology assignment, proves.

## Pipeline

Steps:

0. Create phylogeny somehow.

1. **`findog_s01_ETEnet_v03_27set19.py`** (Python 3.5.5): takes as input a newick phylogeny and identifies speciation events using `ete`. Then, it creates a network-like table of phylogeny nodes (leaves) that appeared via the same speciation event (using various species overlap score thresholds). Each table entry is a pair of `in_seqs` and `out_seqs`, like this:
```
seq1	seq2	S	100.0	0.01
seq1	seq3	S	100.0	0.01
...
```
```
findog_s01_ETEnet_v03_27set19.py <input newick> <output prefix>
```

2. **`findog_s02_tabulate_v01.R`** (R 3.6.1): Assign each component of this network of orthologs to one or more orthogroups. Requires `igraph`.

Required inputs:

* `set_raxml.newick`: phylogeny, newick, includes supports.
* `set_Hsap_names.dict`: dictionary linking specific sequences from one specie (or more?) to gene names (orthogroups that'll be tabulated)
* `TODO`: list of species to create report (so as to include also species that are not present in the phylogeny).

Final outputs (maybe):

* table of ortholog presence/absence per species, or table of counts (including inparalogs?)
* table of ortholog support per species (branch support of the paralog with the highest support per species?)
* lists of genes belonging to each cluster

## ETE3

Two strategies available [here](http://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html#detecting-evolutionary-events).

### Species overlap

See [Huerta-Cepas 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109). 

* Does not require a species tree.

* Only one tweakable parameter, **species overlap score** (`sos_thr`). By default, `sos_thr=0` which means that a single species in common between two node branches will rise a duplication event. This has been shown to perform the best with real data ([Huerta-Cepas 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)).

Script `ete-proves/test_v02_speciesoverlap.py` implements this method with a tree of ADAR enzymes (multigene, paralogy in animals, few outgroups). Observations:

* All identified events make sense (paralogs, orthologs) for ADAR.
* False negatives?


```
D Anocul_ACUA003793-RA_2-103,Anomin_AMIN009728-RA_330-722 <====> Anocul_ACUA021577-RA_55-316
```

### Strict tree reconciliation

Too strict. It'll force the species tree onto every possible subtree and add putative losses; which looks great in theory but requires a degree of compliance between gene trees and species trees that is unreasonable even for single-copy gene families. Won't for large multigene families.

## UPHO

### Setup

1. Download UPHO:

```
git clone https://github.com/ballesterus/UPhO.git
```

2. Run UPHO:

```
python2 /home/xavi/Programes/UPhO/UPhO.py -d _ -in tree/adar.newick -iP -ouT
```

Parameters:

* `-iP` to include inparalogs in the orthplogy groups
* `-m` to specify the minimum number of OTUs in an orthogroup
* `-ouT` Write orthologous branches to newick file
* `-S` minimum support value for orthology evaluation
* `-d` delimiter, default is `|`, but if you use `_` it seems to work fine (even if there are more than one such characters per name)

Problem: highly atomised for some reason! There's no way to tweak granularity.
