# Automatic orthology

Scripts for automatic orthology assignment, proves.

## ETE3

Two strategies available [here](http://etetoolkit.org/docs/latest/tutorial/tutorial_phylogeny.html#detecting-evolutionary-events).

### Species overlap

See [Huerta-Cepas 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109). 

* Does not require a species tree.

* Only one tweakable parameter, **species overlap score** (`sos_thr`). By default, `sos_thr=0` which means that a single species in common between two node branches will rise a duplication event. This has been shown to perform the best with real data ([Huerta-Cepas 2007](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-6-r109)).

Script `ete-proves/test_v02_speciesoverlap.py` implements this method with a tree of ADAR enzymes (multigene, paralogy in animals, few outgroups). Observations:

* All identified events make sense (if it says something is )


```
D Anocul_ACUA003793-RA_2-103,Anomin_AMIN009728-RA_330-722 <====> Anocul_ACUA021577-RA_55-316
```

### Strict tree reconciliation

Too strict. It'll force the species tree onto every single subtree and add putative losses; which looks great in theory but requires a degree of compliance between gene trees and species trees that is unreasonable even for single-copy gene families. Won't for large multigene families.

## UPHO

That:
