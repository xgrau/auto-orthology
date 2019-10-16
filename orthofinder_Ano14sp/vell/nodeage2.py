import ete3
import numpy as np
import pandas as pd

# input files
phs_fn   = "/home/xavi/Documents/auto-orthology/orthofinder_Ano14sp/tree.newick"
ort_fn   = "/home/xavi/Documents/auto-orthology/orthofinder_Ano14sp/outiqt.orthology.csv"
gene_col = "node"
clus_col = "og_cluster"


# define a dictionary of species-two-species relative ages
# for each species in the phyolgeny
def species_age_dict(phs):
		
	# init dict of dicts
	sps_age_dict = dict()
	sps_out_dict = dict()
	sps_list     = phs.get_leaf_names()
	age_root     = 0
	for n,i in enumerate(sps_list):

		# init dict of ages for species i
		sps_age_dict[i] = dict()
		for m,j in enumerate(sps_list):
			if n != m:
				e = phs.get_common_ancestor(i,j) # define species subtree between species i and j
				a = e.get_farthest_leaf()[1]     # distance to the farthest leaf in subtree e (age)
			if n == m:
				a = 0

			# new entry j in dict for species i
			sps_age_dict[i][j] = int(a)

			# is this the most ancient node?
			if sps_age_dict[i][j] > age_root:
				age_root = sps_age_dict[i][j]
		
		# init dict of outgroups for species i
		sps_out_dict[i] = dict()
		for m,j in enumerate(sps_list):
			if n != m:
				e = phs.get_common_ancestor(i,j) # define species subtree between species i and j
				l = e.get_leaf_names()           # which species are present in subtree e?
				d = dict()                       # modified dict of inter-species ages: 
				for k in sps_list:               #
					if np.isin(k, l):             # if species in subtree, use real age
						d[k] = sps_age_dict[i][k]  #
					else:                         # if species not in subtree, age = 0 
						d[k] = 0
				o = e.get_farthest_oldest_leaf(d).name # find farthest oldest leaf by name (outgroup to sps i)
				sps_out_dict[i][j] = i+","+o     # store age as pair that defines a last common ancestor
			if n == m:
				sps_out_dict[i][j] = i+","+i      

	return sps_age_dict, sps_out_dict, sps_list, age_root




# load input tree
phs    = ete3.PhyloTree("%s" % (phs_fn))
phs.set_species_naming_function(lambda node: node.name.split("_")[0] )

# obtain species-to-species dictionary of relative ages
sps_age_dict, sps_out_dict, sps_list, age_root = species_age_dict(phs=phs)

# load orthoclusters
ort = pd.read_csv(ort_fn, sep="\t")
ort = pd.DataFrame(ort).dropna()
ort = ort[[gene_col,clus_col]]

# lister of orthoclusters
clus_lis = np.unique(ort[clus_col])
ages_lis = np.zeros(len(clus_lis))
outg_lis = np.zeros(len(clus_lis), dtype=object)

# loop orthoclusters
for r in ["Anogam"]:
	for n,c in enumerate(clus_lis):
		nod_clu = ort[ort[clus_col] == c][gene_col].values
		sps_clu = np.unique([ m.split("_")[0] for m in nod_clu ])

		a = 0
		for s in sps_clu:
			t=sps_age_dict[r][s]
			
			if t > a:
				a = t
				o = s

		ages_lis[n] = int(a)
		outg_lis[n] = sps_out_dict[r][o]

dat = pd.DataFrame({
	"orthogroup"   : clus_lis,
	"relative_age" : ages_lis,
	"outgroup"     : outg_lis
})