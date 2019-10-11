import markov_clustering
import networkx
import random


numnodes = 100
positions = {i:(random.random() * 2 - 1, random.random() * 2 - 1) for i in range(numnodes)}
network = networkx.random_geometric_graph(numnodes, 0.3, pos=positions)
matrix = networkx.to_scipy_sparse_matrix(network)

result = markov_clustering.run_mcl(matrix, inflation=2)
clusters = markov_clustering.get_clusters(result)    # get clusters
markov_clustering.draw_graph(matrix, clusters, pos=positions, node_size=50, with_labels=False, edge_color="k", cmap="magma")


for inflation in [i / 10 for i in range(15, 26)]:
    result = markov_clustering.run_mcl(matrix, inflation=inflation)
    clusters = markov_clustering.get_clusters(result)
    Q = markov_clustering.modularity(matrix=result, clusters=clusters)
    print("inflation:", inflation, "modularity:", Q)


from sklearn.cluster import AgglomerativeClustering
import numpy as np
X = np.array([[1, 2], [1, 4], [1, 0],
              [4, 2], [4, 4], [4, 0]])
clustering = AgglomerativeClustering().fit(X)
clustering 




clustering.labels_
