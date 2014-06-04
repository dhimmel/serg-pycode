
import hetnet.readwrite.graph
import hetnet.permutation

graph_path = '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/140522-all-assoc-lessmsig/graph/graph.pkl.gz'
graph = hetnet.readwrite.graph.read_pickle(graph_path)

permuted_graph = hetnet.permutation.permute_graph2(graph, verbose=True)


permuted_graph_path = '/home/dhimmels/Documents/serg/gene-disease-hetnet/networks/140522-all-assoc-lessmsig-permuted/graph/graph.pkl.gz'
hetnet.readwrite.graph.write_pickle(permuted_graph, permuted_graph_path)












