import itertools
import os
import csv
import collections

import networkx
import numpy
import numpy.linalg

import bioparser.data




class DiseasomeMaker(object):

    def __init__(self):
        self.ontology = self.get_ontology()

    def set_gene_annotations(self, node_goto_node=dict()):
        node_to_genes = dict()

        gene_dicts = [self.get_gwas_catalog_annotations(),
                      self.get_omim_annotations(),
                      #self.get_ctd_annotations()
                      ]

        nodes = reduce(lambda x, y: set.union(set(x), set(y)), gene_dicts)
        node_to_genes = dict()
        for node in nodes:
            for gene_dict in gene_dicts:
                genes = gene_dict.get(node, set())
                node_to_genes.setdefault(node, set()).update(genes)

        for source_node, sink_node in node_goto_node.items():
            genes = node_to_genes.pop(source_node, set())
            node_to_genes.setdefault(sink_node, set()).update(genes)
        # More efficient, does not enable printing.
        #for node_to_genes_part in node_to_genes_parts:
        #    for key, value in node_to_genes_part.iteritems():
        #        node_to_genes.setdefault(key, set()).update(value)
        self.node_to_genes = node_to_genes
        return node_to_genes

    def get_gwas_catalog_annotations(self):
        # Load gwas_catalog and gene annotations.
        gwas_catalog = bioparser.data.Data().gwas_catalog
        node_to_genes = gwas_catalog.get_doid_id_to_genes(p_cutoff=None, fdr_cutoff=None, mapped_term_cutoff=1, exclude_pmids=set())
        return node_to_genes

    def get_omim_annotations(self):

        morbid_map = bioparser.data.Data().morbid_map
        omim_associations = morbid_map.get_associations()
        mim_to_genes = dict()
        for omim_association in omim_associations:
            mim = omim_association['mim_number']
            gene = omim_association['gene']
            mim_to_genes.setdefault(mim, set()).add(gene)

        doid = bioparser.data.Data().doid
        doid_to_xrefs = doid.get_doid_to_xrefs('OMIM')
        node_to_genes = dict()
        for doid_code, mims in doid_to_xrefs.iteritems():
            genes = set()
            for mim in mims:
                genes |= mim_to_genes.get(mim, set())
            if not genes:
                continue
            node_to_genes[doid_code] = genes

        return node_to_genes

    def get_ctd_annotations(self):

        ctd = bioparser.data.Data().ctd
        symbol_to_gene = bioparser.data.Data().hgnc.get_symbol_to_gene()
        medic_to_genes = dict()
        for ctd_row in ctd.read_gene2disease_filtered():
            if 'marker/mechanism' not in ctd_row['DirectEvidence']:
                continue
            symbol = ctd_row['GeneSymbol']
            gene = symbol_to_gene.get(symbol)
            if not gene:
                continue
            medic_to_genes.setdefault(ctd_row['DiseaseID'], set()).add(gene)

        doid = bioparser.data.Data().doid
        doid_to_xrefs = doid.get_doid_to_xrefs('MSH', 'MESH:')
        doid_to_xrefs.update(doid.get_doid_to_xrefs('OMIM', 'OMIM:'))

        node_to_genes = dict()
        for doid_code, medic_codes in doid_to_xrefs.iteritems():
            genes = set()
            for medic_code in medic_codes:
                genes |= medic_to_genes.get(medic_code, set())
            if not genes:
                continue
            node_to_genes[doid_code] = genes

        return node_to_genes


    def get_ontology(self):
        # Load the disease ontology.
        doid = bioparser.data.Data().doid
        ontology = doid.get_ontology()
        return ontology

    def get_graph(self, gene_minimum=10, exclude=set(), include=set()):
        """
        Returns nodes only
        gene_minimum - nodes with fewer than gene_minimum annotated genes are excluded
        exclude - set of nodes to exclude from the analysis.
        """
        graph = networkx.Graph()
        keep_data_keys = ['name']
        assert not exclude & include
        for node, genes in self.node_to_genes.items():
            if node in exclude:
                continue
            if node not in include and len(genes) < gene_minimum:
                continue
            data = self.ontology.graph.node[node]
            data = {key: data[key] for key in keep_data_keys}
            data['genes'] = genes
            graph.add_node(node, data)

        return graph

    @staticmethod
    def get_gene_counter(graph):
        gene_counter = collections.Counter()
        for node, data in graph.nodes_iter(True):
            gene_counter.update(data['genes'])
        return gene_counter

    @staticmethod
    def get_overlap(set_0, set_1, elem_counter=None):
        assert isinstance(set_0, set)
        assert isinstance(set_1, set)
        intersect = set_0 & set_1
        union = set_0 | set_1
        if elem_counter:
            intersect = sum(elem_counter[elem] ** 0.5 for elem in intersect)
            union = sum(elem_counter[elem] ** 0.5 for elem in union)
        else:
            intersect = len(intersect)
            union = len(union)
        jaccard = float(intersect) / union if union else 0.0
        overlap = {'intersect': intersect, 'union': union, 'jaccard': jaccard}
        return overlap

    @staticmethod
    def connect(graph, weight_metric='jaccard', weight_genes=True):
        if weight_genes:
            gene_counter = DiseasomeMaker.get_gene_counter(graph)
        else:
            gene_counter = None
        for node_0, node_1 in itertools.combinations(graph, 2):
            data_0 = graph.node[node_0]
            data_1 = graph.node[node_1]
            genes_0 = data_0['genes']
            genes_1 = data_1['genes']
            edge_metrics = DiseasomeMaker.get_overlap(genes_0, genes_1, gene_counter)
            weight = edge_metrics[weight_metric]
            if weight:
                graph.add_edge(node_0, node_1, weight=weight)

        # Remove diseases without any gene overlap with other diseases
        unconnected_nodes = {node for node, degree in graph.degree_iter() if not degree}
        graph.remove_nodes_from(unconnected_nodes)
        #graph = graph.subgraph(keep_nodes) # if want to leave original unmodified

        return graph

    @staticmethod
    def get_adjacency_matrix(graph):
        matrix = networkx.adjacency_matrix(graph)
        return matrix

    @staticmethod
    def get_seed_distance(graph, genes, weight_metric='jaccard', weight_genes=True):
        if weight_genes:
            gene_counter = DiseasomeMaker.get_gene_counter(graph)
        else:
            gene_counter = None
        genes = set(genes)
        seed_list = list()
        for data in graph.node.values():
            overlap = DiseasomeMaker.get_overlap(genes, data['genes'], gene_counter)
            seed_list.append(overlap[weight_metric])
        seed = numpy.array(seed_list)
        return seed

    @staticmethod
    def add_node_attribute(graph, key, node_to_value, default=None):
        for node, data in graph.nodes_iter(data=True):
            value = node_to_value.get(node, default)
            data[key] = value

    @staticmethod
    def graph_minus(graph, nodes):
        all_nodes = graph.nodes()
        return graph.subgraph(set(all_nodes) - set(nodes))

    @staticmethod
    def save_as_gml(diseasome, path):
        diseasome = diseasome.copy()
        doid_code_to_name = dict()
        for node, data in diseasome.nodes_iter(data=True):
            data['genes'] = '|'.join(sorted(gene.symbol for gene in data['genes']))
            #del data['genes']
            data['category'] = data['category'] or ''
            data['doid_code'] = node
            doid_code_to_name[node] = data['name']
            del data['name']
        networkx.relabel_nodes(diseasome, doid_code_to_name, copy=False)
        networkx.write_gml(diseasome, path)

    @staticmethod
    def save_flat_txt(diseasome, path):
        write_file = open(path, 'w')
        fieldnames = ['doid_code', 'name', 'category', 'n_genes', 'genes']
        writer = csv.DictWriter(write_file, delimiter='\t', fieldnames=fieldnames)
        writer.writeheader()
        for node, data in diseasome.nodes_iter(data=True):
            genes = data['genes']
            row = {'doid_code': node, 'name': data['name'],
                   'category': data['category'], 'n_genes': len(genes),
                   'genes': '|'.join(sorted(gene.symbol for gene in genes))}
            writer.writerow(row)
        write_file.close()

    @staticmethod
    def get_stats(diseasome):
        print 'Number of nodes: {}'.format(diseasome.order())
        print 'Number of edges: {}'.format(diseasome.size())
        total_annotations = 0
        distinct_genes = set()
        for data in diseasome.node.values():
            genes = data['genes']
            total_annotations += len(genes)
            distinct_genes |= genes
        print 'Number of annotated genes: {}'.format(total_annotations)
        print 'Number of distinct genes: {}'.format(len(distinct_genes))


if __name__ == '__main__':
    import pprint
    dm = DiseasomeMaker()
    #dm.get_gene_annotations()
    #pprint.pprint()
    code = 'DOID:9256'
    code = 'DOID:1520' # colon carcinoma
    print dm.get_gwas_catalog_annotations().get(code)
    print dm.get_omim_annotations().get(code)
