import ast
import argparse
import csv
import os
import collections
import itertools
import logging
import ConfigParser
import pprint
import gzip

import bioparser.gxa
import bioparser.data
import bioparser.gwas_plus
import hetnet
import hetnet.agents
import hetnet.readwrite.graph

from projects.gene_disease_hetnet.data_integration import copub_analysis



def create_graph(associations_path, doidprocess_path, pathophys_path, partition_path, exclude_testing=False):
    data = bioparser.data.Data()
    doid_remove, doid_pop = bioparser.gwas_plus.GwasCatalog.read_ontprocess_info(doidprocess_path)
    exclude_doids = doid_remove | set(doid_pop)

    msigdb = bioparser.data.Data().msigdb
    msig_set_types = msigdb.abbrev_to_name.keys()
    # http://www.broadinstitute.org/gsea/msigdb/collections.jsp
    #msig_set_types = ['c1.all', 'c2.cgp', 'c2.cp.all', 'c3.mir', 'c3.tft',
    #                  'c4.cgn', 'c4.cm', 'c5.bp', 'c5.cc', 'c5.mf', 'c6.all', 'c7.all']

    # Define and initialize networkx graph
    metaedge_tuples = [('disease', 'gene', 'association', 'both'),
                       ('gene', 'gene', 'interaction', 'both'),
                       ('gene', 'tissue', 'expression', 'both'),
                       ('disease', 'tissue', 'cooccurrence', 'both'),
                       ('disease', 'pathophysiology', 'membership', 'both')]
    metaedge_tuples.extend([('gene', set_type, 'membership', 'both') for set_type in msig_set_types])
    metagraph = hetnet.MetaGraph.from_edge_tuples(metaedge_tuples)
    graph = hetnet.Graph(metagraph)

    # Add genes from HGNC
    logging.info('Adding HGNC gene nodes.')
    for gene in data.hgnc.get_genes():
        if not gene.coding:
            continue
        node_data = {'name': gene.name}
        graph.add_node(gene.symbol, 'gene', node_data)

    # Add tissues from BTO
    logging.info('Adding BTO tissue nodes.')
    bto_graph = data.bto.get_animal_tissue_subgraph()
    for bto_id, nx_data in bto_graph.nodes(data=True):
        node_data = {'name': nx_data['name']}
        graph.add_node(bto_id, 'tissue', node_data)

    # Add diseases from DOID
    logging.info('Adding DOID disease nodes.')
    doid_onto = data.doid.get_ontology()
    for doid_id, nx_data in doid_onto.graph.nodes(data=True):
        if doid_id in exclude_doids:
            continue
        node_data = {'name': nx_data['name']}
        graph.add_node(doid_id, 'disease', node_data)

    # Add pathophysiology nodes
    exclude_pathophys = {'unspecific', 'ideopathic'}
    with open(pathophys_path) as read_file:
        reader = csv.DictReader(read_file, delimiter='\t')
        pathophys_rows = list(reader)
    pathophys_rows = [row for row in pathophys_rows
                      if row['pathophysiology'] not in exclude_pathophys]
    pathophys_rows = [row for row in pathophys_rows
                      if row['doid_code'] not in exclude_doids]
    pathophys_terms = {row['pathophysiology'] for row in pathophys_rows}
    for pathophys_term in pathophys_terms:
        graph.add_node(pathophys_term, 'pathophysiology')

    # Add (disease, pathophysiology, membership, both) edges
    for pathophys_row in pathophys_rows:
        doid_code = pathophys_row['doid_code']
        pathophys_term = pathophys_row['pathophysiology']
        graph.add_edge(doid_code, pathophys_term, 'membership', 'both')


    # Add (disease, gene, association, both) edges
    with gzip.open(partition_path) as part_file:
        reader = csv.DictReader(part_file, delimiter='\t')
        part_rows = list(reader)
    assoc_to_part = {(row['disease_code'], row['gene_symbol']): row['part']
                        for row in part_rows if row['status'] != 'negative'}

    logging.info('Adding GWAS catalog disease-gene associations.')
    associations_file = open(associations_path)
    associations_reader = csv.DictReader(associations_file, delimiter='\t')
    doids_with_associations = set()
    for association in associations_reader:
        disease_code = association['disease_code']
        gene_symbol = association['gene_symbol']
        assoc_tuple = disease_code, gene_symbol
        if association['status'] != 'assoc_high':
            continue
        part = assoc_to_part.get(assoc_tuple, 'excluded')
        if exclude_testing and part == 'test':
            continue
        graph.add_edge(disease_code, gene_symbol, 'association', 'both')
        doids_with_associations.add(disease_code)
    associations_file.close()

    # Add (gene, tissue, expression, both) edges
    logging.info('Adding GNF gene-tissue expression.')
    log10_expr_cutoff = 1.3
    logging.info('log10_expression_cutoff: {}'.format(log10_expr_cutoff))
    expressions = data.gnf.expression_generator()
    for expression in expressions:
        if expression['gene'].locus_group != 'protein-coding gene':
            continue
        if expression['log10_expr'] < log10_expr_cutoff:
            continue
        edge_data = {'log10_expr': expression['log10_expr']}
        graph.add_edge(expression['bto_id'], expression['gene'].symbol, 'expression', 'both', edge_data)


    # Add (gene, gene, interaction, both) edges
    logging.info('Adding ppiTrim gene-gene interaction.')
    interactions = data.ppitrim.collapsed_binary_interactions()
    edge_keys = ['pubmed', 'method_id', 'interaction_type_id']
    for interaction in interactions:
        edge_data = {k: v for k, v in interaction.items() if k in edge_keys}
        source = interaction['source'].symbol
        target = interaction['target'].symbol
        try:
            graph.add_edge(source, target, 'interaction', 'both', edge_data)
        except KeyError:
            pass


    # Add (disease, tissue, cooccurrence, both)
    logging.info('Adding CoPub disease-tissue cooccurrence.')
    r_scaled_cutoff = 33
    logging.info('r_scaled_cutoff: {}'.format(r_scaled_cutoff))
    coocc_gen = copub_analysis.doid_bto_cooccurrence_generator()
    for row in coocc_gen:
        doid_id = row['doid_id']
        if doid_id in exclude_doids:
            continue
        if doid_id not in doids_with_associations:
            continue
        bto_id = row['bto_id']
        r_scaled = row['r_scaled']
        if r_scaled < r_scaled_cutoff:
            continue
        edge_data = {'r_scaled': r_scaled}
        graph.add_edge(doid_id, bto_id, 'cooccurrence', 'both', edge_data)

    # Add MSigDB gene sets
    for set_type in msig_set_types:
        for name, description, genes in msigdb.gene_set_generator(set_type):
            unique_name = 'MSigDB_{}:{}'.format(set_type, name)
            node_data = {'description': description}
            graph.add_node(unique_name, set_type, node_data)
            for gene in genes:
                try:
                    graph.add_edge(gene.symbol, unique_name, 'membership', 'both')
                except KeyError:
                    pass

    # Print metanode counter
    logging.info('MetaNode Counts')
    for metanode, nodes in graph.get_metanode_to_nodes().items():
        line = '{}: {}'.format(metanode, len(nodes))
        print line
        logging.info(line)

    # Print metaedge counter
    logging.info('MetaEdge Counts')
    for metaedge, edges in graph.get_metaedge_to_edges(exclude_inverts=True).items():
        line = '{}: {}'.format(metaedge, len(edges))
        print line
        logging.info(line)

    return graph




if __name__ == '__main__':
    # Parse the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--network-dir', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/networks/140514-all-assoc')
    parser.add_argument('--doidprocess-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/data-integration/doid-ontprocess-info.txt')
    parser.add_argument('--pathophys-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/data-integration/pathophysiology.txt')
    parser.add_argument('--partition-path', type=os.path.expanduser, default=
        '~/Documents/serg/gene-disease-hetnet/partitions.txt.gz')
    parser.add_argument('--exclude-testing', action='store_true')
    parser.add_argument('--associations-id', default='processed')
    parser.add_argument('--create', action='store_true')
    args = parser.parse_args()
    network_dir = args.network_dir
    graph_agent = hetnet.agents.GraphAgent(network_dir)
    graph_dir = graph_agent.graph_dir

    associations_path = os.path.join(
        bioparser.data.Data().gwas_plus.directory,
        args.associations_id, 'association-statuses.txt')

    if args.create:

        # Create the graph
        log_path = os.path.join(graph_dir, 'creation.log')
        logging.basicConfig(filename=log_path, level=logging.INFO,
                            filemode='w', format='%(levelname)s:%(message)s')
        graph = create_graph(associations_path=associations_path,
                             doidprocess_path=args.doidprocess_path,
                             pathophys_path=args.pathophys_path,
                             partition_path=args.partition_path,
                             exclude_testing=args.exclude_testing)

        # Save the graph
        graph_agent = hetnet.agents.GraphAgent(network_dir)
        graph_agent.set(graph)
        graph_agent.write_additional_formats()
        sif_path = os.path.join(network_dir, 'graph', 'graph.sif')
        hetnet.readwrite.graph.write_sif(graph, sif_path)
        sif_subset_path = os.path.join(network_dir, 'graph', 'graph-10k.sif')
        hetnet.readwrite.graph.write_sif(graph, sif_subset_path, max_edges=int(1e4), seed=0)
        nodetable_path = os.path.join(network_dir, 'graph', 'node_table.tsv')
        hetnet.readwrite.graph.write_nodetable(graph, nodetable_path)