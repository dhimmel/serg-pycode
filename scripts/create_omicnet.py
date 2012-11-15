# Daniel Himmelstein <daniel.himmelstein@gmail.com>
import os
import shutil
import argparse

import omicnet
import omictools
import gml

node_kinds = omicnet.OmicNet.node_kinds
rel_kinds = omicnet.OmicNet.rel_kinds


parser = argparse.ArgumentParser("Directs the (re)creation or expansion of OmicNet")
parser.add_argument('--neo-db-dir', type=os.path.expanduser, default=
    '/scratch/omicnet2.db')
parser.add_argument('--data-dir', type=os.path.expanduser, default=
    '~/Documents/serg/omicnet/data')
parser.add_argument('--output-dir', type=os.path.expanduser, default=
    '~/Documents/serg/omicnet/output')

parser.add_argument('--umls-metathesaurus-dir', type=os.path.expanduser, default=
    '~/Documents/serg/omicnet/data/umls/2012AA/META')

parser.add_argument('--delete', nargs = '+', default = [],
                    choices = ['graph'] + node_kinds + rel_kinds)
parser.add_argument('--export-gml', default=False)
args = parser.parse_args()


if __name__ == '__main__':
    
    global data_dir
    data_dir = args.data_dir
    
    print ''.center(79, '#')
    
    if 'graph' in args.delete:
        print ' Deleting all neo4j db directory contents'.center(79, '#')
        for root, dirs, files in os.walk(args.neo_db_dir):
            for f in files:
                os.unlink(os.path.join(root, f))
            for d in dirs:
                shutil.rmtree(os.path.join(root, d))

    try:
        
        g = omicnet.OmicNet(args.neo_db_dir)

        # Delete nodes or relationships of kinds specified by --delete
        if args.delete and 'graph' not in args.delete:
            print ' Deleting '.center(79, '#')
            for kind in args.delete:
                if kind in node_kinds:
                    g.delete_all_nodes_of_kind(kind)
                elif kind in rel_kinds:
                    g.delete_all_relationships_of_kind(kind)
        
        # Create gene nodes.
        if not g.get_nodes_of_kind('gene'):
            current = omictools.current_data_dir('hgnc')
            g.create_genes(current)
        
        # Create drug nodes.
        if not g.get_nodes_of_kind('drug'):
            current = omictools.current_data_dir('drugbank')
            g.create_drugs(current)
        
        # Create disease nodes.
        if not g.get_nodes_of_kind('disease'):
            g.create_diseases(args.umls_metathesaurus_dir)
        
        # Create disease to gene relationships from GWAS.
        if not g.get_relationships_of_kind('disease2gene_gwas'):
            current = omictools.current_data_dir('gwas-catalog')
            previous = omictools.preceding_data_dir('gwas-catalog')
            g.create_disease2gene_gwas(current, previous)
        
        # Create drug to gene relationships from drug targers.
        if not g.get_relationships_of_kind('drug2gene_target'):
            current = omictools.current_data_dir('drugbank')
            g.create_drug2gene_target(args.drug_bank_dir)

        # Create disease to drug relationships (indications) from PREDICT Nature paper.
        if not g.get_relationships_of_kind('disease2drug_predict'):
            current = os.path.join(data_dir, 'nature-predict')
            g.create_disease2drug_predict(current)

        # Create gene to gene relationships from IMP predications.
        if not g.get_relationships_of_kind('gene2gene_imp'):
            current = omictools.current_data_dir('imp')
            g.create_gene2gene_imp(current)
        
        g.print_node_counts()
        g.print_relationship_counts()
        
        # Export GML file
        if args.export_gml:
            '/home/dhimmels/Documents/serg/omicnet/omicnet-py2.gml'
            neo_gml = gml.Neo4jWriter(args.export_gml)
            neo_gml.export_graph(g.get_connected_nodes(), g.db.relationships)
            
    finally:
        
        g.shutdown()
