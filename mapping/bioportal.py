import os
import urllib2
import xml.etree.ElementTree
import HTMLParser
import csv

import rdflib

import omictools

def get_ontology_ids():
    """Returns a ontology id to ontology name dictionary."""
    
    """
    url = 'http://bioportal.bioontology.org/ontologies'
    f = urllib2.urlopen(url)
    page_source = f.read()
    usock.close()
    """
    current_dir = omictools.current_data_dir('bioportal')
    path = os.path.join(current_dir, 'bioontology_ids.txt')
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        id_to_name = {row['id']: row['name'] for row in reader}
    return id_to_name
    
    
    
def read_mapping(path, remove_prefix=True):
    print "Reading bioportal mapping file:", path
    rdf = rdflib.graph.Graph()
    rdf.parse(path)
    predicates = ['source', 'target', 'source_ontology_id', 'target_ontology_id', 'mapping_type',
                  'mapping_source', 'mapping_source_name', 'relation']
    prefix = 'http://protege.stanford.edu/ontologies/mappings/mappings.rdfs#'
    subj_to_objs = dict()
    
    for predicate in predicates:
        uri_predicate = rdflib.URIRef(prefix + predicate)
        
        for subj, obj_ in rdf.subject_objects(uri_predicate):
            obj_ = str(obj_)
            
            # Remove unwanted prefix in source and target names.
            if remove_prefix and predicate in ['source', 'target']:
                split_strs = ['/', '#']
                for split_str in split_strs:
                    obj_ = obj_.rsplit(split_str, 1)[-1]
            
            subj_to_objs.setdefault(subj, dict())[predicate] = obj_

    mappings = subj_to_objs.values()
        
    ontology_id_to_name = get_ontology_ids()
    source_ontology_id = mappings[0]['source_ontology_id']
    target_ontology_id = mappings[0]['target_ontology_id']
    print 'source:', ontology_id_to_name[source_ontology_id]
    print 'target:', ontology_id_to_name[target_ontology_id]
    
    source_to_target = {mapping['source']: mapping['target'] for mapping in mappings}
    return source_to_target



if __name__ == '__main__':
    path = '/home/dhimmels/Downloads/meddra_to_nci.rdf'
    source_to_target = read_mapping(path)
    for source, target in source_to_target.items()[:10]:
        print source, '\t', target
        
    
    
"""
for i, ts in enumerate(rdf):
    
    print ts
    print '--------'
    if i > 10: break
    
"""