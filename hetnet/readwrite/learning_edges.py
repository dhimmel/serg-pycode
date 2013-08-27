




def write_text(edge_to_exclusions, path):
    fieldnames = ['edge_str', 'edge_id', 'excluded_edges']
    with open(path, 'w') as write_file:
        writer = csv.writer(write_file, delimiter='\t')
        writer.writerow(fieldnames)
        for edge, exclusions in edge_to_exclusions.iteritems():
            writer.writerow([edge, edge.get_id(), [e.get_id() for e in exclusions]])
    
    
    
def read_text(path, graph):
    pass