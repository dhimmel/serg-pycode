import csv
import json
import os
import urllib
import urllib2

import data

class Querier(object):
    
    def __init__(self, gxa_dir): #, gxa_dir, ontology):
        #self.gxa_dir = gxa_dir
        #self.ontology = None
        self.gxa_dir = gxa_dir
        self.atlas_address = 'www.ebi.ac.uk/gxa/api/v1'
    
    def create_dirs(self):
        if not os.path.isdir(gxa_dir):
            os.mkdir(gxa_dir)
        for directory in 'queries', 'processed':
            path = os.path.join(self.gxa_dir, directory)
            if not os.path.isdir(path):
                os.mkdir(path)
    
    def api_query(self, query_dict):
        """Query the GXA API for a parsed json formatted output."""
        for key, value in query_dict.items():
            #query_dict[key] = value.replace(' ', '+')
            query_dict[key] = urllib.quote_plus(value)
        query = '&'.join('='.join(pair) for pair in query_dict.items())
        query = 'http://' + self.atlas_address + '?' + query + '&format=json'
        results = list()
        
        starting_row = 0
        nrows = 200
        incomplete = True
        while incomplete:
            data = self.api_segment_query(query, starting_row, nrows)
            results.extend(data['results'])
            
            starting_row += nrows
            if starting_row >= data['totalResults']:
                incomplete = False
        return results
    
    def api_segment_query(self, query, starting_row=0, nrows=200):
        """Query GXA API with row subsets."""
        url = query + '&start=' + str(starting_row) + '&rows=' + str(nrows)
        while True:
            try:
                usock = urllib2.urlopen(url)
                data = json.load(usock)
                usock.close()
                return data
            except urllib2.HTTPError:
                print 'urllib2.HTTPError'
    
    def parse_results(self, results):
        """
        Indented JSON for viewing:
        http://www.ebi.ac.uk/gxa/api/v1?updownInEfo=EFO_0003885&species=Homo+Sapiens&format=json&start=0&rows=200&indent
        """
        result_dict = dict()
        for results_by_gene in results:
            #results_by_gene['gene']['id']
            #results_by_gene['gene']['ensemblGeneId'] # KeyError: 'ensemblGeneId' for EFO_0000536 
            gene_name = results_by_gene['gene']['name']
            
            for results_by_factor in results_by_gene['expressions']:
                efo_id = results_by_factor['efoId']
                pval_dict = dict()
                for results_by_experiment in results_by_factor['experiments']:
                    pvalue = results_by_experiment['pvalue']
                    direction = results_by_experiment['updn']
                    pval_dict.setdefault(direction, list()).append(pvalue)
                #pval_dict = self.min_p(pval_dict)
                result_dict.setdefault(efo_id, dict())[gene_name] = pval_dict
        return result_dict
    
    def sidak(self, p, n):
        """Minimum p-value with a Sidak adjustment."""
        return 1 - (1.0 - p) ** n

    def min_p(self, pval_dict, round_digits=None):
        """Minimum p-value with a Sidak adjustment."""
        n = sum(map(len, pval_dict.values()))
        out_dict = dict()
        for direction in 'UP', 'DOWN':
            try:
                p = min(pval_dict[direction])
                p = self.sidak(p, n)
                if round_digits is not None:
                    p = round(p, round_digits)
                out_dict[direction] = p
            except KeyError:
                out_dict[direction] = 1.0
        return out_dict

    def process_result_dict(self, result_dict):
        processed_dict = dict()
        for efo_id, results_by_factor in result_dict.items():
            processed_dict[efo_id] = dict()
            for gene, pval_dict in results_by_factor.items():
                processed_dict[efo_id][gene] = self.min_p(pval_dict, round_digits=5)
        return processed_dict
    
    def query_factors(self, efo_ids):
        self.create_dirs()
        efo_graph = data.Data().efo.get_graph()
        for efo_query_id in efo_ids:
            print 'Querying ', efo_query_id, '-', efo_graph.node[query]['name']
            query_dict = {'anyInEfo': efo_query_id, 'species': 'Homo Sapiens'}
            results = self.api_query(query_dict)
            
            # Save unprocessed but still json queries
            path = os.path.join(gxa_dir, 'queries', efo_query_id + '.json')
            with open(path, 'w') as json_file:
                json.dump(results, json_file, indent=2)
            
            # Save processed files
            result_dict = self.parse_results(results)
            fieldnames = 'gene', 'pval_down', 'pval_up'
            processed_result_dict = self.process_result_dict(result_dict)
            for efo_id, gene_dict in processed_result_dict.iteritems():
                path = os.path.join(gxa_dir, 'processed', efo_id + '.txt')
                f = open(path, 'w')
                writer = csv.writer(f, delimiter='\t')
                writer.writerow(fieldnames)
                for gene, pval_dict in gene_dict.iteritems():
                    row = gene, pval_dict['DOWN'], pval_dict['UP']
                    writer.writerow(row)
                f.close()
        print 'GXA queries complete. Results saved to file.'
        
if __name__ =='__main__':
    gxa_dir = '/home/dhimmels/Documents/serg/data-sources/gxa/130510'
    qq = Querier(gxa_dir)

    compounds = set(data.Data().efo.gxa_query_compounds())
    diseases = set(data.Data().efo.gxa_query_diseases())

    # Remove already queried
    qq.create_dirs()
    path = os.path.join(gxa_dir, 'queries')
    queried = {name.split('.json')[0] for name in os.listdir(path)}
    for term_set in compounds, diseases:
        term_set -= queried
    
    qq.query_factors(compounds)
    qq.query_factors(diseases)
    
    # Save file with names of efo_terms that had differential expression data
    path = os.path.join(gxa_dir, 'processed')
    terms_with_data = {name.split('.txt')[0] for name in os.listdir(path)}
    write_efo_terms = data.Data().efo.write_terms
    compounds = set(data.Data().efo.gxa_query_compounds())
    diseases = set(data.Data().efo.gxa_query_diseases())
    write_efo_terms(terms_with_data & compounds, os.path.join(gxa_dir, 'compounds-with-expression.txt'))
    write_efo_terms(terms_with_data & diseases, os.path.join(gxa_dir, 'diseases-with-expression.txt'))
    
    # See line numbers for processed files in unix from the processed direcetory
    # find . -type f | xargs wc -l
    