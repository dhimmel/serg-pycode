import os
import itertools
import suds.client

import data

#CoPub Web Services
# http://services.nbic.nl/copub5/help/CoPub_web_services.html
# http://services.nbic.nl/copub5/copub.wsdl.html

def load_client():
    global client
    global categories
    wsdl_url = 'http://services.nbic.nl/copub5/copub.wsdl'
    client = suds.client.Client(wsdl_url)
    client.service.version()
    categories = map(str, client.service.get_categories('keyword')[0])

def get_top_keyword(query, category):
    """
    example category: 'tissue', 'disease', 'pathway'
    Returns (None, None) if no match is found.
    """
    assert category in categories
    output = client.service.get_keywords(query=query, category=category, max_results=1)
    try:
        top_result = output[0][0]
    except IndexError:
        return None, None
    bi_id = int(top_result['bi_id'])
    preferred_name = top_result['preferred_name']
    return bi_id, preferred_name

def term_cooccurrence(term_a, term_b):
    """Get r scaled score between two terms. Terms should be bi_ids."""
    try:
        output = client.service.get_references(bi_id1=term_a, bi_id2=term_b, max_results=1)
        r_scaled_score = float(output['r_scaled_score'])
    except Exception:
        r_scaled_score = 0.0
    return r_scaled_score

def term_set_cooccurrences(terms_a, terms_b):
    pair_generator = itertools.product(terms_a, terms_b)
    for term_a, term_b in pair_generator:
        r_scaled_score = term_cooccurrence(term_a, term_b)
        cooccurance_tuple = term_a, term_b, r_scaled_score
        yield cooccurance_tuple
    

if __name__ =='__main__':
    load_client()

#term_a = get_top_keyword('multiple sclerosis', 'disease')[0]
#term_b = get_top_keyword('brain', 'tissue')[0]
#get_r(term_a, term_b)




#client.service.get_keywords(query='brain', category='tissue', max_results=2)
#client.service.get_references(bi_id1=term_a, bi_id2=term_b, max_results=1)
#client.service.get_literature_neighbours(bi_ids={'bi_id':term_a}, categories={'category':'tissue'}, max_results=1000, r_scaled_score_threshold=25)
#scripts.ashg13.calculate_disease_subset()