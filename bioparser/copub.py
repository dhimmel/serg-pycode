import os

import suds.client

import data
#import scripts.ashg13


#CoPub Web Services
# http://services.nbic.nl/copub5/help/CoPub_web_services.html
# http://services.nbic.nl/copub5/copub.wsdl.html
wsdl_url = 'http://services.nbic.nl/copub5/copub.wsdl'
client = suds.client.Client(wsdl_url)
client.service.version()
categories = map(str, client.service.get_categories('keyword')[0])


#class CoPubQuerier(object):


def get_top_keyword(query, category):
    """example category: 'tissue', 'disease', 'pathway'"""
    assert category in categories
    output = client.service.get_keywords(query=query, category=category, max_results=1)
    top_result = output[0][0]
    bi_id = int(top_result['bi_id'])
    preferred_name = top_result['preferred_name']
    return bi_id, preferred_name

def get_r(term_a, term_b):
    """Get r scaled score between two terms. Terms should be bi_ids."""
    output = client.service.get_references(bi_id1=term_a, bi_id2=term_b, max_results=1)
    r_scaled_score = float(output['r_scaled_score'])
    return r_scaled_score

#term_a = get_top_keyword('multiple sclerosis', 'disease')[0]
#term_b = get_top_keyword('brain', 'tissue')[0]
#get_r(term_a, term_b)




#client.service.get_keywords(query='brain', category='tissue', max_results=2)
#client.service.get_references(bi_id1=term_a, bi_id2=term_b, max_results=1)
#client.service.get_literature_neighbours(bi_ids={'bi_id':term_a}, categories={'category':'tissue'}, max_results=1000, r_scaled_score_threshold=25)
#scripts.ashg13.calculate_disease_subset()