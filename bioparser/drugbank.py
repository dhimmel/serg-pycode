import xml.etree.ElementTree
import sys
import os
import re

import data

class DrugBank(object):
    """Reads and parses the drugbank xml file. Extracts drugs and targets."""
    
    """Create drug nodes from a drugbank xml download. Uses xml.etree.ElementTree
    to iteratively parse the input file. Users should ensure that this function
    creates the same number of drugs as listed at http://www.drugbank.ca/stats.
    Each tag in singular_drug_attr indicates an element under the drug in
    which the element text is the desired value. Each tag in plural_drug_attr
    indicates an element where the desired value is a list of the text values
    of all of the element's subelements.
    
    XML parsing stops upon reaching a partners tag. This saves memory when
    parsing Drugbank 3.0 where the xml file contains drugs followed by two
    large partners elements. However future drugbanks may not follow the same
    xml ordering.
    """

    def __init__(self, drugbank_dir=None):
        """Download xml full database: http://www.drugbank.ca/downloads"""
        if not drugbank_dir:
            drugbank_dir = data.current_path('drugbank')
        self.drugbank_dir = drugbank_dir
        self.xml_path = os.path.join(drugbank_dir, 'drugbank.xml')
        self.xml_namespace = '{http://drugbank.ca}'
        self.id_to_drug = dict()
        self.id_to_partner = dict()
        self.drugs = list()
        self.partners = list()
        self.targets = list() # dict describing the link between a drug and partner
    
    def read(self, parse_drugs=True, parse_partners=True, parse_targets=True):
        """Read an xml file to populate self.drugs, self.partners, and self.targets."""
        print "Initiating drug bank reading: ",
        sys.stdout.flush()
        xml_file = open(self.xml_path)
        iter_etree = xml.etree.ElementTree.iterparse(xml_file, ['start', 'end'])
        root = None
        for event, elem in iter_etree:
            if event == "start" and root is None:
                self.root = elem # the first element is root
            if elem.tag == '{http://drugbank.ca}partners' and not parse_partners:
                break
            if event != 'end':
                continue
           
            # Read a drug
            if (elem.tag == '{http://drugbank.ca}drug' and elem.get('type') is not None):
                if parse_drugs:
                    drug = self.parse_drug(elem)
                    self.drugs.append(drug)
                if parse_targets:
                    if not parse_drugs:
                        drug = self.extract_elem_info(elem, singular_tags=['drugbank-id'])
                        drug = self.replace_dashes_in_dict_keys(drug)
                    targets_elem = elem.find('{http://drugbank.ca}targets')
                    for target_elem in targets_elem:
                        target = self.parse_target(drug['drugbank_id'], target_elem)
                        self.targets.append(target)
                self.root.clear()
            
            # Read a partner
            if (elem.tag == '{http://drugbank.ca}partner'):
                if parse_partners:
                    partner = self.parse_partner(elem)
                    self.partners.append(partner)
                self.root.clear()
        xml_file.close()
        print 'complete.'
        
    def parse_drug(self, elem):
        """Parse a drug element and return a dict."""
        
        attributes = ['type']
        singular_tags = ['drugbank-id', 'name', 'cas-number', 'indication']
        plural_tags = ['synonyms', 'brands', 'groups']
        
        name_to_fxn = {'brands': self.remove_bracketed_or_parenthesized_text,
                       'synonyms': self.remove_bracketed_or_parenthesized_text}
        
        drug_dict = self.extract_elem_info(elem, attributes, singular_tags, plural_tags, name_to_fxn)        
        
        try:
            indic = drug_dict['indication']
            indic = indic.replace('"', "'")
            drug_dict['indication'] = indic
        except KeyError:
            pass
        
        # Get All external id entries
        ext_ids = elem.find('{http://drugbank.ca}external-identifiers')
        for ext_id in ext_ids:
            resource = ext_id.findtext('{http://drugbank.ca}resource')
            resource = resource.replace(' ', '_')
            resource = self.remove_bracketed_or_parenthesized_text(resource)
            id = ext_id.findtext('{http://drugbank.ca}identifier')
            drug_dict[resource] = id
        
        drug_dict = self.replace_dashes_in_dict_keys(drug_dict)
        return drug_dict
    
    def parse_target(self, drugbank_id, elem):
        """Parse a target element. drugbank_id is the drugbank-id of the drug
        possessing the target given by elem."""
        attributes = ['partner']
        singular_tags = ['known-action']
        plural_tags = ['actions']
        target = self.extract_elem_info(elem, attributes, singular_tags, plural_tags)
        target['drugbank-id'] = drugbank_id
        target = self.replace_dashes_in_dict_keys(target)
        return target
    
    def parse_partner(self, elem):
        """Parse a partner element and return a dict."""
        
        attributes = ['id']
        singular_tags = ['name', 'gene-name']
        plural_tags = []
        name_to_fxn = dict()
        
        partner_dict = self.extract_elem_info(elem, attributes, singular_tags, plural_tags, name_to_fxn)        
        
        species_elem = elem.find(self.xml_namespace + 'species')
        partner_dict['species'] = species_elem.findtext(self.xml_namespace + 'name')

        partner_dict = self.replace_dashes_in_dict_keys(partner_dict)
        return partner_dict

    def extract_elem_info(self, elem, attributes=list(), singular_tags=list(), plural_tags=list(), name_to_fxn=dict()):
        """Extract text from an xml element and return the desired text
        in a tag to text dicitonary. plural_tags indicates tags which have
        immediate subelements with the desired text. plural tags create a list
        of text from subelements. name_to_fxn indicates links names to a function
        to be applied to individual text strings based on the texts tag. 
        """
        elem_dict = dict()
        for attribute in attributes:
            elem_dict[attribute] = elem.get(attribute)
            
        for tag in singular_tags:
            text = elem.findtext(self.xml_namespace + tag)
            text = text.encode('utf-8')
            if text: elem_dict[tag] = text
            
        for tag in plural_tags:
            plural_elem = elem.find(self.xml_namespace + tag)
            texts = [sub_elem.text.encode('utf-8') for sub_elem in plural_elem]
            # brands include manufacturers in parenthesis. brands and synonyms
            # include language in brackets.
            texts = map(self.remove_bracketed_or_parenthesized_text, texts)
            if texts: elem_dict[tag] = texts
        
        for name in name_to_fxn.keys():
            if name in elem_dict:
                value = elem_dict[name]
                fxn = name_to_fxn[name]
                if isinstance(value, list):
                    value = map(fxn, value)
                else:
                    value = fxn(value)
                elem_dict[name] = value
        
        return elem_dict
    
    @staticmethod
    def remove_bracketed_or_parenthesized_text(s):
        s = re.sub(r'\([^)]*\)', '', s) # parenthesized text
        s = re.sub(r'\[[^]]*\]', '', s) # bracketed text
        s = s.strip()
        s = s.strip('_')
        return s
    
    @staticmethod
    def replace_dashes_in_dict_keys(d, replacement='_'):
        return {key.replace('-', replacement): value for key, value in d.items()}
    
    def get_name_to_drug(self):
        name_to_drug = dict()
        for drug in self.drugs:
            names = [drug['name']] + drug.get('synonyms', []) + drug.get('brands', [])
            name_to_drug.update(dict.fromkeys(names, drug))
        return name_to_drug
    
    def get_id_to_drug(self):
        """Returns a list of dictionaries each representing a drug."""
        return {drug['drugbank_id']: drug for drug in self.drugs}

    def get_id_to_partner(self):
        """Returns a list of dictionaries each representing a partner."""
        return {partner['id']: partner for partner in self.partners}

if __name__ == "__main__":
    db = DrugBank()
    db.read()
    db.drugs
    

