import os

import utilities.omictools

import ctd
import drugbank
import dvd
import efo
import frimp
import gwas_catalog
import hgnc
import ipa
import iref
import meddra
import metathesaurus
import nature_predict
import nci_thesaurus
import omim
import sider
import pdn



global data_dir
data_dir = os.path.expanduser('~/Documents/serg/data-sources/')

def sorted_dir_content(directory, n=0, path_type='directories', require_dated_format=False):
    """
    Retrieves the n (zero indexed) most recent content inside the
    directory given by path. A valid dated format is a six digit 
    numeric string with format YYMMDD. This date should represent the date of
    retrieval for the contained data. If n is invalid, None is returned.
    Otherwise the full directory or file path is returned.
    """
    
    assert path_type in ['directories', 'files', 'both']
    
    contents = os.listdir(directory)
    
    if path_type == 'directories':
        contents = [name for name in contents
                    if os.path.isdir(os.path.join(directory, name))]
    if path_type == 'files':
        contents = [name for name in contents
                    if os.path.isfile(os.path.join(directory, name))]

    if require_dated_format:
        is_dated_format = lambda s: len(s) == 6 and s.isdigit()
        contents = filter(is_dated_format, contents)

    contents.sort(reverse=True)
    
    try:
        name = contents[n]
        return os.path.join(directory, name)
    
    except IndexError:
        print directory, 'does not contain a sufficient number of \
        contents to retrieve file/directory', n
        return None

def source_data_dir(source):
    """source is the name of the folder containing the desired data resource."""
    return os.path.join(data_dir, source)

def current_path(source, **kwargs):
    path = source_data_dir(source)
    return sorted_dir_content(path, n=0, **kwargs)
    
def preceding_path(source, **kwargs):
    path = source_data_dir(source)
    return sorted_dir_content(path, n=1, **kwargs)



@utilities.omictools.singleton
class Data:
    
    def __init__(self):
        #print 'initializing DataObjects'
        self.ctd = ctd.CTD(current_path('ctd'))
        self.drugbank = drugbank.DrugBank(current_path('drugbank'))
        self.dvd = dvd.DvD()
        self.efo = efo.EFO()
        self.frimp = frimp.IMP(current_path('imp'))
        self.gwas_catalog = gwas_catalog.GwasCatalog()
        self.hgnc = hgnc.HGNC()
        self.iref = iref.iRefIndex()
        self.ipa = ipa.IPA()
        self.meddra = meddra.MedDRA(current_path('meddra'))
        self.metathesaurus = metathesaurus.Metathesaurus()
        self.omim = omim.OMIM(current_path('omim'))
        self.nature_predict = nature_predict.NaturePredict(source_data_dir('nature-predict'))
        self.nci = nci_thesaurus.NCIOntology(current_path('nci', path_type='files'))
        self.sider = sider.SIDER(current_path('sider'))
        self.pdn = pdn.PDN()

