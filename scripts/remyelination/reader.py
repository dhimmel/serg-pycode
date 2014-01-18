import os
import csv
import re

import pubchempy

import bioparser


class ScreenReader(object):

    def __init__(self, directory=None):
        if not directory:
            directory = '/home/dhimmels/Documents/serg/remyelination/screen'
        self.directory = directory

    def get_raw_rows(self):
        """
        MBP - Measures differentiation

        """
        if hasattr(self, 'raw_rows'):
            return self.raw_rows
        path = os.path.join(self.directory, 'selleck-screening-results.txt')
        read_file = open(path)
        reader = csv.DictReader(read_file, delimiter='\t')
        cas_pattern = re.compile(r'([0-9]+-[0-9]+-[0-9]+)')
        float_fields = ['mean_PDGFR', 'mean_MBP', 'SEM_PDGFR', 'SEM_MBP']
        rows = list()
        for row in reader:

            # Set empty fields to None
            row = {k: v or None for k, v in row.items()}

            # Convert fields to floats
            for field in float_fields:
                value = row[field]
                if value is None:
                    continue
                row[field] = float(value)

            # Parse cas_numbers
            cas_value = row['cas_numbers']
            if cas_value and cas_value != 'N/A':
                cas_list = re.findall(cas_pattern, cas_value)
            else:
                cas_list = list()
            row['cas_list'] = cas_list
            row['cas_rn'] = cas_list[0] if cas_list else None
            rows.append(row)

        read_file.close()
        self.raw_rows = rows
        return rows

    def add_canonical_smiles(self):
        for row in self.get_raw_rows():
            smiles = ScreenReader.canonical_smiles_from_row(row)
            row['canonical_smiles'] = smiles
            if not smiles:
                print row['name'], ' no smiles found'

    @staticmethod
    def canonical_smiles_from_row(row):
        """
        Returns None if no match is found. Otherwise return the smiles
        """

        try:
            catalog_number = row['catalog_number']
            pubchem_name = '{}_Selleck'.format(catalog_number)
            pubchem_results = pubchempy.get_properties('CanonicalSMILES', pubchem_name, 'name')
            canonical_smiles = list({str(result['CanonicalSMILES']) for result in pubchem_results})
            assert len(canonical_smiles) == 1
            return canonical_smiles[0]
        except (TypeError, pubchempy.NotFoundError, AssertionError):
            # selleck catalog_number mapping failed
            pass

        try:
            cas_rn = row['cas_rn']
            pubchem_results = pubchempy.get_properties('CanonicalSMILES', cas_rn, 'name')
            canonical_smiles = list({str(result['CanonicalSMILES']) for result in pubchem_results})
            assert len(canonical_smiles) == 1
            return canonical_smiles[0]
        except (TypeError, pubchempy.NotFoundError, AssertionError):
            # cas_rn mapping failed
            pass

        return None

    def get_screened_compounds(self):

        if hasattr(self, 'compounds'):
            return self.compounds
        path = os.path.join(self.directory, 'screened_compounds_smiles_results.tdt')
        read_file = open(path)
        reader = csv.DictReader(read_file, delimiter='\t')
        float_fields = ['mean_PDGFR', 'mean_MBP', 'SEM_PDGFR', 'SEM_MBP']
        rows = list()
        for row in reader:
            # Set empty fields to None
            row = {k: v or None for k, v in row.items()}

            # Convert fields to floats
            for field in float_fields:
                value = row[field]
                if value is None:
                    continue
                row[field] = float(value)

            rows.append(row)

        read_file.close()
        self.compounds = rows
        return rows

    def get_smiles_to_compound(self):
        return {compound['canonical_smiles']: compound for compound in self.get_screened_compounds()}

    def SEA_generator(self):
        float_fields = ['Affinity Threshold (nM)', 'p-value', 'Max Tc']
        path = os.path.join(self.directory, 'screened_compounds_predictions.chembl17.ecfp4.csv')
        read_file = open(path)
        reader = csv.DictReader(read_file)
        for row in reader:
            # Convert fields to floats
            for field in float_fields:
                value = row[field]
                if value is None:
                    continue
                row[field] = float(value)
            yield row
        read_file.close()

    def get_SEA_interactions(self):
        if hasattr(self, 'interactions'):
            return self.interactions
        best_row = lambda rows: min(rows, key=lambda row: row['p-value'])
        query_to_rows = dict()
        for row in self.SEA_generator():
            target_id = row['Target ID']
            query_id = row['Query ID']
            query_target_id = query_id, target_id
            query_to_rows.setdefault(query_target_id, list()).append(row)
        query_to_row = {k: best_row(v) for k, v in query_to_rows.iteritems()}
        interactions = query_to_row.values()
        self.interactions = interactions
        return self.interactions

    def get_SEA_targets(self):
        """
        Returns tuples of targets in the form:
        (target_id, target_name, target_description)
        """
        targets = set()
        keys = ['Target ID', 'Target Name/UNIPROT', 'Target Description']
        interactions = self.get_SEA_interactions()
        for interaction in interactions:
            target = tuple(interaction[key] for key in keys)
            targets.add(target)
        return targets

    def write_targets(self):
        path = os.path.join(self.directory, 'targets-chembl17.txt')
        write_file = open(path, 'w')
        writer = csv.writer(write_file, delimiter='\t')
        writer.writerow(['target_id', 'name', 'description'])
        writer.writerows(self.get_SEA_targets())
        write_file.close()


    def binarize_screen(self, protein='MBP'):
        """
        Assumes the control is the last row
        """
        control_mean = 10.667
        control_sem = 1.96
        negatives, positives, omitted = list(), list(), list()

        PDGFR_threshold = 15.0
        MBP_threshold = 5.0
        for compound in self.get_screened_compounds():
            value_mean = compound['mean_{}'.format(protein)]
            value_sem = compound['SEM_{}'.format(protein)]
            if compound['mean_PDGFR'] < PDGFR_threshold and compound['mean_MBP'] < MBP_threshold:
                omitted.append(compound)
                compound['status'] = -1
            elif value_mean > control_mean + control_sem:
                positives.append(compound)
                compound['status'] = 1
            elif value_mean <= control_mean + control_sem:
                negatives.append(compound)
                compound['status'] = 0
            else:
                omitted.append(compound)
                compound['status'] = -1



        summary = 'Positives: {}, Negatives: {}, Omissions: {}'.format(*map(len, (positives, negatives, omitted)))
        print summary
        self.positives = positives
        self.negatives = negatives
        self.omitted = omitted
        return positives, negatives, omitted

    def save(self):
        path = os.path.join(self.directory, 'processed-screening-results.txt')
        write_file = open(path, 'w')
        fieldnames = ['name', 'catalog_number', 'cas_rn', 'canonical_smiles', 'mean_PDGFR', 'mean_MBP', 'SEM_PDGFR', 'SEM_MBP']
        writer = csv.DictWriter(write_file, fieldnames=fieldnames, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(self.rows)
        write_file.close()


def add_drugbank():

    drugbank = bioparser.data.Data().drugbank
    hgnc = bioparser.data.Data().hgnc
    symbol_to_gene = hgnc.get_symbol_to_gene()
    drugbank.read()
    for drug in drugbank:
        drug['cas_number'] = drug.get('cas_number')
        drug['groups'] = drug.get('groups', list())


    id_to_partner = drugbank.get_id_to_partner()
    for target in drugbank.targets:
        partner_id = target['partner']
        partner = id_to_partner[partner_id]
        if partner['species'] != 'Homo sapiens':
            continue
        gene = symbol_to_gene.get(partner.get('gene_name'))
        if not gene:
            continue





if __name__ == '__main__':
    import pprint

    reader = ScreenReader()
    #reader.write_targets()

    #pprint.pprint(reader.SEA_condenser().items()[0:5])
    #reader.get_rows()
    #reader.add_canonical_smiles()
    #reader.save()
    #reader.summarize()
    #reader.save()