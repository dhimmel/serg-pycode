import os
import csv
import re

import pubchempy

import bioparser


class ScreenReader(object):

    def __init__(self, directory):
        self.directory = directory

    def get_rows(self):
        """
        MBP - Measures differentiation

        """
        if hasattr(self, 'rows'):
            return self.rows
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
        self.rows = rows
        return rows

    def add_canonical_smiles(self):
        for row in self.get_rows():
            smiles = ScreenReader.canonical_smiles_from_row(row)
            row['canonical_smiles'] = smiles
            if not smiles:
                print row['name'], ' no smiles found'

    @staticmethod
    def canonical_smiles_from_row(row):

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




    def summarize(self, protein='MBP'):
        """
        Assumes the control is the last row
        """
        rows = self.read_screening_results()
        control_row = rows.pop(-1)
        print control_row
        assert control_row['name'] == 'ctl(15)'
        control_mean = control_row['mean_{}'.format(protein)]
        control_sem = control_row['SEM_{}'.format(protein)]
        negatives, positives, omitted = list(), list(), list()
        for row in rows:
            value_mean = row['mean_{}'.format(protein)]
            value_sem = row['SEM_{}'.format(protein)]
            if value_mean > control_mean + control_sem:
                positives.append(row)
                row['status'] = 1
            elif value_mean <= control_mean - control_sem:
                negatives.append(row)
                row['status'] = 0
            else:
                omitted.append(row)
                row['status'] = None
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

    directory = '/home/dhimmels/Documents/serg/remyelination/'
    reader = ScreenReader(directory)
    reader.get_rows()
    reader.add_canonical_smiles()
    reader.save()
    #reader.summarize()
    #reader.save()