import os
import csv
import math

import data

class GNF(object):

    def __init__(self, directory=None):
        if not directory:
            directory = data.source_data_dir('gnf')
        self.directory = directory
    
    def expression_generator(self):
        symbol_to_gene = data.Data().hgnc.get_symbol_to_gene()

        path = os.path.join(self.directory, 'created', 'expression-bto.txt')
        read_file = open(path)
        reader = csv.reader(read_file, delimiter='\t')
        bto_codes = reader.next()
        for row in reader:
            row = [elem if elem != '' else None for elem in row]
            symbol = row.pop(0)
            gene = symbol_to_gene.get(symbol)
            if not gene:
                continue
            for bto_code, value in zip(bto_codes, row):
                value = float(value)
                log_expr = math.log(value)
                expression = {'gene': gene, 'bto_id': bto_code, 'expr': value, 'log_expr': log_expr}
                yield expression
        read_file.close()

if __name__ == '__main__':
    import pprint
    gnf = GNF()
    for expression in gnf.expression_generator():
        pprint.pprint(expression)
