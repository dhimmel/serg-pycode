import os
import csv
import re


class ScreenReader(object):

    def __init__(self, directory):
        self.directory = directory

    def read_screening_results(self):
        """
        MBP - Measures differentiation

        """
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
            if cas_value:
                cas_list = re.findall(cas_pattern, cas_value)
            else:
                cas_list = list()
                print row
            row['cas_list'] = cas_list

            rows.append(row)

        read_file.close()
        return rows

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
            elif value_mean <= control_mean - control_sem:
                negatives.append(row)
            else:
                omitted.append(row)
        summary = 'Positives: {}, Negatives: {}, Omissions: {}'.format(*map(len, (positives, negatives, omitted)))
        print summary
        return positives, negatives, omitted




if __name__ == '__main__':

    directory = '/home/dhimmels/Documents/serg/remyelination/'
    reader = ScreenReader(directory)
    reader.summarize()