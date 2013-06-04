import gzip
import csv
import os
import re
import subprocess
import argparse

import data

class DbGap(object):
    
    def __init__(self, dbgap_dir=None):
        if dbgap_dir is None:
            dbgap_dir = data.source_data_dir('dbgap')
        self.dbgap_dir = dbgap_dir
    
    def get_analyses_paths(self):
        studies_dir = os.path.join(self.dbgap_dir, 'studies')
        studies = filter(lambda x: os.path.isdir(os.path.join(studies_dir, x)), os.listdir(studies_dir))
        analyses_paths = list()
        for study in studies:
            analyses_dir = os.path.join(studies_dir, study, 'analyses')
            for file_name in os.listdir(analyses_dir):
                if not file_name.startswith(study):
                    continue
                if not file_name.endswith('.txt.gz'):
                    continue
                analysis_path = os.path.join(analyses_dir, file_name)
                analyses_paths.append(analysis_path)
        return analyses_paths
        
    
    def read_analysis_file(self, path):
        anal_file = gzip.open(path)
        line = None
        while not line or line.startswith('#'):
            line = anal_file.readline().strip()
        fieldnames = line.split('\t')
        reader = csv.DictReader(anal_file, fieldnames=fieldnames, delimiter='\t')
        rows = list()
        for row in reader:
            if row['Chr Position'] == '':
                continue
            row['Chr Position'] == int(row['Chr Position'])
            if row['P-value'] == '':
                continue
            row['P-value'] == float(row['P-value'])
            yield row
        anal_file.close()            

    def vegas_format_analyses_files(self):
        analyses_paths = self.get_analyses_paths()
        for path in analyses_paths:
            print path
            self.vegas_format_analysis_file(path)
        
    def vegas_format_analysis_file(self, path):
        anal_file = gzip.open(path)
        row_generator = self.read_analysis_file(path)

        head, tail = os.path.split(path)
        vegas_file_name = tail.replace('.txt.gz', '-vegas-input.txt')
        vegas_path = os.path.join(self.dbgap_dir, 'vegas-input', vegas_file_name)
        vegas_file = open(vegas_path, 'w')
        writer = csv.writer(vegas_file, delimiter='\t')
        for row in row_generator:
            writer.writerow([row['SNP ID'], row['P-value']])
        anal_file.close()
        vegas_file.close()
    
    def kgg_format_analyses_files(self):
        analyses_paths = self.get_analyses_paths()
        for path in analyses_paths:
            self.kgg_format_analysis_file(path)

    def kgg_format_analysis_file(self, path):
        anal_file = gzip.open(path)
        row_generator = self.read_analysis_file(path)

        head, tail = os.path.split(path)
        kgg_file_name = tail.replace('.txt.gz', '-kgg-input.txt')
        kgg_path = os.path.join(self.dbgap_dir, 'kgg-input', kgg_file_name)
        kgg_file = open(kgg_path, 'w')
        writer = csv.writer(kgg_file, delimiter='\t')

        rows = list()
        for row in row_generator:
            if any(x == '' for x in row):
                continue
            rows.append(row)
        rows.sort(key=lambda x: x[2:3])
        writer.writerows(rows)
        
        anal_file.close()
        kgg_file.close()
    
    def run_vegas(self, input_file):
        vegas_path = '/home/dhimmels/Documents/serg/programs/VEGAS/vegas'
        vegas_dir = os.path.join(self.dbgap_dir, 'vegas')
        vegas_input_dir = os.path.join(self.dbgap_dir, 'vegas-input')
        os.chdir(vegas_dir) # vegas writes output and temporary files to the directory its run from
        input_path = os.path.join(vegas_input_dir, input_file)
        output_file = input_file.replace('vegas-input.txt', 'hapmapCEU')
        #output_path = os.path.join(vegas_output_dir, output_file)
        subprocess.call([vegas_path, input_path, '-pop', 'hapmapCEU', '-out', output_file])
        
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run-vegas', nargs=2, default=False, help='Indicates Vegas Should be Run. First arg is study name, second is analysis name.')
    args = parser.parse_args()

    gap = DbGap()
    
    if args.run_vegas:
        input_file = '.'.join(args.run_vegas) + '-vegas-input.txt'
        print input_file
        gap.run_vegas(input_file)
    
    #vegas_input_dir = os.path.join(self.dbgap_dir, 'vegas-input')
    #input_files = os.listdir(vegas_input_dir)

        