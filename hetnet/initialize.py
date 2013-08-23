import os
import shutil
import argparse

def initialize_network(directory, overwrite=False):
    """ """
    if overwrite and os.path.exists(directory):
        shutil.rmtree(directory)
    if not os.path.isdir(directory):
        os.mkdir(directory)
    
    sub_directories = 'graph', 'learning-edges', 'masks', 'metapaths', 'results'
    for sub_directory in sub_directories:
        path = os.path.join(directory, sub_directory)
        if not os.path.isdir(path):
            os.mkdir(path)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--network-dir', type=os.path.expanduser)
    parser.add_argument('--overwrite', action='store_true')
    args = parser.parse_args()

    # Create network directories    
    initialize_network(args.network_dir, args.overwrite)
