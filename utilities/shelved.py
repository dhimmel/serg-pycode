import os
import shelve

class Shelved(object):
    
    def __init__(self, path, names):
        """
        shelve_name_to_path is a dictionary with keys representing shelve
        identifiers and values representing the path to the shelf.
        """
        #self.shelve_dir = path
        self.shelve_name_to_path = {name: os.path.join(path, name + '.shelve') for name in names}
        self.shelves = dict().fromkeys(names)
    
    def shelve_files_exist(self):
        paths = self.shelve_name_to_path.values()
        return all(os.path.isfile(path) for path in paths)
    
    def delete_all_shelves(self):
        paths = self.shelve_name_to_path.values()
        for path in paths:
            if os.path.isfile(path):
                os.remove(path)
    
    def __enter__(self):
        for shelve_name in self.shelves:
            path = self.shelve_name_to_path[shelve_name]
            shelve_ = shelve.open(path)
            self.shelves[shelve_name] = shelve_
        return self
    
    def __exit__(self, *args, **kwargs):
        for shelve_ in self.shelves.values():
            shelve_.close()
    
    def sync(self):
        for shelve_ in self.shelves.values():
            shelve_.sync()

