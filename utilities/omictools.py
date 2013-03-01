import os
import csv
import shelve

#import bioparser.metathesaurus as metathesaurus

class With(object):
    
    def __init__(self, *args):
        self.args = args
        
    def __enter__(self):
        return tuple(arg.__enter__() for arg in self.args)
        
    def __exit__(self, *args, **kwargs):
        for arg in self.args:
            arg.__exit__(*args, **kwargs)


class CommentedFile:
    
    def __init__(self, f, comment_string='#'):
        self.f = f
        self.comment_string = comment_string
        
    def next(self):
        line = self.f.next()
        while line.startswith(self.comment_string):
            line = self.f.next()
        return line
    
    def __iter__(self):
        return self

def singleton(cls):
    """ Decorator function that turns a class into a singleton. """
    import inspect

    # Create a structure to store instances of any singletons that get
    # created.
    instances = {}

    # Make sure that the constructor for this class doesn't take any
    # arguments.  Since singletons can only be instantiated once, it doesn't
    # make any sense for the constructor to take arguments.
    try:
        specification = inspect.getargspec(cls.__init__)
        positional, variable, keyword, default = specification
        message = "Singleton classes cannot accept arguments to the constructor."

        # The constructor should have a self argument, but no others.
        if len(positional) is not 1:
            raise TypeError(message)
        if (variable is not None) or (keyword is not None):
            raise TypeError(message)

    # If the class doesn't have a constructor, that's ok.
    except AttributeError:
        pass

    def get_instance():
        """ Creates and returns the singleton object.  This function is what 
        gets returned by this decorator. """

        # Check to see if an instance of this class has already been
        # instantiated.  If it hasn't, create one.  The `instances` structure
        # is technically a global variable, so it will be preserved between
        # calls to this function.
        if cls not in instances:
            instances[cls] = cls()

        # Return a previously instantiated object of the requested type.
        return instances[cls]

    # Return the decorator function.
    return get_instance

def read_tdt(path):
    """
    Read tab delimited text file with a headers line. Return dictionary of
    column_name to value for each row.
    """
    with open(path) as f:
        dict_reader = csv.DictReader(f, delimiter='\t')
        for row_dict in dict_reader:
            yield row_dict

def write_tdt(path, rows, header=None):
    """Incomplete"""
    with open(path, 'wb') as f:
        writer = csv.writer(f, delimiter='\t')
        if header:
            writer.writerow(header)
        writer.writerows(rows)



def inverse_dict(d):
    """Return a new dictionary with key and values switched."""
    return {value: key for key, value in d.iteritems()}




if __name__ == '__main__':
    print current_data_dir('drugbank')
    print preceding_data_dir('drugbank')
