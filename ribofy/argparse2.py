import argparse


"""
    Inhereted argparse class with few additional functionalities
"""
class argparse2 (argparse.ArgumentParser):
    def __init__ (self, help, *args, **kwargs):        
        super(argparse2, self).__init__(*args, **kwargs)
        self.help=help

    def format_help(self):        
        return (self.help + "\n\n")


