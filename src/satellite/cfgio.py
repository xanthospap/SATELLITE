#! /usr/bin/python

import yaml
import os

""" 
Parse an input SATELLITE configuration file (YAML format) and return 
the corresponding options as a dictionary.
"""
def parseConfigInout(fn: str):
    with open(fn, "r") as fin: 
        ymldata = fin.read()
    return yaml.safe_load(ymldata)

"""
Once we have parsed the YAML config file into a dictionary, we can use 
the following function to get a list of all the elements that are included.

Results should be something like:
['H', 'He', 'N', 'O']
"""
def configElements(dct: dict):
    return [ k['atom'] for k in dct['data_list']['element_list'] ]

"""
Once we have parsed the YAML config file into a dictionary, we can use
the following function to get a list of all the FITS images and related info 
that should be present.

Results should be something like: 
[
  {'element': 'H', 'spectrum' : 'i',  'atomic': 6563, 'fn': '/home/.../Hi_6563.fit'}
  {'element': 'H', 'spectrum' : 'i',  'atomic': 4861, 'fn': '/home/.../Hi_4861.fit'}
  {'element': 'H', 'spectrum' : 'i',  'atomic': 4340, 'fn': '/home/.../Hi_4340.fit'}
  {'element': 'H', 'spectrum' : 'i',  'atomic': 4101, 'fn': '/home/.../Hi_4101.fit'}
  {'element': 'He', 'spectrum': 'i',  'atomic': 5876, 'fn': '/home/.../Hei_5876.fit'}
  {'element': 'He', 'spectrum': 'i',  'atomic': 6678, 'fn': '/home/.../Hei_6678.fit'}
  {'element': 'He', 'spectrum': 'ii', 'atomic': 5412, 'fn': '/home/.../Heii_5412.fit'}
]
"""
def configFitsFileList(dct: dict):
    fns = []
    for element in dct['data_list']['element_list']: 
        for spec in element['spectrum']: 
            for atomic in element['spectrum'][spec]['atomic_list']:
                fn = element['atom'] + spec + '_' + str(atomic) + '.' + d['data_list']['suffix']
                fns.append({'element': element['atom'], 'spectrum': spec, 'atomic': atomic, 'fn': os.path.join(d['data_list']['prefix'], fn)})
    return fns
    
if __name__ == "__main__" :
    import sys
    d = parseConfigInout(sys.argv[1])
    print(d)
    print("Element dictionary is:")
    print(d['data_list']['element_list'])
    print("Element:")
    for element in d['data_list']['element_list']: 
        print(element['atom'])
        print('\tSpectrums:')
        for spec in element['spectrum']: 
            print("\t",spec)
            print("\t\tAtomics:")
            for atomic in element['spectrum'][spec]['atomic_list']:
                fn = element['atom'] + spec + '_' + str(atomic) + '.' + d['data_list']['suffix']
                print("\t\t{:} - expected fn: {:}".format(atomic, os.path.join(d['data_list']['prefix'], fn)))
    print("------------------------------------------------------------------")
    fns = configFitsFileList(d)
    for i in fns: print(i)
