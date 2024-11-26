#! /usr/bin/python

import yaml
import os

def parseConfigInout(fn: str):
    """ 
    Parse an input SATELLITE configuration file (YAML format) and return 
    the corresponding options as a dictionary.
    """
    with open(fn, "r") as fin: 
        ymldata = fin.read()
    return yaml.safe_load(ymldata)

def configElements(dct: dict):
    """
    Once we have parsed the YAML config file into a dictionary, we can use 
    the following function to get a list of all the elements that are included.

    Results should be something like:
    ['H', 'He', 'N', 'O']
    """
    return [ k['atom'] for k in dct['data_list']['element_list'] ]

def configFitsFileList(dct: dict):
    """
    Once we have parsed the YAML config file into a dictionary, we can use
    the following function to get a list of all the FITS images and related info 
    that should be present.

    Results should be something like: 
    [
      {'type': 'data', 'element': 'H', 'spectrum' : 'i',  'atomic': 6563, 'fn': '/home/.../Hi_6563s.fit'}
      {'type': 'data', 'element': 'H', 'spectrum' : 'i',  'atomic': 4861, 'fn': '/home/.../Hi_4861s.fit'}
      {'type': 'data', 'element': 'H', 'spectrum' : 'i',  'atomic': 4340, 'fn': '/home/.../Hi_4340s.fit'}
      {'type': 'data', 'element': 'H', 'spectrum' : 'i',  'atomic': 4101, 'fn': '/home/.../Hi_4101s.fit'}
      {'type': 'data', 'element': 'He', 'spectrum': 'i',  'atomic': 5876, 'fn': '/home/.../Hei_5876s.fit'}
      {'type': 'error', 'element': 'He', 'spectrum': 'i',  'atomic': 6678, 'fn': '/home/.../Hei_6678e.fit'}
      {'type': 'error', 'element': 'He', 'spectrum': 'ii', 'atomic': 5412, 'fn': '/home/.../Heii_5412e.fit'}
    ]
    """
    fns = []
    did = dct['data_list']['data_img_id']
    eid = dct['data_list']['error_img_id']
    for element in dct['data_list']['element_list']: 
        for spec in element['spectrum']: 
            for atomic in element['spectrum'][spec]['atomic_list']:
                fn = element['atom'] + spec + '_' + str(atomic) + did +'.' + dct['data_list']['suffix']
                fns.append({'type': 'data', 'element': element['atom'], 'spectrum': spec, 'atomic': atomic, 'fn': os.path.join(dct['data_list']['prefix'], fn)})
                fn = element['atom'] + spec + '_' + str(atomic) + eid +'.' + dct['data_list']['suffix']
                fns.append({'type': 'data', 'element': element['atom'], 'spectrum': spec, 'atomic': atomic, 'fn': os.path.join(dct['data_list']['prefix'], fn)})
    return fns

def configSpecificSlitAnalysis(dct: dict):
    d = dct['analysis']['specific_slit_analysis']
# no specific-slit analysis
    if d['skip'] in ['1', 'True', 'true']: return None
# return a list of dictionaries, one for each slit
    slits = []
    for slit in d['slits']:
        ps = [ int(x) for x in slit.split(',') ]
        slits.append({'PA': ps[0], 'w': ps[1], 'h': ps[2], 'x': ps[3], 'y': ps[4]})
    return slits
    
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
