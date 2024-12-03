#! /usr/bin/python

import yaml
import os
import re
import sys

import pyneb as pn
from satellite import roman


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
    return [ k['element'] for k in dct['data_list']['element_list'] ]

def configFitsFileList(dct: dict):
    """
    Once we have parsed the YAML config file into a dictionary, we can use
    the following function to get a list of all the FITS images and related info 
    that should be present.

    Results should be something like: 
    [
      {'element': 'H', 'spectrum' : 'i',  'atomic': 6563, 'fns': '/home/.../Hi_6563s.fit', 'fne': '/home/.../Hi_6563e.fit'}
      {'element': 'H', 'spectrum' : 'i',  'atomic': 4861, 'fns': '/home/.../Hi_4861s.fit', 'fne': '/home/.../Hi_4861e.fit'}
      {'element': 'H', 'spectrum' : 'i',  'atomic': 4340, 'fns': '/home/.../Hi_4340s.fit', 'fne': '/home/.../Hi_4340e.fit'}
      {'element': 'H', 'spectrum' : 'i',  'atomic': 4101, 'fns': '/home/.../Hi_4101s.fit', 'fne': '/home/.../Hi_4101e.fit'}
      {'element': 'He', 'spectrum': 'i',  'atomic': 5876, 'fns': '/home/.../Hei_5876s.fit','fne': '/home/.../Hei_5876e.fit'}
    ]
    """
    fns = []
    did = dct['data_list']['data_img_id']
    eid = dct['data_list']['error_img_id']
    for element in dct['data_list']['element_list']: 
        for spec in element['spectrum']: 
            for atomic in element['spectrum'][spec]['atomic_list']:
                sfn = element['atom'] + spec + '_' + str(atomic) + did +'.' + dct['data_list']['suffix']
                efn = element['atom'] + spec + '_' + str(atomic) + eid +'.' + dct['data_list']['suffix']
                fns.append({'element': element['atom'], 'spectrum': spec, 'atomic': atomic, 'fns': os.path.join(dct['data_list']['prefix'], sfn), 'fne': os.path.join(dct['data_list']['prefix'], efn)})
    return fns

def elementInputDict(fitsd: list, reference_element: dict):
    for obj in fitsd:
        matched = True
        for key in ['element', 'spectrum', 'atomic']:
            if obj[key] != reference_element[key]: matched = False
        if matched == True: return obj
    raise RuntimeError("[ERROR] Failed matching element {:}/{:}/{:} to available input file".format(reference_element['element'], reference_element['spectrum'], reference_element['atomic']))

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

def checkInputFits(fitsd: list):
    missing_files = []
    for idx, fits in enumerate(fitsd):
        for ftype in ['fns', 'fne']:
            file_is_missing = False
            fitsfn = fits[ftype]
            if not os.path.isfile(fitsfn):
                print("[WRNNG] Missing Fits file {:}".format(fitsfn), file=sys.stderr)
                bn = os.path.basename(fitsfn)
                file_is_missing = True
# find spectrum in filename
                sstr = fits['spectrum']
# replace roman spectrum with int and see if file exists
                if bn.find(sstr) >= 0:
                    gfn = bn.replace(sstr, str(roman.roman2int(sstr)), 1)
                    guess = os.path.join(os.path.dirname(fitsfn), gfn)
                    if os.path.isfile(guess):
# change the name in the return list
                        print("[DEBUG] Fits filename {:} is missing; using {:}".format(os.path.basename(fitsfn), gfn), file=sys.stderr)
                        fitsd[idx][ftype] = guess
                        file_is_missing = False
                if file_is_missing and bn.find(sstr) >= 0:
                    gfn = bn.replace(sstr, sstr.upper(), 1)
                    guess = os.path.join(os.path.dirname(fitsfn), gfn)
                    if os.path.isfile(guess):
                        print("[DEBUG] Fits filename {:} is missing; using {:}".format(os.path.basename(fitsfn), gfn), file=sys.stderr)
                        fitsd[idx][ftype] = guess
                        file_is_missing = False
            if file_is_missing: missing_files.append(fitsfn)
    return fitsd, missing_files

def indexOf(atom: str, spectrum: str, atomic_number: int, fitsd: list):
    for idx, lst in enumerate(fitsd):
        if lst['element'] == atom and lst['spectrum'] == spectrum and lst['atomic'] == atomic_number:
            return idx
    return -1

def closestPyNebElement(key: str, atomic_number: int):
    min_diff = 1000000
    pn_wl = None
    for wl in pn.LINE_LABEL_LIST[key]:
        if re.fullmatch('[0-9]*A', wl):
            diff = abs(int(wl[0:-1]) - atomic_number)
            if diff < min_diff: 
                pn_wl = wl
                min_diff = diff
    print("[WRNNG] PyNeb is missing wavelength {:} for element {:}; using {:} instead".format(atomic_number, key, pn_wl))
    return pn_wl

def objectIntensityPyNebCode(atom: str, spectrum: str, atomic_number: int):
#
# A full list can be obtained as:
# import pyneb as pn
# pn.LINE_LABEL_LIST
#
    try:
        spectrum = int(spectrum)
    except:
        spectrum = roman.roman2int(spectrum)
    if atom in ['H', 'He']:
        pnatom = "{:}{:}r".format(atom, spectrum)
    else:
        pnatom = "{:}{:}".format(atom, spectrum)
    
    try:
        wls = pn.LINE_LABEL_LIST[pnatom]
    except:
        wls = None
        print('[ERROR] Failed matching object {:}/{:} (aka {:}) to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, pnatom), file=sys.stderr)
    if wls is None: 
        raise RuntimeError('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))

    pnatomic = '{:}A'.format(atomic_number)
    if pnatomic not in wls:
         pnatomic = closestPyNebElement(pnatom, atomic_number)

    pnstr = '_'.join([pnatom, pnatomic])

# Validate
    try:
        pn.LINE_LABEL_LIST[pnstr.split('_')[0]].index(pnstr.split('_')[1])
        return pnstr
    except:
        print('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number), file=sys.stderr)
    raise RuntimeError('[ERROR] Failed matching object {:}/{:}/{:} to PyNeb (see LINE_LABEL_LIST)'.format(atom, spectrum, atomic_number))

def satellite_str2pyneb_str(sstr: str):
    parts = sstr.split('_')
    assert(len(parts) == 2)
    atomic_number = int(parts[1])
    g = re.fullmatch("([A-Za-z]*)([0-9])", parts[0].strip())
    if g: 
        return objectIntensityPyNebCode(g[1], g[2], atomic_number)
    g = re.fullmatch("([A-Za-z]*)([iIvV]*)(r+)", parts[0].strip())
    if g: 
        return objectIntensityPyNebCode(g[1], roman.roman2int(g[2]), atomic_number)
    raise RuntimeError("[ERROR] Failed resolving satellite string {:}".format(sstr))

def partialResolve(pstr: str):
    lstpl = []; lst = [];
    for i in pstr.split('+'):
        lstpl += [1e0, i.strip()]
    for idx, s in enumerate(lstpl):
        try:
            1 + s
            lst += [s]
        except:
            if '-' in s:
                lst += [s.split('-')[0]]
                for ss in s.split('-')[1:]:
                    lst += [-1,ss]
            else:
                lst += [s]
    return lst

def resolveRatioStr(rstr: str):
    [nom, denom] = rstr.split('/')
    nom = nom.lstrip('(').rstrip(')')
    denom = denom.lstrip('(').rstrip(')')
    return partialResolve(nom), partialResolve(denom)

def configElementRatiosList(dct: dict):
    return dct['analysis']['specific_slit_analysis']['log_ratios']
