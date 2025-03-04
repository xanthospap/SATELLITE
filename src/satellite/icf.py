import re
import pyneb as pn
import numpy as np
import satellite.roman as sr


def ctrue(*args): return True


def omega(elemspec_abundancies):
    return elemspec_abundancies['O3'] / (elemspec_abundancies['O2'] + elemspec_abundancies['O3'])


def U(elemspec_abundancies):
    return elemspec_abundancies['He2'] / (elemspec_abundancies['He1'] + elemspec_abundancies['He2'])


def omega1435(ea): return omega(ea) <= .5
def omega1436(ea): return omega(ea) > .5
def omega1414(ea): return omega(ea) <= .95
def omega1414b(ea): return omega(ea) > .95
def omega1429b(ea): return omega(ea) > .02 and omega(ea) < .95
def omega1432(ea): return omega(ea) > .95
def u1417c(ea): return U(ea) < .015


elemspec2ionstr_dict = {
    'HeI': 'He+',
    'HeII': 'He2+',
    'HI': 'H+',
    'NI': 'N0',
    'NII': 'N+',
    'OI': 'O0',
    'OII': 'O+',
    'OIII': 'O2+',
    'SII': 'S+',
    'SIII': 'S2+',
    'ArIII': 'Ar2+',
    'ArIV': 'Ar3+',
    'ClII': 'Cl+',
    'ClIII': 'Cl2+',
    'CI': 'C0',
    'CII': 'C+',
    'NeIII': 'Ne2+',
    'NiII': 'Ni+',
    'CaII': 'Ca+',
    'FeII': 'Fe+',
    'FeIII': 'Fe2+',
}


def elemspec2ionstr(element, spectrum):
    # if spectrum is an int, turn it to a roman letter
    try:
        int(spectrum)
        spectrum = sr.int2roman(spectrum)
    except:
        pass
    return elemspec2ionstr_dict['{:}{:}'.format(element, spectrum.upper())]


icf_dict = {
    'Ar': {'kb': [{'name': 'KB94_A32', 'cond': ctrue}, {'name': 'KB94_A30.10', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_35', 'cond': omega1435}, {'name': 'DIMS14_36', 'cond': omega1436}]},
    'N': {'kb': [{'name': 'KB94_A1.10', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_14', 'cond': omega1414}, {'name': 'DIMS14_14b', 'cond': omega1414b}]},
    'O': {'kb': [{'name': 'KB94_A10', 'cond': ctrue}, {'name': 'KB94_A8', 'cond': ctrue}, {'name': 'KB94_A6', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_12', 'cond': ctrue}]},
    'S': {'kb': [{'name': 'KB94_A36.10', 'cond': ctrue}, {'name': 'KB94_A38.10', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_23', 'cond': ctrue}, {'name': 'DIMS14_26', 'cond': ctrue}]},
    'Cl': {'kb': [], 'dims': [{'name': 'DIMS14_29b', 'cond': omega1429b}, {'name': 'DIMS14_32b', 'cond': omega1432}]},
    'Ne': {'kb': [{'name': 'KB94_A28.10', 'cond': ctrue}, {'name': 'KB94_A27', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_17a', 'cond': ctrue}, {'name': 'DIMS14_17b', 'cond': ctrue}, {'name': 'DIMS14_17c', 'cond': u1417c}]},
    'C': {'kb': [{'name': 'KB94_A12', 'cond': ctrue}, {'name': 'KB94_A13.10', 'cond': ctrue}, {'name': 'KB94_A16', 'cond': ctrue}, {'name': 'KB94_A19', 'cond': ctrue}, {'name': 'KB94_A21', 'cond': ctrue}, {'name': 'KB94_A26', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_39', 'cond': ctrue}]},
    'Fe': {'kb': [{'name': 'RR052', 'cond': ctrue}, {'name': 'RR05_3', 'cond': ctrue}, {'name': 'RR05_4', 'cond': ctrue}], 'dims': []}
}


def computeIcfs(abundancies, logger):
    # Returned dictionary, in the form:
    # {
    # 'Ne': {'KB94_A28.10': val, 'DIMS14_17a': val, ...},
    # 'Fe': { ...}
    # }
    ret = {}
# Initialize PyNeb's ICF
    icf = pn.ICF()
    # abundancies holds values pairs of values, like 'O2': (1e-2, 1e-3)
    # Reduce it only hold abundancies values (i.e. drop the errors)
    abundancies = {a[0]: a[1][0] for a in abundancies.items()}
    # get list of unique elements, assuming keys in abundancies are of type e.g. 'O3'
    for element in list(set([re.match(r"([A-Za-z]+)([0-9]?)", x)[1] for x in abundancies])):
        if element not in ['H', 'He']:
            # get the entry of icf_dict for the element
            entry = icf_dict[element]
            icfs_element_dict = {}
            # compute abundancy from KBs
            for kb in entry['kb']:
                if kb['cond'](abundancies):
                    # print(abundancies,
                    icfs_element_dict.update(icf.getElemAbundance(
                        abundancies, icf_list=kb['name']))
            for kb in entry['dims']:
                if kb['cond'](abundancies):
                    icfs_element_dict.update(icf.getElemAbundance(
                        abundancies, icf_list=kb['name']))
            ret[element] = icfs_element_dict
    return ret

def printIcfs(dict_of_icfs, fn, logger):
    # Extract all unique combinations of inner keys
    inner_combinations = set()
    for outer_dict in dict_of_icfs.values():
        for key, sub_dict in outer_dict.items():
            for inner_key in sub_dict.keys():
                inner_combinations.add((key, inner_key))
    
    # Sort the combinations for a consistent output
    inner_combinations = sorted(inner_combinations)
    
    with open(fn, 'w') as fout:
        # Print header
        headers = ["#Slit "] + list(dict_of_icfs.keys())
        print("{:<20} {:>13} {:>13} {:>13} {:>13}".format(*headers), file=fout)
        
        # Print each row
        for key, inner_key in inner_combinations:
            row = [f"{key}_{inner_key}"]
            for outer_key in dict_of_icfs.keys():
                value = dict_of_icfs[outer_key].get(key, {}).get(inner_key, np.nan)
                if np.isnan(value):
                    row.append("nan")
                else:
                    row.append("{:13.12e}".format(value))
            print("{:<20} {:>13} {:>13} {:>13} {:>13}".format(*row), file=fout)
    return fn
