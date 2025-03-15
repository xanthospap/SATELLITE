import re
import pyneb as pn
import numpy as np
import satellite.roman as sr


def ctrue(*args): return True


def omega(elemspec_abundancies):
    return elemspec_abundancies['O3'] / (elemspec_abundancies['O2'] + elemspec_abundancies['O3'])


def U(elemspec_abundancies):
    return elemspec_abundancies['He3'] / (elemspec_abundancies['He2'] + elemspec_abundancies['He3'])


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
    """ This function actually applies the dictionary elemspec2ionstr_dict.
        Given an element and a spectrum, we construct a (possible) key entry 
        of the dictionary, and return the corresponding value.

        Parameters
        -----------
        element: str 
            The element, e.g. 'He', 'H', 'O', etc
        spectrum:
            The spectrum, which can be either an integer value (i.e. 2) or 
            an (upper- or lower-case) roman letter; e.g. ('i', 'II', etc)

        Returns
        -------
        str:
            The corresponding entry of the elemspec2ionstr_dict dictionary.

        Examples
        --------
        >>> elemspec2ionstr('He', 1)     -> 'He+'
        >>> elemspec2ionstr('Fe', 'iii') -> 'Fe2+'
    """
    # if spectrum is an int, turn it to a roman letter
    # print('elemspec2ionstr got {:} and {:}'.format(element, spectrum.upper()))
    try:
        # print('\ttrying int spectrum ...')
        int(spectrum)
        # print('\tint ok ...')
        spectrum = sr.int2roman(int(spectrum))
        # print('\tint2roman ok ...')
    except:
        pass
    print('Note: {:}{:} -> {:}'.format(element, spectrum,
          elemspec2ionstr_dict['{:}{:}'.format(element, spectrum.upper())]))
    return elemspec2ionstr_dict['{:}{:}'.format(element, spectrum.upper())]


def ionstr2pyneb(ion):
    """ Given an ion string (e.g. 'Fe2+') the function will return the 
        corresponding pyneb literal.

        Parameters
        ----------
        ion: str
            The ion string (e.g 'Fe2+', 'He+') etc.

        Returns
        -------
        str:
            The corresponding pyneb entry/string

        Examples:
    """
    g = re.match(r'([A-Za-z]+)(0)', ion)
    try:
        return g[1]
    except:
        pass
    # print('-->matching [{:}]'.format(ion))
    g = re.match(r'([A-Za-z]+)([0-9]?)\+', ion)
    try:
        print('Note {:} -> {:}'.format(ion,
              g[1] + (str(int(g[2])+1) if g[2] != '' else '2')))
        return g[1] + (str(int(g[2])+1) if g[2] != '' else '2')
    except:
        raise RuntimeError(
            'ERROR Failed transforming Ion {:} to PyNeb string'.format(ion))


icf_dict = {
    'Ar': {'kb': [{'name': 'KB94_A32', 'cond': ctrue}, {'name': 'KB94_A30.10', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_35', 'cond': omega1435}, {'name': 'DIMS14_36', 'cond': omega1436}]},
    'N': {'kb': [{'name': 'KB94_A1.10', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_14', 'cond': omega1414}, {'name': 'DIMS14_14b', 'cond': omega1414b}]},
    'O': {'kb': [{'name': 'KB94_A10', 'cond': ctrue}, {'name': 'KB94_A8', 'cond': ctrue}, {'name': 'KB94_A6', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_12', 'cond': ctrue}]},
    'S': {'kb': [{'name': 'KB94_A36.10', 'cond': ctrue}, {'name': 'KB94_A38.10', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_23', 'cond': ctrue}, {'name': 'DIMS14_26', 'cond': ctrue}]},
    'Cl': {'kb': [], 'dims': [{'name': 'DIMS14_29b', 'cond': omega1429b}, {'name': 'DIMS14_32', 'cond': omega1432}]},
    'Ne': {'kb': [{'name': 'KB94_A28.10', 'cond': ctrue}, {'name': 'KB94_A27', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_17a', 'cond': ctrue}, {'name': 'DIMS14_17b', 'cond': ctrue}, {'name': 'DIMS14_17c', 'cond': u1417c}]},
    'C': {'kb': [{'name': 'KB94_A12', 'cond': ctrue}, {'name': 'KB94_A13.10', 'cond': ctrue}, {'name': 'KB94_A16', 'cond': ctrue}, {'name': 'KB94_A19', 'cond': ctrue}, {'name': 'KB94_A21', 'cond': ctrue}, {'name': 'KB94_A26', 'cond': ctrue}], 'dims': [{'name': 'DIMS14_39', 'cond': ctrue}]},
    'Fe': {'kb': [{'name': 'RR05_2', 'cond': ctrue}, {'name': 'RR05_3', 'cond': ctrue}, {'name': 'RR05_4', 'cond': ctrue}], 'dims': []}
}

def IcfNameList():
    """ Return a list with all values 'name' from the icf_dict dictionary.
    """
    return [entry['name'] for element in icf_dict.values() for section in element.values() for entry in section]


def renameIons(abundancies):
    acopy = {}

    def es_split(elemsp):
        g = re.match(r"([A-Za-z]+)([0-9]?)", elemsp)
        return g[1], g[2]
    for k, v in abundancies.items():
        acopy[ionstr2pyneb(elemspec2ionstr(*es_split(k)))] = v
        print('--> renaming {:} to {:}'.format(k,
              ionstr2pyneb(elemspec2ionstr(*es_split(k)))))
    return acopy


def computeIcfs(abundancies, logger):
    # Returned dictionary, in the form:
    # {
    # 'Ne': {'KB94_A28.10': val, 'DIMS14_17a': val, ...},
    # 'Fe': { ...}
    # }
    ret = {}
    # Initialize PyNeb's ICF
    icf = pn.ICF()
    allIcf = icf.getElemAbundance(renameIons({a[0]: a[1][0] for a in abundancies.items()}), IcfNameList())
    print(allIcf)  
    return allIcf

def computeIcfs_obsolete(abundancies, logger):
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
    abundancies = renameIons({a[0]: a[1][0] for a in abundancies.items()})
    print('----------+++Abundancies={:}'.format(abundancies))
    # get list of unique elements, assuming keys in abundancies are of type e.g. 'O3'
    # print(abundancies)
    # if element in list(set([re.match(r"([A-Za-z]+)([0-9]?)", x)[1] for x in abundancies])):
    for element in ['O', 'N', 'Ne', 'S', 'Cl', 'Ar', 'C']:
        if element not in ['H', 'He']:
            # get the entry of icf_dict for the element
            entry = icf_dict[element]
            icfs_element_dict = {}
            # compute abundancy from KBs
            for kb in entry['kb']:
                if kb['cond'](abundancies):
                    print(
                        "icf.getElemAbundance(abundancies, icf_list={:}".format(kb['name']))
                    icfs_element_dict.update(icf.getElemAbundance(
                        abundancies, icf_list=kb['name']))
                    # print('icf.getElemAbundance(abundancies, icf_list={:})'.format(kb['name']))
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

    # Get all outermost keys dynamically
    outer_keys = sorted(dict_of_icfs.keys())

    with open(fn, 'w') as fout:

        # Print header
        headers = ["#"] + [f"Slit {key}" for key in outer_keys]
        header_format = "{:<20}" + " {:>13}" * len(outer_keys)
        print(header_format.format(*headers), file=fout)

        # Print each row
        for key, inner_key in inner_combinations:
            row = [f"{key}_{inner_key}"]
            for outer_key in outer_keys:
                # Access value or use nan if not present
                value = dict_of_icfs[outer_key].get(
                    key, {}).get(inner_key, np.nan)
                if np.isnan(value):
                    row.append("nan")
                else:
                    row.append("{:13.12e}".format(value))

            # Print row with dynamic column count
            print(header_format.format(*row), file=fout)
