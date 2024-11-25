import numpy as np
import pyneb as pn

def ionAbundancies(dct, atomic_data_set):
# load PyNeb's atomic data list
    pn.atomicData.includeFitsPath()
    pn.atomicData.setDataFileDict(atomic_data_set)

# iterate atoms
    for element in dct:
        atom = element['element']
        num  = element[]
        pn.atom(atom, 
