#! /usr/bin/python3

from satellite import fitsutils
import numpy as np
import sys

if __name__ == "__main__":

# Testing setCenterPixel
    mat1 = np.array([[1,2,3],[4,5,6]])
    print(mat1)

    cm = fitsutils.getVerticalSlit(mat1, 0, 1, 1, 1)
    print(cm)
    assert(cm == np.array([2]))
    
    cm = fitsutils.getVerticalSlit(mat1, 1, 2, 1, 1)
    print(cm)
    assert(cm == np.array([6]))
    
    try:
        cm = fitsutils.getVerticalSlit(mat1, 2, 2, 1, 1)
        print("Expected exception but function did not throw!", file=sys.stderr)
        sys.exit(9)
    except:
        pass
    
    mat1 = np.array([[1,2,3,4,5,6],[4,5,6,7,8,9],[0,9,8,7,6,5],[0,1,2,3,4,5]])
    print(mat1)

    cm = fitsutils.getVerticalSlit(mat1, 0, 1, 1, 1)
    print(cm)
    assert(cm == np.array([2]))
    
    cm = fitsutils.getVerticalSlit(mat1, 1, 2, 1, 1)
    print(cm)
    assert(cm == np.array([6]))
    
    cm = fitsutils.getVerticalSlit(mat1, 0, 4, 1, 4)
    print(cm)
    assert(np.array_equal(np.array([5,8,6,4]), cm))
    
    cm = fitsutils.getVerticalSlit(mat1, 0, 4, 3, 4)
    print(cm)
    #assert(np.array_equal(np.array([]), cm))
    
    cm = fitsutils.getVerticalSlit(mat1, 0, 4, 2, 3)
    print(cm)
    
    try:
        cm = fitsutils.getVerticalSlit(mat1, 0,4,1,5)
        print("Expected exception but function did not throw!", file=sys.stderr)
        sys.exit(9)
    except:
        pass
