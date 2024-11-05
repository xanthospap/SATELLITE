#! /usr/bin/python3

from satellite import fitsutils
import numpy as np

if __name__ == "__main__":

# Testing setCenterPixel
    mat1 = np.array([[1,2,3],[4,5,6]])
    print(mat1)
    print('\nCenter pixel: {:},{:}={:}'.format(1,1,mat1[1][1]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 1)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][1])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(1,2,mat1[1][2]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 2)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(0,0,mat1[0][0]))
    cm   = fitsutils.setCenterPixel(mat1, 0, 0)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[0][0])
    print(cm)
    
    mat1 = np.array([[1,2,3],[4,5,6],[7,8,9]])
    print(mat1)
    print('\nCenter pixel: {:},{:}={:}'.format(1,1,mat1[1][1]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 1)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][1])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(1,2,mat1[1][2]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 2)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(0,0,mat1[0][0]))
    cm   = fitsutils.setCenterPixel(mat1, 0, 0)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[0][0])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(2,2,mat1[2][2]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 2)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][2])
    print(cm)
    
    mat1 = np.array([[1,2,3,-1],[4,5,6,-2],[7,8,9,-3]])
    print(mat1)
    print('\nCenter pixel: {:},{:}={:}'.format(1,1,mat1[1][1]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 1)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][1])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(1,2,mat1[1][2]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 2)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(0,0,mat1[0][0]))
    cm   = fitsutils.setCenterPixel(mat1, 0, 0)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[0][0])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(2,2,mat1[2][2]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 2)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:}'.format(2,3,mat1[2][3]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 3)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][3])
    print(cm)
    
# Using the rectangular option
    mat1 = np.array([[1,2,3],[4,5,6]])
    print(mat1)
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(1,1,mat1[1][1]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 1, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][1])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(1,2,mat1[1][2]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 2, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(0,0,mat1[0][0]))
    cm   = fitsutils.setCenterPixel(mat1, 0, 0, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[0][0])
    print(cm)
    
    mat1 = np.array([[1,2,3],[4,5,6],[7,8,9]])
    print(mat1)
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(1,1,mat1[1][1]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 1, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][1])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(1,2,mat1[1][2]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 2, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(0,0,mat1[0][0]))
    cm   = fitsutils.setCenterPixel(mat1, 0, 0, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[0][0])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(2,2,mat1[2][2]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 2, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][2])
    print(cm)
    
    mat1 = np.array([[1,2,3,-1],[4,5,6,-2],[7,8,9,-3]])
    print(mat1)
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(1,1,mat1[1][1]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 1, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][1])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(1,2,mat1[1][2]))
    cm   = fitsutils.setCenterPixel(mat1, 1, 2, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[1][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(0,0,mat1[0][0]))
    cm   = fitsutils.setCenterPixel(mat1, 0, 0, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[0][0])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(2,2,mat1[2][2]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 2, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][2])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(2,3,mat1[2][3]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 3, True)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][3])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(2,3,mat1[2][3]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 3, True, np.inf)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][3])
    print(cm)
    
    print('\nCenter pixel: {:},{:}={:} out->rectangular'.format(2,3,mat1[2][3]))
    cm   = fitsutils.setCenterPixel(mat1, 2, 3, True, -999)
    k,l=cm.shape
    assert(cm[k//2][l//2]==mat1[2][3])
    print(cm)
