#! /usr/bin/python

from astropy.io import fits
import numpy as np
import sys

"""
Transform a (row, column) index to a (x, y) index, given the shape of a matrix.
Index (0,0) is at the top left.
"""
def rc2car(r,c): return c,r
def car2rc(x,y): return y,x

def loadFitsImageData(fn: str):
    """ 
    Dead-simple loading of FITS image data to a numpy 2D array

    No information is extracted; Also, we suppose we are only interested in 
    the first HDU (in case more than one exist).
    """
    with fits.open(fn) as hdul: return hdul[0].data

def getVerticalSlit(mat, row: int, col: int, width: int, height: int):
    n, m = mat.shape
    if row >= n or col >= m:
        raise RuntimeError("[ERROR] Base coordinates of slit requested are out of image (dim:{:}x{:} requested:({:},{:}))!\n".format(n,m,row,col))
    if width%2 == 0:
        print("[WRNNG] Slit width requested ({:} pixels) is even; truncating to nearest odd number (i.e. {:})".format(width, width+1), file=sys.stderr)
        width = width + 1
    w = width // 2
    if (w!=1) and (col - w < 0 or col + w >= m):
        print("[ERROR] Invalid slit width! Requested width {:} centered at {:} but matrix width is {:}".format(width, col, m), file=sys.stderr)
        raise RuntimeError("[ERROR] Invalid slit width! The slit requested would fall outside the image\n")
    h = height // 2
    if (h!=1) and (row - h < 0 or col + h >= n):
        print("[ERROR] Invalid slit height! Requested height {:} centered at {:} but matrix height is {:}".format(height, row, n), file=sys.stderr)
        raise RuntimeError("[ERROR] Invalid slit height! The slit requested would fall outside the image\n")
    l, r = (col-w, col+w+1)
    t, b = (row-h, row+h+1)
    return mat[t:b, l:r].flatten() if (width<=1 or height<=1) else mat[t:b, l:r]

def setCenterPixel(mat, row :int, col :int, rectangular=False, missing_vals=0):
    """
    Given a matrix of dimensions nxm and index of one pixel, i.e. (row, col), 
    the function will return a new matrix, of size kxl, such that the pixel 
    (row, col) is placed exactly at the center of the new matrix.

    The new (returned) matrix, will most probably be of larger dimensions than 
    the input matrix.

    Parameters:
        mat: Original 2d-array (matrix) of dimensions (n,m)
        row: Row of central pixel, in range [0,n)
        col: Column of central pixel [0,m)
        rectangular: If set to true, the resulting matrix will be rectangular
        missing_vals: Set value for elements of the resulting matrix to be added

    Returns:
        A 2D-array (matrix) of dimensions (k,l), where k>=n and l>=m, with the 
        centered around the pixel (row, col) of the original matrix.
    """
    n,m = mat.shape
    k,l = 2*max(row, n-row-1)+1, 2*max(col, m-col-1)+1
    if rectangular: k=l=max(k,l)
    zer = np.full((k, l), missing_vals)
    tl = (k//2-row, l//2-col)
    zer[tl[0]:tl[0]+n, tl[1]:tl[1]+m] = mat
    return zer
