#! /usr/bin/python

from astropy.io import fits
import numpy as np
from scipy.ndimage import rotate, affine_transform
import sys

"""
Transform a (row, column) index to a (x, y) index, given the shape of a matrix.
Index (0,0) is at the top left.
"""


def rc2car(r, c):
    return c, r


def car2rc(x, y):
    return y, x


def loadFitsImageData(fn: str):
    """
    Dead-simple loading of FITS image data to a numpy 2D array

    No information is extracted; Also, we suppose we are only interested in
    the first HDU (in case more than one exist).
    """
    with fits.open(fn) as hdul:
        return hdul[0].data


def getVerticalSlit(mat, row: int, col: int, width: int, height: int, logger):
    """WARNING
    Keep this function consistent with
    def _slit_box_in_rotated_rc(row, col, width, height)
    defined below
    """
    n, m = mat.shape
    if row >= n or col >= m:
        raise RuntimeError(
            "[ERROR] Base coordinates of slit requested are out of image (dim:{:}x{:} requested:({:},{:}))!\n".format(
                n, m, row, col
            )
        )
    if width % 2 == 0:
        if logger:
            logger.warning(
                "Slit width requested ({:} pixels) is even; truncating to nearest odd number (i.e. {:})".format(
                    width, width + 1
                )
            )
        width = width + 1
    w = width // 2
    if (w != 1) and (col - w < 0 or col + w >= m):
        if logger:
            logger.error(
                "Invalid slit width! Requested width {:} centered at {:} but matrix width is {:}".format(
                    width, col, m
                )
            )
        raise RuntimeError(
            "[ERROR] Invalid slit width! The slit requested would fall outside the image\n"
        )
    h = height // 2
    if (h != 1) and (row - h < 0 or row + h >= n):
        if logger:
            logger.error(
                "Invalid slit height! Requested height {:} centered at {:} but matrix height is {:}".format(
                    height, row, n
                )
            )
        raise RuntimeError(
            "[ERROR] Invalid slit height! The slit requested would fall outside the image\n"
        )
    l, r = (col - w, col + w + 1)
    t, b = (row - h, row + h + 1)
    return mat[t:b, l:r].flatten() if (width <= 1 or height <= 1) else mat[t:b, l:r]


def getVerticalSlitPolygon(mat, row: int, col: int, width: int, height: int, logger):
    n, m = mat.shape
    if row >= n or col >= m:
        raise RuntimeError(
            "[ERROR] Base coordinates of slit requested are out of image (dim:{:}x{:} requested:({:},{:}))!\n".format(
                n, m, row, col
            )
        )
    if width % 2 == 0:
        if logger:
            logger.warning(
                "Slit width requested ({:} pixels) is even; truncating to nearest odd number (i.e. {:})".format(
                    width, width + 1
                )
            )
        width = width + 1
    w = width // 2
    if (w != 1) and (col - w < 0 or col + w >= m):
        if logger:
            logger.error(
                "Invalid slit width! Requested width {:} centered at {:} but matrix width is {:}".format(
                    width, col, m
                )
            )
        raise RuntimeError(
            "[ERROR] Invalid slit width! The slit requested would fall outside the image\n"
        )
    h = height // 2
    if (h != 1) and (row - h < 0 or col + h >= n):
        if logger:
            logger.error(
                "Invalid slit height! Requested height {:} centered at {:} but matrix height is {:}".format(
                    height, row, n
                )
            )
        raise RuntimeError(
            "[ERROR] Invalid slit height! The slit requested would fall outside the image\n"
        )
    l, r = (col - w, col + w + 1)
    t, b = (row - h, row + h + 1)
    # return mat[t:b, l:r].flatten() if (width <= 1 or height <= 1) else mat[t:b, l:r]
    return (row, col), (t, l), (t, r - 1), (b - 1, l), (b - 1, r - 1)


def setCenterPixel(mat, row: int, col: int, rectangular=False, missing_vals=0):
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
    n, m = mat.shape
    k, l = 2 * max(row, n - row - 1) + 1, 2 * max(col, m - col - 1) + 1
    if rectangular:
        k = l = max(k, l)
    zer = np.full((k, l), missing_vals)
    tl = (k // 2 - row, l // 2 - col)
    zer[tl[0] : tl[0] + n, tl[1] : tl[1] + m] = mat
    return zer


#  rotate around (0,0) in counterclockwise fashion by angle angle_deg
def rotate2d_00(array, angle_deg):
    return rotate(
        array, angle=angle_deg, reshape=True
    )  # reshape=True preserves full content


def _rotation_matrix(angle_deg):
    a = np.deg2rad(angle_deg)
    ca, sa = np.cos(a), np.sin(a)
    # Row/col convention: vectors are [row, col]
    return np.array([[ca, -sa], [sa, ca]], dtype=float)


#  rotate around (centerx, centery) in counterclockwise fashion by angle angle_deg
def rotate2d(array, angle_deg, centerx=0, centery=0):
    if centerx == centery and centerx == 0:
        return rotate2d_00(array, -1.0 * angle_deg)

    # angle_rad = np.deg2rad(angle_deg)
    # cos_a = np.cos(angle_rad)
    # sin_a = np.sin(angle_rad)

    # Rotation matrix
    # rotation_matrix = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    R = _rotation_matrix(angle_deg)

    # Compute the offset to keep the rotation centered at `center`
    center = np.asarray([centerx, centery])
    offset = center - R @ center

    # Perform affine transformation
    rotated = affine_transform(
        array,
        R,
        offset=offset,
        order=1,  # linear interpolation
        mode="constant",
        cval=0.0,
    )
    return rotated


def _slit_box_in_rotated_rc(row, col, width, height):
    """
    Returns the four corners (row, col) of the slit rectangle
    in the *rotated-image* coordinate frame (pixel centers).
    Corner order: TL, TR, BR, BL.

    WARNING!
    --------
    This should do exactly what getVerticalSlit does, but not operate on a
    matrix.
    """
    # enforce odd width like your getVerticalSlit
    if width % 2 == 0:
        width += 1
    w = width // 2
    h = height // 2

    # getVerticalSlit uses [t:b, l:r] with r/b exclusive
    t, b = row - h, row + h + 1
    l, r = col - w, col + w + 1

    # Corners at pixel centers (inclusive indices use -1 on the r/b side)
    TL = np.array([t, l], dtype=float)
    TR = np.array([t, r - 1], dtype=float)
    BR = np.array([b - 1, r - 1], dtype=float)
    BL = np.array([b - 1, l], dtype=float)
    return np.vstack([TL, TR, BR, BL])


def corners_in_original_from_rotated(corners_rc_rot, center_rc, angle_used_deg):
    """
    Map points from the *rotated image* back to *original image* pixel coords.
    If the rotated image was produced by rotating ORIGINAL by angle_used_deg
    about center_rc, the inverse mapping uses -angle_used_deg.
    """
    C = np.asarray(center_rc, dtype=float)
    R_inv = _rotation_matrix(-angle_used_deg)  # inverse rotation
    # (points - C) @ R_inv^T  + C   since our points are row vectors
    return (corners_rc_rot - C) @ R_inv.T + C


def testGetVerticalSlit(fits_fn, slit):
    def getFitsSlit(fits_fn, slit):
        mat = rotate2d(
            loadFitsImageData(fits_fn), slit["PA"], slit["y"] - 1, slit["x"] - 1
        )
        c, tl, tr, bl, br = getVerticalSlitPolygon(
            mat, slit["y"] - 1, slit["x"] - 1, slit["w"], slit["h"], None
        )
        return c, tl, tr, bl, br, mat

    c, tl, tr, bl, br, mat = getFitsSlit(fits_fn, slit)
    mark_value = np.max(mat) * 1.2
    for p in [c, tl, tr, bl, br]:
        mat[p[0], p[1]] = mark_value

    # Create a Primary HDU (Header/Data Unit)
    hdu = fits.PrimaryHDU(mat)

    # Create an HDU list and write to a new FITS file
    hdu.writeto("test_rotated.fits", overwrite=True)

    # rotate back
    mat = rotate2d(mat, -slit["PA"], slit["y"] - 1, slit["x"] - 1)
    hdu = fits.PrimaryHDU(mat)
    hdu.writeto("test_rotated_and_reset.fits", overwrite=True)
