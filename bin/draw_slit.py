#!/usr/bin/env python

from satellite import version
import satellite.fitsutils as fs
import os
import argparse
import re
from collections import defaultdict
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

CORNER_ORDER = ["TL", "TR", "BR", "BL"]

LINE_RE = re.compile(
    r"^\s*(\d+)\s+([A-Za-z]{2})\s+([+-]?\d+(?:\.\d*)?)\s*,\s*([+-]?\d+(?:\.\d*)?)\s*$"
)


def parse_corners_file(path, one_indexed=False):
    """
    Returns dict: slit_id -> dict(corner_name -> (row_float, col_float))
    """
    slits = defaultdict(dict)
    with open(path, "r") as f:
        for lineno, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = LINE_RE.match(line)
            if not m:
                raise ValueError(f"Parse error on line {lineno}: {line!r}")
            slit_id = int(m.group(1))
            corner = m.group(2).upper()
            if corner not in CORNER_ORDER:
                raise ValueError(f"Unknown corner {corner!r} on line {lineno}")
            row = float(m.group(3))
            col = float(m.group(4))
            if one_indexed:
                row -= 1.0
                col -= 1.0
            slits[slit_id][corner] = (row, col)
    # sanity: each slit should have all 4 corners
    for sid, corners in slits.items():
        missing = [c for c in CORNER_ORDER if c not in corners]
        if missing:
            raise ValueError(f"Slit {sid} missing corners: {missing}")
    return slits


def bresenham_line(r0, c0, r1, c1):
    """
    Integer line rasterization (inclusive endpoints).
    Returns arrays (rows, cols).
    """
    r0, c0, r1, c1 = int(r0), int(c0), int(r1), int(c1)
    dr = abs(r1 - r0)
    dc = abs(c1 - c0)
    sr = 1 if r0 < r1 else -1
    sc = 1 if c0 < c1 else -1
    r, c = r0, c0
    rows, cols = [], []
    if dc > dr:
        err = dc // 2
        while c != c1:
            rows.append(r)
            cols.append(c)
            err -= dr
            if err < 0:
                r += sr
                err += dc
            c += sc
        rows.append(r1)
        cols.append(c1)
    else:
        err = dr // 2
        while r != r1:
            rows.append(r)
            cols.append(c)
            err -= dc
            if err < 0:
                c += sc
                err += dr
            r += sr
        rows.append(r1)
        cols.append(c1)
    return np.asarray(rows), np.asarray(cols)


def draw_poly_outline(data, corners_rc, value, width=1):
    """
    Draw 1+ pixel wide outline connecting corners in order and closing the loop.
    corners_rc: list/array of (row, col) (floats ok; will be rounded).
    """
    n, m = data.shape
    pts = np.asarray(corners_rc, dtype=float).round().astype(int)
    pts[:, 0] = np.clip(pts[:, 0], 0, n - 1)
    pts[:, 1] = np.clip(pts[:, 1], 0, m - 1)

    edges = list(zip(pts, np.vstack([pts[1:], pts[:1]])))
    # offsets for thickness
    t = max(0, int((width - 1) // 2))
    offsets = [(dr, dc) for dr in range(-t, t + 1) for dc in range(-t, t + 1)]

    for (r0, c0), (r1, c1) in edges:
        rr, cc = bresenham_line(r0, c0, r1, c1)
        for dr, dc in offsets:
            rrt = rr + dr
            cct = cc + dc
            mask = (rrt >= 0) & (rrt < n) & (cct >= 0) & (cct < m)
            data[rrt[mask], cct[mask]] = value


def pick_outline_value(data):
    if np.issubdtype(data.dtype, np.integer):
        return np.iinfo(data.dtype).max
    # floating
    vmax = np.nanmax(data)
    if not np.isfinite(vmax):
        return 1.0
    return vmax


def mark_slits_on_fits(
    corners_txt,
    fits_in,
    fits_out="koko.fits",
    one_indexed=False,
    outline_width=1,
    outline_value=None,
):
    slits = parse_corners_file(corners_txt, one_indexed=one_indexed)

    with fits.open(fits_in, mode="readonly") as hdul:
        # work on a copy of primary image
        data = hdul[0].data
        if data is None:
            raise RuntimeError("Primary HDU has no image data.")
        data_out = data.copy()
        val = (
            outline_value if outline_value is not None else pick_outline_value(data_out)
        )

        for sid in sorted(slits.keys()):
            corners = [slits[sid][k] for k in CORNER_ORDER]  # TL, TR, BR, BL
            draw_poly_outline(data_out, corners, value=val, width=outline_width)

        # write out: preserve header, add HISTORY
        hdu = fits.PrimaryHDU(data_out, header=hdul[0].header)
        hdu.header.add_history("Slit outlines drawn from file: {}".format(corners_txt))
        hdul_out = fits.HDUList([hdu])
        hdul_out.writeto(fits_out, overwrite=True)


def save_preview_png(corners_txt, fits_in, out_png="koko.png", one_indexed=False):
    """
    Make a quicklook image with slit outlines + indices.
    - Draws polygons from TL,TR,BR,BL points
    - Labels each slit at its centroid with the slit index
    """
    slits = parse_corners_file(corners_txt, one_indexed=one_indexed)

    with fits.open(fits_in) as hdul:
        img = hdul[0].data
        if img is None:
            raise RuntimeError("Primary HDU has no image data.")
        img = img.astype(float)

    # simple contrast stretch
    vmin, vmax = np.nanpercentile(img, [0.5, 99.5])

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(img, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)

    for sid in sorted(slits.keys()):
        # corners in (row,col) -> matplotlib needs (x=col, y=row)
        corners_rc = np.array([slits[sid][k] for k in CORNER_ORDER], dtype=float)
        corners_xy = np.column_stack([corners_rc[:, 1], corners_rc[:, 0]])

        poly = Polygon(corners_xy, fill=False, lw=1.5, ec="yellow")
        ax.add_patch(poly)

        center_rc = corners_rc.mean(axis=0)
        ax.text(
            center_rc[1],
            center_rc[0],
            str(sid),
            color="yellow",
            fontsize=10,
            ha="center",
            va="center",
            bbox=dict(boxstyle="round,pad=0.2", fc="black", ec="none", alpha=0.5),
        )

    ax.set_xlabel("Col")
    ax.set_ylabel("Row")
    ax.set_title(os.path.basename(fits_in))
    plt.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="Mark slit outlines on a FITS image.")
    p.add_argument("corners_txt", help="Text file with 'Slit Corner Row, Col' lines")
    p.add_argument("fits_in", help="Input FITS image")
    p.add_argument("-o", "--output", default="koko.fits", help="Output FITS filename")
    p.add_argument(
        "--one-indexed",
        action="store_true",
        help="Interpret Row,Col as 1-indexed (FITS-style).",
    )
    p.add_argument("--width", type=int, default=1, help="Outline width in pixels")
    p.add_argument(
        "--value",
        type=float,
        default=None,
        help="Pixel value to use for the outline (default: dtype max for ints, max(data) for floats)",
    )
    args = p.parse_args()

    mark_slits_on_fits(
        args.corners_txt,
        args.fits_in,
        fits_out=args.output,
        one_indexed=args.one_indexed,
        outline_width=args.width,
        outline_value=args.value,
    )

    save_preview_png(
        args.corners_txt, args.fits_in, out_png="koko.png", one_indexed=args.one_indexed
    )


if __name__ == "__main__":
    main()
