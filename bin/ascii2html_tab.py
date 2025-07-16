#!/usr/bin/env python

from satellite import plotters

import argparse
import os
import sys
import logging

satellite_version = "2.r1"

class myFormatter(
    argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter
):
    pass


parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description="Spectroscopic Analysis Tool for intEgraL FieLd unIt daTacubEs",
    epilog=(
        """National Observatory of Athens,
    Institute for Astronomy, Astrophysics, Space Applications and Remote Sensing\n
    Send bug reports to:
    Stavros Akras, xanthos@mail.ntua.gr
    May, 2024"""
    ),
)

parser.add_argument(
    "-b",
    "--total-abundancies",
    metavar="TOTAL_ABUNDANCIES",
    dest="totalabun",
    default=None,
    required=True,
    help="Total abundancies data file.",
)

if __name__ == "__main__":

    # parse cmd
    args = parser.parse_args()
    columns, dct = plotters.parseTotalAbundancies(args.totalabun)

    nr_columns = len(columns) + 1
    offset = None
    try:
        start = int(columns[0])
        offset = 1 if (start == 0) else 0
    except:
        pass    

    print("""
            <!DOCTYPE html>
            <html lang="en">
            <head>
              <meta charset="UTF-8" />
              <title>Total Abundancies Table</title>
              <style>
              table {
                width: 50%;
                border-collapse: collapse;
                margin: 20px 0;
              }
              table, th, td {
                border: 1px solid #333;
              }
              thead th {
                background-color: #ffe6e6; /* light pink */
                text-align: center;
                font-weight: bold; /* default, but explicit */
              }
              tbody td {
                font-weight: normal; /* make sure tbody is not bold */
                text-align: left; /* or whatever you prefer */
                background-color: transparent; /* optional, explicit */
              }
              th, td {
                padding: 8px 12px;
              }
              </style>
            </head>
            <body>
    """)

    def referenceAbundance(element, slit_nr):
        for k, v in dct.items():
            if k == element:
                for ki, vi in v.items():
                    if ki == 'TA':
                        return vi[slit_nr][0]
        raise RuntimeError(f'Failed getting reference abundance (TA) for element {element} and slit {slit}')

    print('<table>')
    for k, v in dct.items():
        print(f'<thead><tr><th colspan="{nr_columns*2+1}">Element {k}</th></tr>')
        print('<tr><th>Slit Nr.</th>')
        for cname in columns:
            try:
                idx = int(cname) + offset
            except:
                idx = cname
            print(f'<th colspan="2">{idx}</th>')
        print('</tr></thead><tbody>')
        for ki, vi in v.items():
# TA or ICF
            print(f'<tr> <th>{ki}</th>')
            for ve in vi:
                print(f'<td>{ve[0]:.3e}</td><td>{ve[1]:.3e}</td>')
            print('</tr>')
# Ratio w.r.t reference abundance
            if ki != 'TA':
                print(f'<tr> <th>ICF</th>')
                for j, ve in enumerate(vi):
                    print(f'<td colspan="2">{ve[0]/referenceAbundance(k, j):.2f}</td>')
                print('</tr>')
        print('</tbody>')
    print('</table>')

    print('</body></html>')
