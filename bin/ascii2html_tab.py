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

    print('<table>')
    for k, v in dct.items():
        print(f'<thead><tr><th colspan="{nr_columns*2+1}">Element {k}</th></tr>')
        print('<tr><th>Slit Nr.</th>')
        for cname in columns: print(f'<th colspan="2">{cname}</th>')
        print('</tr></thead>')
        print('<tbody>')
        for ki, vi in v.items():
            print('<tr>')
            print(f'<th>{ki}</th>')
            for ve in vi:
                print(f'<td>{ve[0]}</td><td>{ve[1]}</td>')
            print('</tr>')
        print('</tbody>')
    print('</table>')

    print('</body></html>')
