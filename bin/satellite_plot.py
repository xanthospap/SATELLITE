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
    "-a",
    "--line-abundancies",
    metavar="LINE_ABUNDANCIES",
    dest="line_abundancies",
    default=None,
    required=False,
    help="Line abundancies data file.",
)

parser.add_argument(
    "-i",
    "--line-intensities",
    metavar="LINE_INTENCITIES",
    dest="line_intensities",
    default=None,
    required=False,
    help="Line intensities data file.",
)

parser.add_argument(
    "-r",
    "--line-ratios",
    metavar="LINE_RATIOS",
    dest="line_ratios",
    default=None,
    required=False,
    help="Line ratios data file.",
)

parser.add_argument(
    "-t",
    "--tene",
    metavar="TENE_DIAGNOSTICS",
    dest="tene",
    default=None,
    required=False,
    help="Temperature/Density data file.",
)

parser.add_argument(
    "-b",
    "--total-abundancies",
    metavar="TOTAL_ABUNDANCIES",
    dest="totalabun",
    default=None,
    required=False,
    help="Total abundancies data file.",
)

parser.add_argument(
    "--bar-plots",
    action="store_true",
    dest="barplot",
    help="Plot using barplots instead of scatter plots.",
)

if __name__ == "__main__":

    # parse cmd
    args = parser.parse_args()

    if args.line_abundancies:
        plotters.plotLineAbundancies(
            args.line_abundancies, "line_abundancies.pdf", args.barplot
        )
    if args.line_intensities:
        plotters.plotLineIntensities(
            args.line_intensities, "line_intensities.pdf", args.barplot
        )
    if args.line_ratios:
        plotters.plotLineRatios(args.line_ratios, "line_ratios.pdf", args.barplot)
    if args.tene:
        plotters.PlotTeNeDiagnostics(
            args.tene, "temperature.pdf", "density.pdf", args.barplot
        )
    if args.totalabun:
        plotters.plotTotalAbundancies(
            args.totalabun, "total_abundancies.pdf", args.barplot
        )
