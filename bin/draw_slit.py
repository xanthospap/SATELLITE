#!/usr/bin/env python

from satellite import version
import satellite.fitsutils as fs

import os
import sys

slit = {'PA': 79, 'w': 8, 'h': 25, 'x': 37, 'y': 148}

fs.testGetVerticalSlit(sys.argv[1], slit)
