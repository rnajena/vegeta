#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import glob

configFile = sys.argv[1]
rfamResults = sys.argv[2]

foundFamily = {}

for file in glob.glob(f"{rfamResults}/*out"):
  with open(file, 'r') as inputStream:
    for line in inputStream:
      if line.startswith('#'):
        continue
      splitLine = line.rstrip().split()
      if splitLine[-2] == '!':
        foundFamily[splitLine[3]] = (splitLine[7], splitLine[8])


with open(configFile, 'a') as outputStream:
  outputStream.write("#=GC RFam_families\n")
  for family, position in foundFamily.items():
    start, stop = position
    outputStream.write(f"{family} {start} {stop}\n")