#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import re
import os

try:
  stkAln = sys.argv[1]
  configFile = sys.argv[2]
except IndexError:
  print()
  print("Please give both files, the alignment and the config-file.")
  print("python3 vegeta_config2meta.py <STK-ALN> <CONFIG>")
  print()
  sys.exit(1)

if os.path.isfile(configFile):
  print()
  print("Careful! The config file already exists. Be sure you are not overwriting stuff.")
  print()
  exit(1)

metaInformation = {}

try:
  with open(stkAln, 'r') as inputStream:

    regionRegex = re.compile(r'\|-*[^\|\.]*\|')

    for line in inputStream:
      if not line.startswith('#=GC'):
        continue
      
      splitLine = line.rstrip().split()
      metaRow = splitLine[1]
      
      if metaRow.startswith('SS_cons'):# or metaRow.startswith('Anno'):
        continue
      

      annotation = [x.replace('|','') for x in splitLine[2].split('.') if x]
      regions = []
      for x in annotation:
        if metaRow.startswith('Anno'):
          regions = [y for y in x.split('-') if y]
          regions = sorted(set(regions), key=regions.index)

        else:
          regions.append([y for y in x.split('-') if y][0])
      
      
      matches = regionRegex.finditer(splitLine[2])
      
      regionPositions = []
      for match in matches:
        regionPositions.append((match.start()+1, match.end()))
      metaInformation[f'#=GC {metaRow}'] = zip(regions, regionPositions)
        
      
except FileNotFoundError:
  print()
  print("The path to your alignment file is invalid. Please check it!")
  print()
  sys.exit(2)


with open(configFile, 'w') as outputStream:
  for metaRow, content in metaInformation.items():
    outputBlock = f"{metaRow}\n"
    for region, span in content:
      start, stop = span
      outputBlock += f"{region} {start} {stop}\n"
    outputBlock += '\n'
    outputStream.write(outputBlock)
