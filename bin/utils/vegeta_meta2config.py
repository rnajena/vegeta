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

sequences = {}
metaInformation = {}

def extract_structure(string):
  openBr = 0
  start = 0
  stop = 0
  structures = []
  for idx, char in enumerate(string):
    if char in ['[','(','<','{']:
      openBr += 1
      if openBr == 1:
        start = idx
    
    if char in [']',')','>','}']:
      openBr -= 1
      if openBr == 0:
        stop = idx
        structures.append((string[start:stop+2], (start, stop+1)))
  return structures


try:
  with open(stkAln, 'r') as inputStream:

    regionRegex = re.compile(r'\|-*[^\|.]*\|')

    for line in inputStream:
      if line.startswith('//') or not line.strip():
        continue

      if not line.startswith('#'):
        line = line.strip().split()
        #print(line)
        virus = line[0]
        sequence = line[1].replace('.','-')
        sequences[virus] = sequence
        continue


      if not line.startswith('#=GC'):
        continue
    
      splitLine = line.rstrip().split()
      metaRow = splitLine[1]
      
      if metaRow.startswith('SS_cons'):# or metaRow.startswith('Anno'):
        continue
      
      if any([x in line for x in ['[','(','<','{']]):
        metaInformation[f'#=GC {metaRow}'] = extract_structure(splitLine[2])
        continue
      
      splitLine[2] = splitLine[2].replace('.','-')
      annotation = [x.replace('|','') for x in splitLine[2].split('.') if x]
      regions = []
      for x in annotation:
        if not x:
          continue
        if metaRow.startswith('Anno'):
          regions = [y for y in x.split('-') if y]
          regions = sorted(set(regions), key=regions.index)

        else:
          try:
            regions.extend([y for y in x.split('-') if y])
            
          except IndexError:
            print(annotation)
            print(x)
      
      
      matches = regionRegex.finditer(splitLine[2])
      #print(regions)
      regionPositions = []
      for match in matches:
        regionPositions.append((match.start()+1, match.end()))
      metaInformation[f'#=GC {metaRow}'] = list(zip(regions, regionPositions))
        
      
except FileNotFoundError:
  print()
  print("The path to your alignment file is invalid. Please check it!")
  print()
  sys.exit(2)


#print(sequences)


for header, alnSequence in sequences.items():
  with open(f"{header}_{configFile}", 'w') as outputStream:
    outputStream.write(f"#Accession\n{header}\n")
    #print()
    #print(header)
    for metaRow, content in metaInformation.items():
      outputBlock = f"{metaRow}\n"
      for region, span in content:
        #print(metaRow, region, span)
        start, stop = span
        relativeStart = start - alnSequence[:start].count('-')
        relativeStop = stop - alnSequence[:stop].count('-')
        outputBlock += f"{region} {relativeStart} {relativeStop}\n"
      outputBlock += '\n'
      outputStream.write(outputBlock)


# with open(configFile, 'w') as outputStream:
  # for metaRow, content in metaInformation.items():
    # outputBlock = f"{metaRow}\n"
    # for region, span in content:
      # start, stop = span
      # outputBlock += f"{region} {start} {stop}\n"
    # outputBlock += '\n'
    # outputStream.write(outputBlock)
