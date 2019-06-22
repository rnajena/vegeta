#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# god damn, this script looks so ghetto... 
# ToDo: Make it nice and clean at some point.

import sys

try:
  gbFile = sys.argv[1]
  stkAln = sys.argv[2]
except IndexError:
  print("Please give both files, the annotation and the alignment.")
  print("python3 gff2anno.py <GB-File> <STK-Aln>")
  exit(1)


#regionsToAnnotate = ["5'UTR", "3'UTR", "mat_peptide"]
accession = ""
features = {}

alignmentRow = ""

try:
  with open(gbFile, 'r') as inputStream:
    parsing = False
    for line in inputStream:
      if line.startswith("ACCESSION"):
        accession = line.strip().split()[1]
      if line.strip().startswith("FEATURES"):
        
        parsing = True

      if line.strip().startswith("ORIGIN"):
        break

      if not parsing:
        continue

      line = line.strip().split()
      if "UTR" in line[0]:
        features[line[0]] = tuple(map(int, line[1].split('..')))

      if "peptide" in line[0]:
        geneSymbol = inputStream.readline().strip().split("=")[1].strip('"')
        features[geneSymbol] = tuple(map(int, line[1].split('..')))
except FileNotFoundError:
  print("Wasn't able to find the GenBank file. Please check your input")
  exit(2)

print(accession)
print(features)

with open(stkAln, 'r') as inputStream:
  for line in inputStream:
    if line.startswith('#'):
      continue
    if line.startswith(accession):
      alignmentRow = line.strip().split()[1]
      break
      
#print(alignmentRow)

annoString = ["-"]*len(alignmentRow)

nuclCounter = 0
for idx, nt in enumerate(alignmentRow):
  if nt == '-':
    continue
  #if nt != '-':
  nuclCounter += 1
  #else:
  
  for region, positions in features.items():
    start, stop = positions
    if nuclCounter in range(start+20, stop-20, len(region)+50):
      annoString[idx:idx+len(region)] = region
    if nuclCounter in positions:
      annoString[idx] = '|'
      

print(''.join(annoString))