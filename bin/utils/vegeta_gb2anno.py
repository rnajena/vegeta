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
  print("python3 vegeta_gb2anno.py <GB-File> <STK-Aln>")
  exit(1)


accession = ""
features = {}

ranges = []
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

      try:
        if "peptide" in line[0]:
          geneSymbol = ''
          while not "product" in geneSymbol:
            geneSymbol = inputStream.readline().strip().split("=")[1].strip('"').replace(' ','_')
          startAnno = int(line[1].split('..')[0])
          stopAnno = int(line[1].split('..')[1])
          if not any([startAnno in x for x in ranges]):
            features[geneSymbol] = tuple(map(int, line[1].split('..')))
            ranges.append(range(startAnno,stopAnno))
            

      except ValueError:
        continue
        #print(print(line))
        #exit(1)
except FileNotFoundError:
  print("Wasn't able to find the GenBank file. Please check your input")
  exit(2)

print(accession)
print(features)
wholeAlignment = ''

alnLength = 0
alnStart = 0

with open(stkAln, 'r') as inputStream:
  wholeAlignment = ''.join(inputStream.readlines()).rstrip('/\n')
  inputStream.seek(0)
  for line in inputStream:
    if line.startswith('#'):
      continue
    if line.startswith(accession):
      alignmentRow = line.strip().split()[1]
      alnLength = len(line.rstrip().split()[1])
      alnStart = len(line.rstrip()) - alnLength
      break
      
annoString = ["-"]*len(alignmentRow)

nuclCounter = 0
for idx, nt in enumerate(alignmentRow):
  if nt == '-':
    continue
  nuclCounter += 1
  
  for region, positions in features.items():
    start, stop = positions

    if nuclCounter in positions:
      annoString[idx] = '|'

    if nuclCounter in range(start+20, stop-20, len(region)+50):
      annoString[idx:idx+len(region)] = region

annoString[-1] = '|'
      
with open(stkAln, 'w') as outputStream:
  annoString = f"#=GC Annotation{' '*(alnStart-len('#=GC Annotation'))}{''.join(annoString)}"
  outputStream.write(f"{wholeAlignment}\n{annoString}\n//")