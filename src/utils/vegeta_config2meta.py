#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys

try:
  stkAln = sys.argv[1]
  configFile = sys.argv[2]
except IndexError:
  print()
  print("Please give both files, the alignment and the config-file.")
  print("python3 vegeta_config2meta.py <STK-ALN> <CONFIG>")
  print()
  sys.exit(1)


accession = ''
metaRow = {}

try:
  with open(configFile, 'r') as inputStream:
    rowName = ''
    rowContent = {}

    for line in inputStream:

      if line.startswith('##') or not line.strip():
        continue

      if line.startswith('#Accession'):
        accession = inputStream.readline().rstrip()
        continue

      if line.startswith('#=GC'):
        if rowName and not rowContent:
          metaRow[rowName] = {}
        if rowContent:
          metaRow[rowName] = rowContent
          rowContent = {}
        rowName = line.rstrip()
        continue
      splitLine = line.rstrip().split()
      symbol = splitLine[0]
      positions = tuple(map(int,splitLine[1:]))
      rowContent[symbol] = positions
    metaRow[rowName] = rowContent

except FileNotFoundError:
  print()
  print("The path to your config file is invalid. Please check it!")
  print()
  sys.exit(2)


alnLength = 0
alnStart = 0
aln = ''
#print(accession)

try:
  with open(stkAln, 'r') as inputStream:
    wholeAln = inputStream.readlines()
    wholeAln = ''.join([x for x in wholeAln if not any([x.startswith(y) for y in metaRow])]).rstrip('/\n')
    inputStream.seek(0)
    for line in inputStream:
      if line.startswith('#') or not line.strip():
        continue
      
      if line.startswith(accession):
        aln = line.rstrip().split()[1]
        alnLength = len(aln)
        alnStart = len(line.rstrip()) - alnLength
        break  
except FileNotFoundError:
  print()
  print("The path to your alignment file is invalid. Please check it!")
  print()
  sys.exit(2)


printRow = ''

for row, d_annotation in metaRow.items():
  metaAnno = ['.']*alnLength
  nuclCounter = 0
  for idx, nt in enumerate(aln):
    if nt == '-':
      continue

    nuclCounter += 1

    for symbol, positions in d_annotation.items():
      start, stop = positions

      #print(nuclCounter)
      if nuclCounter in positions:
        metaAnno[idx] = '|'
      #if start < nuclCounter < stop:
      #  metaAnno[idx+1:idx+stop-start] = '-'*(stop-start)

      j = len(symbol)
      
      if nuclCounter in range(start+2, stop-2, len(symbol)+20):
      #for i in range(idx+2, idx+start-stop-2, j):
        metaAnno[idx:idx+len(symbol)] = symbol

  printRow += row + ' '*(alnStart-len(row)) + ''.join(metaAnno) + "\n"

with open(stkAln, 'w') as outputStream:
  outputStream.write(f"{wholeAln}\n{printRow}\n//\n")


