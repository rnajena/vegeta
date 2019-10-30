#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import math
from collections import Counter

import numpy as np
from Bio import AlignIO

class Aligner(object):
  """
  """

  def __init__(self, logger, inputFile, proc, outdir, seedSize, shannon):
    """
    """

    self.logger = logger
    self.inputFile = inputFile
    self.outdir = outdir
    self.proc = proc
    self.seedSize = seedSize
    self.shannon = shannon

    self.mafft_scaffold()



  def __getstate__(self):
    self_dict = self.__dict__.copy()
    del self_dict['pool']
    return self_dict

  def mafft_scaffold(self):
    """
    """
    with open(f"{self.outdir}/scaffold.aln", 'w') as outputStream:
      cmd = f"mafft --quiet --clustalout --thread {self.proc} {self.inputFile}"
      subprocess.run(cmd.split(), stdout=outputStream, check=True)

  def find_seeds_in_scaffold(self):
    """
    """
    alignment = AlignIO.read(f"{self.outdir}/scaffold.aln", 'clustal')
    alnLength = alignment.get_alignment_length()
    holyEntropies = {}
    for start in range(0, alnLength-self.seedSize):
      entropies = []
      for colIdx in range(start, start+self.seedSize):
        column = alignment[:, colIdx].upper().replace('T','U')
        gaps = (x for x in range(len(alignment)))
        gaplessColumn = [x if x != '-' else gaps.__next__() for x in column]
        
        freqs = {nt : 0 for nt in "ACGU"}
        freqs.update(Counter(gaplessColumn))
        entropy = sum(x*math.log2(x) for x in map(lambda x: x/len(alignment), freqs.values()) if x != 0) * -1
        entropies.append(entropy)
        #print(column, entropy)
      seedEntropy = np.average(entropies)
      holyEntropies.update({start : seedEntropy})
      #if seedEntropy < self.shannon:
      #  print(start, start+self.seedSize)
      #  print(alignment[:, start:start+self.seedSize])
      #  print(seedEntropy)
      # print()
    if self.shannon == -1:
      self.shannon = np.percentile(list(holyEntropies.values()), 5)
    seeds = {k : v for k,v in holyEntropies.items() if v < self.shannon}
    for start, entropy in seeds.items():
      print(start, start+self.seedSize)
      print(alignment[:, start:start+self.seedSize])
      print(entropy)
      print()
      
    

  def read_sequences(self):
    """
    """

    fastaContent = {}
    
    with open(self.inputFile, 'r') as inputStream:
      header = ''
      seq = ''

      for line in inputStream:
        if line.startswith(">"):
          if header:
              fastaContent[header] = seq
          header = line.rstrip("\n").replace(':','_').replace(' ','_').lstrip(">")
          seq = ''
        else:
          seq += line.rstrip("\n").upper().replace('U','T')

      fastaContent[header] = seq
    return fastaContent