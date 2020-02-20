#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
import os
import subprocess
import math
from collections import Counter
import glob

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
    self.alignment = ""
    self.seeds = {}
    self.nonSeeds = {}
    
  def __getstate__(self):
    self_dict = self.__dict__.copy()
    del self_dict['pool']
    return self_dict

  def mafft_scaffold(self):
    """
    """
    with open(f"{self.outdir}/scaffold.aln", 'w') as outputStream:
      cmd = f"mafft --quiet --reorder --clustalout --thread {self.proc} {self.inputFile}"
      subprocess.run(cmd.split(), stdout=outputStream, check=True)
    self.alignment = AlignIO.read(f"{self.outdir}/scaffold.aln", 'clustal')
    self.nonSeeds[0] = self.alignment.get_alignment_length()

  def find_seeds_in_scaffold(self):
    """
    """
    aln = self.alignment
    alnLength = aln.get_alignment_length()
    holyEntropies = {}
    for start in range(0, alnLength-self.seedSize):
      entropies = []
      for colIdx in range(start, start+self.seedSize):
        column = aln[:, colIdx].upper().replace('T','U')
        gaps = (x for x in range(len(aln)))
        gaplessColumn = [x if x != '-' else gaps.__next__() for x in column]
        
        freqs = {nt : 0 for nt in "ACGU"}
        freqs.update(Counter(gaplessColumn))
        entropy = sum(x*math.log2(x) for x in map(lambda x: x/len(aln), freqs.values()) if x != 0) * -1
        entropies.append(entropy)
      holyEntropies.update({start : np.average(entropies)})

    if self.shannon == -1:
      cutoff = np.percentile(list(holyEntropies.values()), 10)
    else:
      cutoff = self.shannon
  
    seedCandidates = {k: k+self.seedSize for k,v in holyEntropies.items() if v < cutoff}
    
    # print(len(seedCandidates))
    # print(len(self.seeds))

    for start, stop in seedCandidates.items():
      self.seeds[start] = stop

    # print(len(self.seeds))
    deleteMe = []
    for idx, start in enumerate(sorted(list(self.seeds))):
      for secondStart in sorted(list(self.seeds))[:idx]:
        if secondStart - self.seedSize < start <= self.seeds[secondStart] + self.seedSize:
          self.seeds[secondStart] = self.seeds[start]
          deleteMe.append(start)
          break

    #print(len(self.seeds), len(deleteMe), len(seedCandidates))
    self.seeds = {start : stop for start, stop in self.seeds.items() if start not in deleteMe}
    
    


    nonSeedStart = 0
    self.nonSeeds = {}
    for idx, start in enumerate(sorted(list(self.seeds))):
      stop = self.seeds[start]
      self.nonSeeds[idx] = (start-1, stop)
      nonSeedStart = stop+1
    self.nonSeeds[idx+1] = (nonSeedStart, self.alignment.get_alignment_length())

  def extract_non_seeds(self):
    """
    """
    for idx, (start, stop) in self.nonSeeds.items():
      alnFragment = self.alignment[:, start:stop]
      with open(f"{self.outdir}/tmpSequences/diverseFragment_{idx}.fasta", 'w') as outputStream:
        for record in alnFragment:
          outputStream.write(f">{record.id}\n{str(record.seq).replace('-','')}\n")
      
      
      
  def refine_fragments(self):
    """
    """
    for idx in self.nonSeeds:
      file = f"{self.outdir}/tmpSequences/diverseFragment_{idx}.fasta"
      start, stop = self.nonSeeds[idx] 
      if stop - start <= 300:
        cmd = f"mlocarna --quiet --stockholm -s 400 --threads {self.proc} {file}"
        subprocess.run(cmd.split(), check=True)

  def merge_fragments(self):
    """
    """
    finalAlignment = {record.id : "" for record in self.alignment}
    print(finalAlignment)


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