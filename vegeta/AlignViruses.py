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
import itertools

import numpy as np
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


from StructureViruses import StructCalculator

class Aligner(object):
  """
  """

  def __init__(self, logger, inputFile, proc, outdir, seedSize, shannon, structureParameter, prefix):
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
    #self.structureParameter = structureParameter 
    self.refinedAlignment = []
    self.prefix = prefix
    
  def __getstate__(self):
    self_dict = self.__dict__.copy()
    del self_dict['pool']
    return self_dict

  def mafft_scaffold(self):
    """
    """
    with open(f"{self.outdir}/{self.prefix}_scaffold.aln", 'w') as outputStream:
      cmd = f"mafft --quiet --reorder --clustalout --thread {self.proc} {self.inputFile}"
      subprocess.run(cmd.split(), stdout=outputStream, check=True)
    self.alignment = AlignIO.read(f"{self.outdir}/{self.prefix}_scaffold.aln", 'clustal')
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
        
        freqs = {nt : 0 for nt in "ACGUN"}
        freqs.update(Counter(gaplessColumn))
        entropy = sum(x*math.log2(x) for x in map(lambda x: x/len(aln), freqs.values()) if x != 0) * -1
        entropies.append(entropy)
      holyEntropies.update({start : np.average(entropies)})

    # if self.shannon == -1:
    #   cutoff = np.percentile(list(holyEntropies.values()), 10)
    # else:
    #   cutoff = self.shannon
    cutoff = np.percentile(list(holyEntropies.values()), self.shannon*100)
  
    seedCandidates = {k: k+self.seedSize for k,v in holyEntropies.items() if v <= cutoff}
    
    for start, stop in seedCandidates.items():
      self.seeds[start] = stop

    
    deleteMe = []
    for idx, start in enumerate(sorted(list(self.seeds))):
      for secondStart in sorted(list(self.seeds))[:idx]:
        if start-1 <= self.seeds[secondStart] < self.seeds[start]:
        #if secondStart < start <= self.seeds[secondStart] + self.seedSize:
          self.seeds[secondStart] = self.seeds[start]          
          deleteMe.append(start)
          break

    #print(len(self.seeds), len(deleteMe), len(seedCandidates))
    self.seeds = {start : stop for start, stop in self.seeds.items() if start not in deleteMe}
    #print(len(self.seeds), len(deleteMe), len(seedCandidates))
    #print(self.seeds)
    #exit(0)


    nonSeedStart = 0
    self.nonSeeds = {}
    for idx, start in enumerate(sorted(list(self.seeds))):
      stop = self.seeds[start]
      self.nonSeeds[idx] = (nonSeedStart, start-1)
      nonSeedStart = stop+1
    self.nonSeeds[idx+1] = (nonSeedStart, self.alignment.get_alignment_length())

  def extract_non_seeds(self):
    """
    """
    for idx, (start, stop) in self.nonSeeds.items():
      alnFragment = self.alignment[:, start:stop+1]
      with open(f"{self.outdir}/tmpSequences/{self.prefix}_diverseFragment_{idx}.fasta", 'w') as outputStream:
        for record in alnFragment:
          outputStream.write(f">{record.id}\n{str(record.seq).replace('-','')}\n")
      
      
      
  def refine_fragments(self, windowSize, stepSize):
    """
    """
    for idx in self.nonSeeds:
      file = f"{self.outdir}/tmpSequences/{self.prefix}_diverseFragment_{idx}.fasta"
      start, stop = self.nonSeeds[idx]
      if all( [ len(str(x.seq).replace('-','')) <= 300 for x in self.alignment[:, start:stop] ] ):
        cmd = f"mlocarna --quiet --stockholm --threads {self.proc} {file}"
        subprocess.run(cmd.split(), check=True)

      else:
        cmd = f"mafft --clustalout --quiet --thread {self.proc} {file}"
        bn = os.path.splitext(os.path.basename(file))[0]
        path = f"{self.outdir}/tmpSequences/{bn}.out/results/result.aln"
        try:
          os.makedirs(f"{self.outdir}/tmpSequences/{bn}.out/results/")
        except FileExistsError:
          self.logger.warning(f"{self.outdir}/tmpSequences/{bn}.out exists! Will overwrite content.")
        with open(path, 'w') as outputStream:
          subprocess.run(cmd.split(), check=True, stdout=outputStream)

  def merge_fragments(self):
    """
    """
    flexible = [(k,*v) for k,v in self.nonSeeds.items()]
    static = [(k,v) for k,v in self.seeds.items()]
    
    if self.nonSeeds[0][0] == 0:
      order = itertools.zip_longest(flexible, static)
    else:
      order = itertools.zip_longest(static, flexible)
    
    finalAlignment = {record.id : "" for record in self.alignment}
    for fragment in order:  
      for element in fragment:
        if not element or element[-1] == -1:
          continue
        if len(element) == 3:
          if element[1] == element[2]:
            continue
          alnFrag = AlignIO.read(f"{self.outdir}/tmpSequences/{self.prefix}_diverseFragment_{element[0]}.out/results/result.aln", 'clustal')
          for record in alnFrag:
            finalAlignment[record.id] += str(record.seq).upper().replace('U','T')
        elif len(element) == 2:
          for record in self.alignment:
            finalAlignment[record.id] += str(record.seq)[element[0]:element[1]+1].upper().replace('U','T')

    #print(finalAlignment)
    for name, sequence in finalAlignment.items():
      self.refinedAlignment.append(SeqRecord(Seq(sequence), id=name))
    self.refinedAlignment = MultipleSeqAlignment(self.refinedAlignment)
    
    AlignIO.write(self.refinedAlignment, f"{self.outdir}/{self.prefix}_refinedAlignment.aln", "clustal")

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
