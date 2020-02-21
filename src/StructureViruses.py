# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
import os
import subprocess

from Bio import AlignIO
import RNA

class StructCalculator(object):
  """
  """

  def __init__(self, path, logger, outdir, windowSize, stepSize, proc):
    """
    """

    self.logger = logger
    self.outdir = outdir
    self.windowSize = windowSize 
    self.stepSize = stepSize
    self.proc = proc
    self.path = path
    self.alignment = self.__read_alignment(path)
    self.finalStructure = '.' * self.alignment.get_alignment_length()
    self.overlappingStructures = {}

  def __read_alignment(self, path):
    """
    """

    return(AlignIO.read(path, 'fasta'))

  def apply_lalifold(self):
    """
    """

    cmd = f"RNALalifold -L {self.windowSize} {self.path}"
    lalifoldResult = subprocess.getoutput(cmd)
    structureWindows = {}
    windows = lalifoldResult.split("\n")[1:-1]
    for structure in windows:
      structure = structure.split()
      start = int(structure[-3])
      structureWindows[start] = structure[0]

    nonOverlap = {}
    lastStart = -1
    lastStop = -1

    for start in sorted(list(structureWindows)):
      if lastStart < start <= lastStop:
        nonOverlap[lastStart] = start + len(structureWindows[start]) -1
      else:
        nonOverlap[start] = start + len(structureWindows[start]) -1
        lastStart = start
      lastStop = start + len(structureWindows[start]) -1

    for start, stop in nonOverlap.items():
      if stop == len(structureWindows[start])+start-1:
        self.finalStructure = self.finalStructure[:start] + structureWindows[start] + self.finalStructure[stop+1:]
      else:
        self.overlappingStructures[start] = stop
    
  def resolve_conflicts(self):
    """
    """
    for start, stop in self.overlappingStructures.items():
      alignmentFragment = [str(x.seq) for x in self.alignment[:, start:stop]]
      alifoldObject = RNA.fold_compound(alignmentFragment)
      structure, _ = alifoldObject.mfe()
      self.finalStructure = self.finalStructure[:start] + structure + self.finalStructure[stop+1:]