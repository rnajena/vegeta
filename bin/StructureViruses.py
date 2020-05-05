#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
"""

import sys
import os
import subprocess
import shutil
import glob
import math
import itertools

import numpy as np
from Bio import AlignIO

import RNA

class StructCalculator(object):
  """
  """

  def __init__(self, path, logger, outdir, windowSize, stepSize, proc, allowLP, tbpp, prefix):
    """
    """
    self.logger = logger
    self.prefix = prefix
    self.outdir = outdir
    self.windowSize = windowSize
    self.allowLP = allowLP
    self.stepSize = stepSize
    self.proc = proc
    self.path = path
    self.tbpp = tbpp
    self.alignment = self.__read_alignment(path)
    self.alnLength = self.alignment.get_alignment_length()
    self.finalStructure = '.' * self.alnLength
    self.windows = self.__create_sliding_window()
    self.overlappingStructures = {}
    self.nonOverlap = {}
    self.bpps = {}
    self.used_basepairs = []
    

  def __read_alignment(self, path):
    """
    """

    return(AlignIO.read(path, 'clustal'))

  def __create_sliding_window(self):
    """
    """
    windows = [(start, start+self.windowSize) for start in range(0, self.alnLength-self.windowSize, self.stepSize)] + [(self.alnLength-self.windowSize,self.alnLength)]
    for idx,(start,stop) in enumerate(windows):
      with open(f"{self.outdir}/tmpSequences/{self.prefix}_window_{idx}.aln" ,'w') as outputStream:
        fragment = self.alignment[:, start:stop]
        for record in fragment:
          outputStream.write(f">{record.id}\n{record.seq}\n")

    return(windows)
  

  def apply_alifold(self):
    """
    """
    TRASH = open(os.devnull, 'w')
    for idx, _ in enumerate(self.windows):
      file = f"{self.outdir}/tmpSequences/{self.prefix}_window_{idx}.aln"
      cmd = f"RNAalifold --noLP -p --cfactor 0.6 --nfactor 0.5 -r --id-prefix={self.prefix}_window_{idx} {file}"
      subprocess.run(cmd.split(), check=True, stderr=TRASH, stdout=TRASH)
      os.remove(f"{self.prefix}_window_{idx}_0001_ali.out")
      os.remove(f"{self.prefix}_window_{idx}_0001_ss.ps")
      shutil.move(f"{self.prefix}_window_{idx}_0001_dp.ps", f"{self.outdir}/tmpSequences/")
    TRASH.close()


  def calculate_avg_bpp(self):
    """
    """
    
    for file in glob.glob(f"{self.outdir}/tmpSequences/{self.prefix}_window_*.ps"):
      
      idx = int(os.path.basename(file).split('_')[-3])
      currentWindow = self.windows[idx]
      with open(file, 'r') as inputStream:
        for line in inputStream:
          line = line.strip()
          if not line.endswith("ubox") or line.startswith("%"):
            continue
          line = line.split()
          start = int(line[3]) + currentWindow[0]
          stop = int(line[4]) + currentWindow[0]
          bpp = math.pow(float(line[5]), 2)

          if bpp >= self.tbpp:
            self.__update_bpps(start, stop, bpp)
            self.__update_bpps(stop, start, bpp)

    for start, values in self.bpps.items():
      for stop, probabilites in values.items():
        self.bpps[start][stop] = np.average(probabilites)

  def __update_bpps(self, x, y, bpp):
    """
    """
    if x not in self.bpps:
      self.bpps[x] = {y : [bpp]}
    elif y not in self.bpps[x]:
      self.bpps[x].update({y : [bpp]})
    else:
      self.bpps[x][y].append(bpp)   

  def generate_ilp(self):
    """
    """
    ilp = ILP(self.bpps, self.outdir, self.prefix)
    self.used_basepairs = ilp.used_basepairs
    print(sorted(self.used_basepairs))
    exit(0)

  def finalize_structure(self):
    """
    """
    allStarts = [x[0] for x in self.used_basepairs]
    allStops = [x[1] for x in self.used_basepairs]
    structure = list(self.finalStructure)

    for start, stop in self.used_basepairs:
      if not self.allowLP:
        startRange = [start-1, start+1]
        stopRange = [stop-1, stop+1]
        if not (any([x in allStarts for x in startRange]) or any([x in allStops for x in stopRange])):
          continue
      if stop-start < 4: 
        continue
      structure[start] = '('
      structure[stop] = ')'
      self.finalStructure = self.finalStructure[:start] + '(' + self.finalStructure[start+1:stop] + ')' + self.finalStructure[stop+1:]
    
    localStructures = []
    currentStructure = []

    for idx, char in enumerate(structure):
      if char == '(':
        #print(currentStructure)
        if not currentStructure:
          start = idx
        currentStructure.append(idx)
      if char == ')':
        #print(currentStructure)
        currentStructure.pop()
        if not currentStructure:
          localStructures.append((start, idx))

    for start, stop in localStructures:
      print(start, stop)
      print(self.finalStructure[start:stop+1])
      print(self.finalStructure[start:stop+1].count('('), self.finalStructure[start:stop+1].count(')'))
    
    print()
    print(self.finalStructure[stop+1:])
    print(self.finalStructure[stop+1:].count('('), self.finalStructure[stop+1:].count(')'))
    exit(0)
    

class ILP(object):
  """
  """

  def __init__(self, bpp_dict, outdir, prefix):
    """
    """
    self.bpp_dict = bpp_dict
    self.outdir = outdir
    self.used_basepairs = []
    self.prefix = prefix
    isSolved = self.generate_ilp()
    if not isSolved:
      self.solve_ilp()
      self.extract_used_basepairs()

  def generate_ilp(self):
    """
    """

    trivialCases = {}
    bpp_iterator = sorted(list(self.bpp_dict))
    #print(bpp_iterator)
    for idx, start in enumerate(bpp_iterator):
      values = self.bpp_dict[start]
      if len(values) == 1:
        interactionPartner = list(values.keys())[0]
        if len(self.bpp_dict[interactionPartner]) == 1:
          posBetween = [x for x in range(start+1, interactionPartner) if x in self.bpp_dict]
          if not any([ y < start or y > interactionPartner for x in posBetween for y in self.bpp_dict[x].keys()]):
            trivialCases[start] = interactionPartner
            trivialCases[interactionPartner] = start
            if start < interactionPartner:
              self.used_basepairs.append((start, interactionPartner))
    
    filteredBPPs = {pos : values for pos, values in self.bpp_dict.items() if pos not in trivialCases}
    if not filteredBPPs:
      return(True)
    bpp_iterator = sorted(list(filteredBPPs))
    

    connectedComponents = []
    newComponent = set()

    while bpp_iterator:
      left = bpp_iterator[0]  
      maxRight = max(self.bpp_dict[left])
      #inBetween = [left] + [x for x in range(left, maxRight) if x in self.bpp_dict]
      inBetween = [left] + [x for x in range(left, maxRight) if x in filteredBPPs]
      currentLength = len(inBetween)
      #print(inBetween)
      while 1:
        maxBPPs = [max(list(self.bpp_dict[x].keys())) for x in inBetween]
        newMaxRight = max(maxBPPs)
        #inBetween = [left] + [x for x in range(left, newMaxRight) if x in self.bpp_dict]
        inBetween = [left] + [x for x in range(left, newMaxRight) if x in filteredBPPs]

        if currentLength == len(inBetween):
          for x in inBetween:
            newComponent.add(x)

          bpp_iterator = [x for x in bpp_iterator if x not in newComponent]
          connectedComponents.append(sorted(newComponent))
          newComponent = set()  
          break
        else:
          currentLength = len(inBetween)
    
    #print(connectedComponents)
    for idx, component in enumerate(connectedComponents):
      if len(component) == 1:
        continue
      bpp_iterator = list(component)
      with open(f"{self.outdir}/tmpSequences/{self.prefix}_structure_{idx}.ilp", 'w') as outputStream:
        outputStream.write("Maximize\n")
        outputStream.write("obj: ")
        edges = []
        #bpp_iterator = sorted(list(self.bpp_dict))
        for start in bpp_iterator:
          values = self.bpp_dict[start]
          for stop, probability in values.items():
            if stop <= start or stop-start < 4:
              continue
            outputStream.write(f"+ {probability} e_{start}_{stop} ")
            edges.append(f"e_{start}_{stop}")
      
        outputStream.write("\nSubject To\n")
        variableCounter = 1
        conflictCounter = 0
        for idx, start in enumerate(bpp_iterator):
          values = self.bpp_dict[start]
          if len(values) != 1:
            constraint = ' + '.join([f'e_{start}_{x}' for x in map(str, values)])
            outputStream.write(f"c{variableCounter}: {constraint} <= 1\n")
            variableCounter += 1
        
          highest_partner = sorted(list(values), reverse=True)[0]
        
          for potentialConflictStart in bpp_iterator[idx:]:
            if potentialConflictStart >= highest_partner:
              break
            potentialConflictValues = self.bpp_dict[potentialConflictStart]
            for stop in values:
              if stop <= start:
                continue
              for potentialConflictStop in potentialConflictValues:
                if start < potentialConflictStart < stop < potentialConflictStop:
                  conflictCounter += 1
                  outputStream.write(f"c{variableCounter}: e_{start}_{stop} + e_{potentialConflictStart}_{potentialConflictStop} <= 1\n")
                  variableCounter += 1

        if variableCounter == 1:
          outputStream.write(f"c{variableCounter}: {edges[0]} <= 1\n")

        outputStream.write("\nBinary\n")
        for edge in edges:
          outputStream.write(f"{edge}\n")
        outputStream.write("End\n")

    return(False)
        
      
  def solve_ilp(self):
    """
    """
    TRASH = open(os.devnull, 'w')
    for file in glob.glob(f"{self.outdir}/tmpSequences/{self.prefix}_*.ilp"):
      bn = os.path.basename(file)
      cmd = f"glpsol --lp {file} --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim 14400 -o {self.outdir}/tmpSequences/{bn}.sol"
      subprocess.run(cmd.split(), check=True, stderr=TRASH, stdout=TRASH)
    TRASH.close()
    

  def extract_used_basepairs(self):
    """
    """
    for file in glob.glob(f"{self.outdir}/tmpSequences/{self.prefix}_*.sol"):
      with open(file, 'r') as inputStream:
        for line in inputStream:
          line = line.strip().split()
          if len(line) < 4: 
            continue
          if line[1].startswith("e_") and int(line[3]) == 1:
            start, stop = line[1].split('_')[1:]
            self.used_basepairs.append((int(start),int(stop)))
          