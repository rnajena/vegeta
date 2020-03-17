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

import numpy as np
from Bio import AlignIO

import RNA

class StructCalculator(object):
  """
  """

  def __init__(self, path, logger, outdir, windowSize, stepSize, proc, allowLP, tbpp):
    """
    """
    self.logger = logger
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
      with open(f"{self.outdir}/tmpSequences/window_{idx}.aln" ,'w') as outputStream:
        fragment = self.alignment[:, start:stop]
        for record in fragment:
          outputStream.write(f">{record.id}\n{record.seq}\n")

    return(windows)
  

  def apply_alifold(self):
    """
    """
    for idx, _ in enumerate(self.windows):
      file = f"{self.outdir}/tmpSequences/window_{idx}.aln"
      cmd = f"RNAalifold --noLP -p --cfactor 0.6 --nfactor 0.5 -r --id-prefix=window_{idx} {file}"
      subprocess.run(cmd.split(), check=True)
    for file in os.listdir():
      if file.endswith("dp.ps"):
        shutil.move(file, f"{self.outdir}/tmpSequences/{file}")
      elif file.endswith("out") or file.endswith("ss.ps"):
        os.remove(file)


  def calculate_avg_bpp(self):
    """
    """
    
    for file in glob.glob(f"{self.outdir}/tmpSequences/window_*.ps"):
      
      idx = int(os.path.basename(file).split('_')[1])
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
    ilp = ILP(self.bpps, self.outdir)
    self.used_basepairs = ilp.used_basepairs

  def finalize_structure(self):
    """
    """
    allStarts = [x[0] for x in self.used_basepairs]
    allStops = [x[1] for x in self.used_basepairs]

    for start, stop in self.used_basepairs:
      if not self.allowLP:
        startRange = [start-1, start+1]
        stopRange = [stop-1, stop+1]
        if not (any([x in allStarts for x in startRange]) or any([x in allStops for x in stopRange])):
          if stop-start >= 4:
            continue

      self.finalStructure = self.finalStructure[:start] + '(' + self.finalStructure[start+1:stop] + ')' + self.finalStructure[stop+1:]
    
    #print(self.used_basepairs)
    print(self.finalStructure)

class ILP(object):
  """
  """

  def __init__(self, bpp_dict, outdir):
    """
    """
    self.bpp_dict = bpp_dict
    self.outdir = outdir
    self.used_basepairs = []
    self.generate_ilp()
    self.solve_ilp()
    self.extract_used_basepairs()

  def generate_ilp(self):
    """
    """

    trivialCases = {}
    bpp_iterator = sorted(list(self.bpp_dict))
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
    bpp_iterator = sorted(list(filteredBPPs))

    with open(f"{self.outdir}/tmpSequences/structure.ilp", 'w') as outputStream:
      outputStream.write("Maximize\n")
      outputStream.write("obj: ")
      edges = []
      #bpp_iterator = sorted(list(self.bpp_dict))
      for start in bpp_iterator:
        values = self.bpp_dict[start]
        for stop, probability in values.items():
          if stop <= start:
            continue
          outputStream.write(f"+ {probability} e_{start}_{stop} ")
          edges.append(f"e_{start}_{stop}")
      
      outputStream.write("\nSubject To\n")
      variableCounter = 1
      #bpp_iterator = sorted(list(self.bpp_dict))
      conflictCounter = 0
      for idx, start in enumerate(bpp_iterator):
        values = self.bpp_dict[start]
        #if len(values) == 1:
          #interactionPartner = list(values.keys())[0]
          #if len(self.bpp_dict[interactionPartner]) == 1:
            #posBetween = [x for x in range(start+1, interactionPartner) if x in self.bpp_dict]
            #for x in posBetween:
              #if x in self.bpp_dict:
              #for y in self.bpp_dict[x].keys():
                #print(start, x, y, interactionPartner)
            #if not any([ y < start or y > interactionPartner for x in posBetween for y in self.bpp_dict[x].keys()]):
              # outputStream.write(f"c{variableCounter}: e_{start}_{interactionPartner} = 1\n")
              #variableCounter += 1
              #continue
              #print(start, interactionPartner)
          #exit(0)
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
              #if not (start < potentialConflictStart < potentialConflictStop < stop or start < stop < potentialConflictStart < potentialConflictStop):
                #print(f"c{variableCounter}: e_{start}_{stop} + e_{potentialConflictStart}_{potentialConflictStop} <= 1")
                conflictCounter += 1
                outputStream.write(f"c{variableCounter}: e_{start}_{stop} + e_{potentialConflictStart}_{potentialConflictStop} <= 1\n")
                variableCounter += 1
      
      outputStream.write("\nBinary\n")
      for edge in edges:
        outputStream.write(f"{edge}\n")
      outputStream.write("End\n")
        
      
  def solve_ilp(self):
    """
    """
    cmd = f"glpsol --lp {self.outdir}/tmpSequences/structure.ilp --mipgap 0.01 --pcost --cuts --memlim 16834 --tmlim 14400 -o  {self.outdir}/tmpSequences/structure.ilp.sol"
    subprocess.run(cmd.split(), check=True)

  def extract_used_basepairs(self):
    """
    """
    
    with open(f"{self.outdir}/tmpSequences/structure.ilp.sol", 'r') as inputStream:
      for line in inputStream:
        line = line.strip().split()
        if len(line) < 4: 
          continue
        if line[1].startswith("e_") and int(line[3]) == 1:
          start, stop = line[1].split('_')[1:]
          self.used_basepairs.append((int(start),int(stop)))
          


  # def finalize_structure(self):
  #   """
  #   """
  #   print([len(x) for _,x in self.bpps.items()])


  # def apply_lalifold(self):
  #   """
  #   """

  #   cmd = f"RNALalifold -L {self.windowSize} {self.path}"
  #   lalifoldResult = subprocess.getoutput(cmd)
  #   structureWindows = {}
  #   windows = lalifoldResult.split("\n")[1:-1]
  #   for structure in windows:
  #     structure = structure.split()
  #     start = int(structure[-3])
  #     structureWindows[start] = structure[0]

  #   #nonOverlap = {}
  #   lastStart = -1
  #   lastStop = -1

  #   for start in sorted(list(structureWindows)):
  #     if lastStart < start <= lastStop:
  #       self.nonOverlap[lastStart] = start + len(structureWindows[start]) -1
  #     else:
  #       self.nonOverlap[start] = start + len(structureWindows[start]) -1
  #       lastStart = start
  #     lastStop = start + len(structureWindows[start]) -1

  #   for start, stop in self.nonOverlap.items():
  #     if stop == len(structureWindows[start])+start-1:
  #       self.finalStructure = self.finalStructure[:start] + structureWindows[start] + self.finalStructure[stop+1:]
  #     else:
  #       self.overlappingStructures[start] = stop
    
  #   self.resolve_conflicts()

  # def resolve_conflicts(self):
  #   """
  #   """
  #   for start, stop in self.overlappingStructures.items():
  #     alignmentFragment = [str(x.seq) for x in self.alignment[:, start:stop]]
  #     alifoldObject = RNA.fold_compound(alignmentFragment)
  #     structure, _ = alifoldObject.mfe()
  #     self.finalStructure = self.finalStructure[:start] + structure + self.finalStructure[stop+1:]