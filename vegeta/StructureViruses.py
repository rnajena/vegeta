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
import random
import itertools

import numpy as np
from scipy.stats import norm
from scipy.stats.mstats import zscore

from Bio import AlignIO
import RNA

class StructCalculator(object):
  """
  """

  def __init__(self, path, logger, outdir, windowSize, stepSize, proc, allowLP, tbpp, prefix, shuffle, pvalue):
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
    self.shuffle = shuffle
    self.pvalue = pvalue
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
          if idx == 0 or idx == len(windows)-1:
            groupsInSeq = itertools.groupby(str(record.seq))
            result = max([0] + [sum(1 for _ in group) for gap, group in groupsInSeq if gap == '-'])
            if result >= 15:
              continue
          outputStream.write(f">{record.id}\n{record.seq}\n")

    return(windows)
  

  def apply_alifold(self):
    """
    """
    TRASH = open(os.devnull, 'w')
    
    for idx, _ in enumerate(self.windows):
      OUTPUT = open(f"{self.outdir}/tmpSequences/{self.prefix}_window_{idx}_alifold.out", 'w')
      file = f"{self.outdir}/tmpSequences/{self.prefix}_window_{idx}.aln"
      cmd = f"RNAalifold --noLP -p2 --MEA --cfactor 0.6 --nfactor 0.5 -r --id-prefix={self.prefix}_window_{idx} {file}"
      subprocess.run(cmd.split(), check=True, stdout=OUTPUT, stderr=TRASH)
      os.remove(f"{self.prefix}_window_{idx}_0001_ss.ps")
      os.remove(f"{self.prefix}_window_{idx}_0001_ali.out")
      shutil.move(f"{self.prefix}_window_{idx}_0001_dp.ps", f"{self.outdir}/tmpSequences/")
      OUTPUT.close()
    TRASH.close()

  def extract_structures(self):
    """
    """
    
    for file in glob.glob(f"{self.outdir}/tmpSequences/{self.prefix}*_alifold.out"):
      idx = int(os.path.basename(file).split('_')[-2])
      currentWindow = self.windows[idx]
      offset = currentWindow[0]
      with open(file, 'r') as inputStream:
        currentLine = 0
        for line in inputStream:
          if not line or not line.startswith('.') and not line.startswith('('):
            continue
          currentLine += 1
          if currentLine == 1:# or currentLine == 4: # MFE and MEA structure, respecitvely.
            bracketStructure = line.split()[0]
            start = 0
            stack = 0
            for idx, char in enumerate(bracketStructure):
              if char == '(':
                if not stack:
                  start = idx
                stack += 1
              elif char == ')':
                stack -= 1
                if not stack:
                  closedStructure = (start  + offset, idx + offset)
                  self.overlappingStructures[closedStructure] = bracketStructure[start:idx+1]
    print(self.get_all_nonoverlapping_structures())


  def get_all_nonoverlapping_structures(self, searchList=None, result=[]):
    """
    """
    if not searchList:
      searchList = sorted(self.overlappingStructures.keys())
    #searchList = sorted(self.overlappingStructures.keys())
    # stopping condition #1 ... last element was overlapped and is a singleton
    if len(searchList) == 1:
        result.append(searchList[0])
        return result
    temp = [searchList[0]]  # to hold intermediate result, starting with first input element
    for idx, element in enumerate(searchList[1:], 1):
        if element[0] >= temp[-1][1]+1:
            temp.append(element)
        #else:  # found a new starting point
         #   self.get_all_nonoverlapping_structures(searchList[idx:], result)
    result.append(tuple(temp))
    return(result)


  def calculate_avg_bpp(self):
    """
    """
    
    for file in glob.glob(f"{self.outdir}/tmpSequences/{self.prefix}_window_*.ps"):
      
      idx = int(os.path.basename(file).split('_')[-3])
      #if idx != 0:
      #  continue
      currentWindow = self.windows[idx]
      with open(file, 'r') as inputStream:
        for line in inputStream:
          line = line.strip()
          if not line.endswith("ubox") or line.startswith("%"):
            continue
          line = line.split()
          start = int(line[3]) -1 + currentWindow[0] # -1, since the dot.ps file counts nucleotides starting at 1
          stop = int(line[4]) -1 + currentWindow[0] # -1, since the dot.ps file counts nucleotides starting at 1
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
    with open(f"{self.outdir}/tmpSequences/parsed_bpp.txt", 'w') as outStream:
      outStream.write("\n".join(map(str, sorted(list(self.bpps.items())))))
    ilp = ILP(self.bpps, self.outdir, self.prefix)
    self.used_basepairs = ilp.used_basepairs
    #print(sorted(self.used_basepairs))
    #exit(0)

  def __rna_alifold(self, query):
    """
    Wrapper for the ViennaRNA API method fold_compound().
    The ribosum_scoring and noLP parameters are set to 1 first, then the sequence (or alignment)
    is folded. Furthermore, the covariance of the mfe structure is calculated.
    rna_alifold() returns the mfe, the covariance score and the structure itself.

    Keyword arguments:
    ary -- the query (either sequence or alignment) that has to be folded.
    """
    if all([len(x) == 0 for x in query]):
        #logger.warn("Found an empty alignment window.")
        return[10,0,'']

    RNA.cvar.ribo = 1
    RNA.cvar.noLP = 1
    fc = RNA.fold_compound(query)

    structure, mfe = fc.mfe()
    covar = fc.eval_covar_structure(structure)
    return [mfe, covar, structure]

  def __soft_shuffle(self, aln, shuffle):
    """
    Soft shuffles a given alignment aln and calculates the z-score based p-value
    for the energy of the consensus sequence. The soft shuffle applies at least 0.1*len(aln) and
    at most 0.4*len(aln) changes to the alignment.

    Keyword arguments:
    aln -- query alignment that has to be shuffled
    mfe -- the mfe of the corresponding consensus secondary structure
    covar -- covariance score of the corresponding consensus secondary structure
    shuffle -- determines how many different sequences are generated
    """
    aln = [str(x.seq) for x in aln]
    mfe, covar, structure = self.__rna_alifold(aln)
    aln = list(map(list, aln))

    wi = len(aln[0])
    min_shuffle = int(wi * 0.1)
    
    mfes = [mfe - covar]

    for i in range(0, shuffle):
        # print("z: {}".format(i))
        k = random.randint(min_shuffle, min_shuffle + int((wi - min_shuffle) * 0.45))
        ary = np.array(aln).T

        for j in range(0, k):
            p = random.sample(range(wi), 2)
            tmp = np.copy(ary[p[0]])
            ary[p[0]] = ary[p[1]]
            ary[p[1]] = tmp

        new_ary = list(map("".join, ary.T))
        mfe, covar, structure = self.__rna_alifold(new_ary)
        mfes.append(mfe - covar)

    a = np.array(mfes)
    z = zscore(a)[0]
    p_values = norm.sf(abs(z)) * 2

    return z, p_values

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
      #self.finalStructure = self.finalStructure[:start] + '(' + self.finalStructure[start+1:stop] + ')' + self.finalStructure[stop+1:]
    
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
    np.seterr(all='ignore')
    for start, stop in localStructures:
      fragment = self.alignment[:, start:stop+1]
      #print(self.finalStructure[start:stop+1])


      ############################
      # zscore, pvalue = self.__soft_shuffle(fragment, self.shuffle)
      # if pvalue > self.pvalue:
      #   for i in range(start,stop+1):
      #     structure[i] = '.'
      ############################

      #exit(0)
      #print(start, stop)
      #print(self.finalStructure[start:stop+1].count('('), self.finalStructure[start:stop+1].count(')'))
    
    #print()
    #print(self.finalStructure[stop+1:])
    #print(self.finalStructure[stop+1:].count('('), self.finalStructure[stop+1:].count(')'))
    #print()
    self.finalStructure = ''.join(structure)
    #print(self.finalStructure)
    #exit(0)

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
        interactionPartner = {}
        for x in inBetween:
          interactionPartner.update({y : bpp for y,bpp in self.bpp_dict[x].items() if y>maxRight and bpp>=0.7})
        #print(inBetween)
        #print(interactionPartner)
        if interactionPartner:
          maxBPPs = [max(list(interactionPartner))]
          newMaxRight = max(maxBPPs)
        else:
          newMaxRight = maxRight
        #maxBPPs = [max(list(self.bpp_dict[x].keys())) for x in inBetween if ]
        #newMaxRight = max(maxBPPs)
        
        
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
    print(len(filteredBPPs))
    print(len(connectedComponents))
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
              #print(f"e_{start}_{stop}")
              continue
            outputStream.write(f"+ {probability} e_{start}_{stop} ")
            edges.append(f"e_{start}_{stop}")
      
        outputStream.write("\nSubject To\n")
        variableCounter = 1
        conflictCounter = 0
        positions = set()
        for edge in edges:
          start = edge.split('_')[1]
          stop = edge.split('_')[2]
          positions.add(int(start))
          positions.add(int(stop))
          outputStream.write(f"c{variableCounter}: e_{start}_{stop} - e_{stop}_{start} = 0\n")
          variableCounter += 1
          # outputStream.write(f"c{variableCounter}: p_{start} + p_{int(start)-1} + p_{int(start)+1} > 1\n")
          # variableCounter += 1
          # outputStream.write(f"c{variableCounter}: p_{stop} + p_{int(stop)-1} + p_{int(stop)+1}  > 1\n")
          # variableCounter += 1
          

        for idx, start in enumerate(bpp_iterator):
          values = self.bpp_dict[start]
          if len(values) != 1:
            constraint = ' + '.join([f'e_{start}_{x}' for x in map(str, values)]) #if f'e_{start}_{x}' in edges])
            if constraint:
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

        for pos in positions:
          outputStream.write(f"c{variableCounter}: p_{pos} + p_{pos-1} + p_{pos+1} > 1\n")
          variableCounter += 1

        if variableCounter == 1:
          outputStream.write(f"c{variableCounter}: {edges[0]} <= 1\n")

        outputStream.write("\nBinary\n")
        for edge in edges:
          outputStream.write(f"{edge}\n")
        for pos in positions:
           outputStream.write(f"p_{pos}\n")
           outputStream.write(f"p_{pos-1}\n")
           outputStream.write(f"p_{pos+1}\n")
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
          