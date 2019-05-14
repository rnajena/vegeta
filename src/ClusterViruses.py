#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
from multiprocessing import Pool
import itertools
import math
import numpy as np
import os

class Clusterer(object):
  """
  """

  d_profiles = {}
  matrix = {}
  id2header = {}
  header2id = {}
  mclMatrix = ''
  allCluster = []

  def __init__(self, sequenceFile, k, cutoff):
    self.sequenceFile = sequenceFile
    self.k = k
    self.cutoff = cutoff
    self.nucleotides = set(["A","C","G","T"])
    self.allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(self.nucleotides, repeat=self.k))}

    self.d_sequences = self.read_sequences()
    

  def read_sequences(self):
    """
    """
    fastaContent = {}
    idHead = 0
    with open(self.sequenceFile, 'r') as inputStream:
      header = ''
      seq = ''

      for line in inputStream:
        if line.startswith(">"):
          if header:
            fastaContent[header] = seq
            self.id2header[str(idHead)] = header
            self.header2id[header] = str(idHead)
          header = line.rstrip("\n").replace(':','_').replace(' ','_').lstrip(">")
          seq = ''
          idHead += 1
        else:
          seq += line.rstrip("\n").upper().replace('U','T')
    return fastaContent


  def determine_profile(self):
    """
    """
    for header, sequence in self.d_sequences.items():    
      self.d_profiles[header] = [0]*len(self.allKmers)
      kmer = [sequence[start : start + self.k] for start in range(len(sequence) - self.k)]
      for k in kmer:
        try:
          self.d_profiles[header][self.allKmers[k]] += 1 
        except KeyError:
          continue
      self.matrix[header] = []
    

  def calc_pd(self, seqs):
    seq1, seq2 = seqs
    
    profile1 = np.array(self.d_profiles[seq1])
    profile2 = np.array(self.d_profiles[seq2])
    distance = np.sqrt(np.sum((profile1 - profile2)**2))
      
    return {seq1 : (seq2, distance), seq2 : (seq1, distance)}


  def pairwise_distances(self, proc):
    """
    """
    p = Pool(processes=proc)
    distances = (p.map(self.calc_pd, itertools.combinations(self.d_profiles, 2)))
    
    for pairwiseDist in distances:
      for first,second in pairwiseDist.items():
        self.matrix[first].append(second)
      
  def normalize_function(self):
    """
    """
    allDistances = []
    
    for row in self.matrix.values():
      for _, distance in row:
        allDistances.append(distance)
    
    def normalize(x):
      return 1-((x - min(allDistances)) / (max(allDistances) - min(allDistances)))

    return normalize


  def create_matrix(self):
    """
    """
    self.mclMatrix = f"(mclheader\nmcltype matrix\ndimensions {len(self.d_sequences)}x{len(self.d_sequences)}\n)\n"
    self.mclMatrix += f"(mcldoms\n{' '.join(self.id2header)} $\n)\n"
    self.mclMatrix += f"(mclmatrix\nbegin\n"
    normalize = self.normalize_function()
  
    for header, distances in self.matrix.items():
      row = self.header2id[header]
      consideredDistances = ' '.join([f'{self.header2id[node]}:{normalize(dist)}' for node, dist in distances if normalize(dist) > self.cutoff])
      
      self.mclMatrix += f"{row} {consideredDistances} $\n"
    self.mclMatrix += ")\n"
      
  def perform_mcl(self, outdir):
    """
    """
    with open(f"{outdir}/mclInput.txt", 'w') as outputStream:
      outputStream.write(self.mclMatrix)
    os.system(f"mcl {outdir}/mclInput.txt -I 9 -o {outdir}/mclClustered.txt 2>/dev/null")
    
  def extract_cluster(self, outdir):
    """
    """
    with open(f'{outdir}/mclClustered.txt', 'r') as inputStream:
      while not inputStream.readline().startswith('begin'):
        continue

      newCluster = []

      for line in inputStream:
        if line.startswith(')'):
          self.allCluster.append(newCluster)
          break
  
        cluster = line.rstrip("$\n").split(' ')
        line = line.rstrip("$\n").split()

        if line[-1] == '$':
          line = line[:-1]

        if cluster[0]:
          if newCluster:
            self.allCluster.append(newCluster)
          newCluster = [self.id2header[x] for x in line[1:]]
        else:
          newCluster.extend([self.id2header[x] for x in line])

    #print(self.allCluster)
    

  def get_centroids(self):
    """
    """
    pass

