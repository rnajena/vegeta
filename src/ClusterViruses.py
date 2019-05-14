#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
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

  def __init__(self, sequenceFile, k):
    self.sequenceFile = sequenceFile
    self.k = k
    nucleotides = ["A","C","G","T"]
    self.allKmers = [''.join(kmer) for kmer in itertools.product(nucleotides, repeat=self.k)]
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
      self.d_profiles[header] = {kmer : 0 for kmer in self.allKmers}
      kmer = [sequence[start : start + self.k] for start in range(len(sequence) - self.k)]
      self.d_profiles[header].update(Counter(kmer))
      self.matrix[header] = {}
    

  def pairwise_distances(self):
    """
    """
    for seq1, seq2 in itertools.combinations(self.d_profiles, 2):
      profile1 = np.array([self.d_profiles[seq1][kmer] for kmer in self.allKmers])
      profile2 = np.array([self.d_profiles[seq2][kmer] for kmer in self.allKmers])
      
      distance = np.sqrt(np.sum((profile1 - profile2)**2))
      
      self.matrix[seq1].update({seq2:distance})
      self.matrix[seq2].update({seq1:distance})
      
  def create_matrix(self):
    """
    """
    self.mclMatrix = f"(mclheader\nmcltype matrix\ndimensions {len(self.d_sequences)}x{len(self.d_sequences)}\n)\n"
    self.mclMatrix += f"(mcldoms\n{' '.join(self.id2header)} $\n)\n"
    self.mclMatrix += f"(mclmatrix\nbegin\n"
    
    for header, distances in self.matrix.items():
      self.mclMatrix += f"{self.header2id[header]} {' '.join([f'{self.header2id[node]}:{dist}' for node, dist in distances.items()])} $\n"
    self.mclMatrix += ")\n"
      
  def perform_mcl(self, outdir):
    """
    """
    with open(f"{outdir}/mclInput.txt", 'w') as outputStream:
      outputStream.write(self.mclMatrix)
    os.system(f"mcl {outdir}/mclInput.txt -I 5 -o {outdir}/mclClustered.txt 2>/dev/null")
    
  def extract_cluster(self, outdir):
    """
    """
    startParsing = False
    with open(f'{outdir}/mclClustered.txt', 'r') as inputStream:
      for line in inputStream:

        if startParsing:
          if line.startswith(')'):
            break
          cluster = line.rstrip("$\n").split()
          header = [self.id2header[x] for x in cluster[1:]]
          self.allCluster.append(header)

        if line.startswith("begin"):
          startParsing = True
          continue   
    

  def get_centroids(self):
    """
    """
    pass

