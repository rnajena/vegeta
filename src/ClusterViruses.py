#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
from multiprocessing import Pool
import itertools
import math
import numpy as np
import os

import umap.umap_ as umap
import hdbscan

class Clusterer(object):
  """
  """

  d_profiles = {} 
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
    self.dim = len(self.d_sequences)
    self.matrix = np.zeros(shape=(self.dim, self.dim),dtype=float)

  def read_sequences(self):
    """
    """
    fastaContent = {}
    idHead = -1
    uniqueSeqs = set()
    with open(self.sequenceFile, 'r') as inputStream:
      header = ''
      seq = ''

      for line in inputStream:
        if line.startswith(">"):
          if header:
            if not seq in uniqueSeqs:
              fastaContent[idHead] = seq
              self.id2header[idHead] = header
              self.header2id[header] = idHead
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
    
  def calc_pd(self, seqs):
    seq1, seq2 = seqs
    
    profile1 = np.array(self.d_profiles[seq1])
    profile2 = np.array(self.d_profiles[seq2])
    distance = np.sqrt(np.sum((profile1 - profile2)**2))
    return (seq1, seq2, distance)

  def pairwise_distances(self, proc):
    """
    """
    p = Pool(processes=proc)
    for seq1, seq2, dist in p.map(self.calc_pd, itertools.combinations(self.d_profiles, 2)):
      self.matrix[seq1][seq2] = dist
      self.matrix[seq2][seq1] = dist

    self.matrix[self.matrix == 0] = np.nan
    normalize = self.normalize_function()
    self.normalize_matrix(normalize)
      
  def normalize_function(self):
    """
    """
    
    minimum = np.nanmin(self.matrix)
    maximum = np.nanmax(self.matrix)
    #cutoff = self.cutoff

    def normalize(x):
      normalizedDist = 1-((x - minimum) / (maximum - minimum))
      #return normalizedDist
      return normalizedDist if normalizedDist != np.nan else 1
      
    return normalize

  def normalize_matrix(self,f):
    """
    """
    newMatrix =  np.zeros(shape=(self.dim, self.dim),dtype=float)
    for i, row in enumerate(self.matrix):
      for j, _ in enumerate(row):
        if i != j:
          newMatrix[i][j] = f(self.matrix[i][j])
    self.matrix = newMatrix    


  def create_matrix(self):
    """
    """
    self.mclMatrix = f"(mclheader\nmcltype matrix\ndimensions {len(self.d_sequences)}x{len(self.d_sequences)}\n)\n"
    self.mclMatrix += f"(mcldoms\n{' '.join(map(str,self.id2header))} $\n)\n"
    self.mclMatrix += f"(mclmatrix\nbegin\n"
    
    for row, entries in enumerate(self.matrix):
      consideredDistances = ' '.join([f'{col}:{dist}' for col, dist in enumerate(entries) if row != col and dist != 0.0])  
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
        line = list(map(int, line))

        if cluster[0]:
          if newCluster:
            self.allCluster.append(newCluster)
          newCluster = line[1:]
        else:
          newCluster.extend(line)

  def get_centroids(self, outdir):
    """
    """
    centroids = []
    for cluster in self.allCluster:

      tmpMinimum = 5
      centroidOfCluster = -1

      if len(cluster) == 1:
        centroidOfCluster = cluster[0]
        centroids.append(centroidOfCluster)
        break

      for sequence in cluster:
        averagedDistance = 0
        for neighborSequence in cluster:
          if sequence == neighborSequence: 
            continue
          averagedDistance += self.matrix[sequence][neighborSequence]
        averagedDistance /= len(cluster)-1

    
        if averagedDistance < tmpMinimum:
          tmpMinimum = averagedDistance
          centroidOfCluster = sequence
        
      centroids.append(centroidOfCluster)
    
    with open(f'{outdir}/representative_viruses.fa', 'w') as outStream:
      for centroid in centroids:
        outStream.write(f">{self.id2header[centroid]}\n{self.d_sequences[centroid]}\n")


  def apply_umap(self):
    clusterable_embedding = umap.UMAP(
          n_neighbors=30,
          min_dist=0.0,
          n_components=10,
          random_state=42,
      ).fit_transform(self.matrix)
    
    clusterer = hdbscan.HDBSCAN()
    clusterer.fit(clusterable_embedding)

    print(clusterer.labels_)