#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from collections import Counter
from multiprocessing import Pool
import itertools
import math
import numpy as np
import os

import scipy
import random

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
  dim = 0
  probabilities = []

  def __init__(self, sequenceFile, k, cutoff, proc):
    self.sequenceFile = sequenceFile
    self.k = k
    self.cutoff = cutoff
    self.nucleotides = set(["A","C","G","T"])
    self.allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(self.nucleotides, repeat=self.k))}
    #self.proc = proc
    

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

      fastaContent[idHead] = seq
      self.id2header[idHead] = header
      self.header2id[header] = idHead
    return fastaContent

  def profile(self, entry):
    header, sequence = entry
    profile = [0]*len(self.allKmers)
    for k in iter([sequence[start : start + self.k] for start in range(len(sequence) - self.k)]):
        try:
          profile[self.allKmers[k]] += 1 
        except KeyError:
          continue
    return (header, profile)

  def determine_profile(self, proc):
    p = Pool(proc)
    for header, profile in p.map(self.profile, self.d_sequences.items()):
      self.d_profiles[header] = profile
    p.close()
    p.join()
    
  def calc_pd(self, seqs):
    seq1, seq2 = seqs
    
    profile1 = np.array(self.d_profiles[seq1])
    profile2 = np.array(self.d_profiles[seq2])

    distance = scipy.spatial.distance.cosine(profile1, profile2)
    return (seq1, seq2, distance)

  #de#f pairwise_distances(self):
    #"""
    #p = Pool(processes=self.proc)
    #for seq1, seq2, dist in p.map(self.calc_pd, itertools.combinations(self.d_profiles, 2)):
    #  self.matrix[seq1][seq2] = dist
    #  self.matrix[seq2][seq1] = dist  
    #p.close()
      
  def get_centroids(self, outdir, proc):
    """
    """
    centroids = []
    seqCluster = { x : [] for x in set(self.allCluster)}
    

    for idx, cluster in enumerate(self.allCluster):
      seqCluster[cluster].append(idx)

    p = Pool(proc)    
    for cluster, sequences in seqCluster.items():
      # the -1 cluster from HDBSCAN is unclassified and thus ignored
      if cluster == -1:
        continue
      subProfiles = {seq : profile for seq,profile in self.d_profiles.items() if seq in sequences}

      
      for seq1, seq2, dist in p.map(self.calc_pd, itertools.combinations(subProfiles, 2)):
        self.matrix[seq1][seq2] = dist
        self.matrix[seq2][seq1] = dist
      
      

      tmpMinimum = math.inf
      centroidOfCluster = -1
      if len(sequences) == 1:
        centroidOfCluster = cluster[0]
        centroids.append(centroidOfCluster)
        continue
      for sequence in sequences:
        averagedDistance = 0
        for neighborSequence in sequences:
          if sequence == neighborSequence:
            continue
          averagedDistance += self.matrix[sequence][neighborSequence]
        averagedDistance /= len(sequences)-1
        if averagedDistance < tmpMinimum:
          tmpMinimum = averagedDistance
          centroidOfCluster = sequence
        
      centroids.append(centroidOfCluster)
    p.close()
    p.join()
    with open(f'{outdir}/representative_viruses.fa', 'w') as outStream:
      for centroid in centroids:
        outStream.write(f">{self.id2header[centroid]}\n{self.d_sequences[centroid]}\n")

  def apply_umap(self):
    profiles = []
    for idx, _ in enumerate(self.d_profiles):
      profiles.append(self.d_profiles[idx])
  
    clusterable_embedding = umap.UMAP(
          n_neighbors=20,
          min_dist=0.0,
          n_components=10,
          random_state=42,
      ).fit_transform(profiles)
    
    clusterer = hdbscan.HDBSCAN()
    clusterer.fit(clusterable_embedding)

    self.allCluster = clusterer.labels_
    self.probabilities = clusterer.probabilities_