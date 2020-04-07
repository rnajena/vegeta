#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

from collections import Counter
from multiprocessing import Pool
import multiprocessing as mp
import itertools
import math
import numpy as np
import os
import subprocess
import re

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
  allCluster = []
  dim = 0
  probabilities = []



  codon2aminoacid = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W', 
    } 

  def __init__(self, logger, sequenceFile, k, cutoff, proc):
    """
    """


    #self.logger = logger
    self.sequenceFile = sequenceFile
    self.reducedSequences = f"{os.path.realpath(os.path.dirname(self.sequenceFile))}/reduced.fasta"
    self.k = k
    self.cutoff = cutoff
    self.nucleotides = set(["A","C","G","T"])
    self.allKmers = {''.join(kmer):x for x,kmer in enumerate(itertools.product(self.nucleotides, repeat=self.k))}
    self.proc = proc
    self.d_comp = {"A" : "T", "C" : "G", "G" : "C", "T" : "A"}
    self.d_sequences = {}
    self.dim = 0
    
    self.centroids = []
    #TODO: der regex muss lazy evaluieren... see whether it does!
    #self.regex_orf = re.compile(r'((ATG|CTG)(...){25,}?(TAG|TAA|TGA))')
    self.regex_orf = re.compile(r'M[^*]{25,}?\*')

  def rev_comp(self, sequence):
    """
    """
    return ("".join([self.d_comp[nt] if nt in self.d_comp else nt for nt in sequence[::-1]]))
  
  def remove_redundancy(self):
    """
    """
    cmd = f"cd-hit-est -i {self.sequenceFile} -o {self.reducedSequences} -c 1"
    TRASH = open(os.devnull, 'w')
    subprocess.run(cmd.split(), check=True, stderr=TRASH, stdout=TRASH)
    self.d_sequences = self.read_sequences()
    self.dim = len(self.d_sequences)
    self.matrix = np.zeros(shape=(self.dim, self.dim),dtype=float)

  def read_sequences(self):
    """
    """
    fastaContent = {}
    idHead = -1
    
    with open(self.reducedSequences, 'r') as inputStream:
      header = ''
      seq = ''

      for line in inputStream:
        if line.startswith(">"):
          if header:
            seq = seq + "X"*10 + self.rev_comp(seq)
            #if not seq in uniqueSeqs:
            fastaContent[idHead] = seq
            self.id2header[idHead] = header
            self.header2id[header] = idHead
          header = line.rstrip("\n").replace(':','_').replace(' ','_').lstrip(">")
          seq = ''
          idHead += 1
        else:
          seq += line.rstrip("\n").upper().replace('U','T')

      #seq = seq + "X"*10 + self.rev_comp(seq)
      fastaContent[idHead] = seq
      self.id2header[idHead] = header
      self.header2id[header] = idHead
    return fastaContent

  def profile(self, entry):
    """
    """
    header, sequence = entry
    profile = [0]*len(self.allKmers)
    for k in iter([sequence[start : start + self.k] for start in range(len(sequence) - self.k)]):
        try:
          profile[self.allKmers[k]] += 1 
        except KeyError:
          continue
    kmerSum = sum(profile)
    profile = list(map(lambda x: x/kmerSum, profile))
    return (header, profile)

  def determine_profile(self, proc):
    p = Pool(self.proc)
    allProfiles = p.map(self.profile, self.d_sequences.items())
    p.close()
    p.join()
    for header, profile in allProfiles:
    #for header, profile in p.map(self.profile, self.d_sequences.items()):
      self.d_profiles[header] = profile

  def calc_pd(self, seqs):
    #seq1, seq2, subProfiles = seqs
    
    
    #profile1 = np.array(subProfiles[seq1])
    #profile2 = np.array(subProfiles[seq2])
    
    for element in seqs:
      try:
        stuff = (element[0], element[1])
      except TypeError:
        #print(element)
        #print(seqs)
        return None
    #return (0)
    
    seq1, profile1 = seqs[0]
    seq2, profile2 = seqs[1]
    #except TypeError:
    #  print(seqs)
    #  return (-1, -1, 0)

    distance = scipy.spatial.distance.cosine(profile1, profile2)
    #return distance
    return (seq1, seq2, distance)
    #       

  def apply_umap(self, outdir):
    """
    """
    profiles = []
    for idx, _ in enumerate(self.d_profiles):
      profiles.append(self.d_profiles[idx])
  
    clusterable_embedding = umap.UMAP(
          n_neighbors=20,
          min_dist=0.0,
          n_components=10,
          random_state=42,
          metric='cosine',
      ).fit_transform(profiles)
    
    clusterer = hdbscan.HDBSCAN()
    clusterer.fit(clusterable_embedding)

    self.allCluster = clusterer.labels_
    self.probabilities = clusterer.probabilities_

    print()
    print()

    with open(f'{outdir}/cluster.txt', 'w') as outStream:
      for i in set(self.allCluster):
        with open(f'{outdir}/cluster{i}.fa', 'w') as fastaOut:
          outStream.write(f"Cluster: {i}\n")
          for idx, label in enumerate(self.allCluster):
            if label == i:
              outStream.write(f"{self.id2header[idx]}\n")
              fastaOut.write(f">{self.id2header[idx]}\n{self.d_sequences[idx]}\n")
        outStream.write("\n")


  def get_centroids(self, outdir, proc):
    """
    """
    
    seqCluster = { x : [] for x in set(self.allCluster)}

    for idx, cluster in enumerate(self.allCluster):
      seqCluster[cluster].append(idx)

    p = Pool(self.proc)

    for cluster, sequences in seqCluster.items():
      if cluster == -1:
        continue

      subProfiles = {seq : profile for seq, profile in self.d_profiles.items() if seq in sequences}


      for result in p.map(self.calc_pd, itertools.combinations(subProfiles.items(), 2)):
    #for seq1, seq2 in itertools.combinations(subProfiles, 2):
        #print(result)
        seq1, seq2, dist = result
        #dist = self.calc_pd((seq1, seq2, subProfiles))
        self.matrix[seq1][seq2] = dist
        self.matrix[seq2][seq1] = dist
      
      #self.logger.info("Distances are calculated!")
      tmpMinimum = math.inf
      centroidOfCluster = -1

      if len(sequences) == 1:
        centroidOfCluster = cluster[0]
        self.centroids.append(centroidOfCluster)
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

      self.centroids.append(centroidOfCluster)
    #self.centroids = p.map(self.centroid_from_cluster, seqCluster.items())
    p.close()
    p.join()


  def split_centroids(self, outdir):
    """
    """
    centroids = { centroid : self.d_sequences[centroid].split("X"*10) for centroid in self.centroids }

     

    # with open(f'{outdir}/representative_viruses.fa', 'w') as outStream:
    #   for centroidID, sequence in centroids.items():
    #     outStream.write(f">{self.id2header[centroidID]}\n{sequence}\n")
    reprSeqs = {}
    for centroidID, strands in centroids.items():
      positiveStrand = ""
      longestCDS = 0
      for strand in strands:
        for frame in range(3):
          proteinSequence = ""
          for fragment in range(frame, len(strand), 3):
            codon = strand[fragment:fragment+3]
            if len(codon) != 3:
              continue
            try:
              proteinSequence += self.codon2aminoacid[codon]
            except KeyError:
              proteinSequence += 'X'
          matches = self.regex_orf.findall(proteinSequence)
          allORFs = "".join([x for x in matches if x])
          if len(allORFs) > longestCDS:
            longestCDS = len(allORFs)
            positiveStrand = strand
      reprSeqs[centroidID] = positiveStrand 
    
    with open(f'{outdir}/representative_viruses.fa', 'w') as outStream:
      for centroidID, sequence in reprSeqs.items():
        outStream.write(f">{self.id2header[centroidID]}\n{sequence}\n")
      

      #for strand in strands:
      #  matches = self.regex_orf.findall(strand)
        #matches = [self.regex_orf.findall(strand) for strand in strands]
        
      #  print(centroidID, len(matches))
      #  allORFs = "".join([x[0] for x in matches if x])
      #  print(len(allORFs)/len(strand))
      
      #print(f"Start Codon {[x.count('ATG') for x in strands]}")
      #print(f"Stop Codon {[len(self.regex_stop.findall(strand)) for strand in strands]}")
    