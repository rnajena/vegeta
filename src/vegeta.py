#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VeGETA -- Viral GEnome sTructure Alignments

VeGETA is a python program that takes several genome sequences
from different viruses as an input. In the very first step,
VeGETA will cluster these sequences into groups (clades) based 
on their sequence similarity. For each clade, the centroid sequence is
determined as representative genome, i.e. the sequence with the lowest
distance to all other sequences of this clade. 

In a second step, VeGETA calculates a multiple sequence alignment of 
the representative genomes. This alignment is progressively refined
using 'MAFFT' and 'LocARNA'. In the end, we provide a local-structure guided
MSA of the representative viruses.

Python Dependencies:
  docopt

Other Dependencies:
  ViennaRNA package
  MAFFT
  LocARNA

Contact:
  kevin.lamkiewicz@uni-jena.de

Citation:
  Lamkiewicz, K., et al. (20xx), "Structure-guided multiple sequence alignments of complete viral genomes.", Journal, Issue, Volume.

Usage:
  vegeta.py [options] <inputSequences> [<genomeOfInterest>]

Options:
  -h, --help                          Show this help message and exits.
  -v, --verbose                       Get some extra information from VeGETA during calculation. [Default: False]
  --version                           Prints the version of VeGETA and exits.
  -o DIR, --output DIR                Specifies the output directory of VeGETA. [Default: pwd]

  -k KMER, --kmer KMER                Length of the considered kmer. [Default: 7]
  --cutoff CUTOFF                     Cutoff threshold for the initial graph during clustering. The larger the value the more relationships are
                                      neglected for clustering, despite being closely related. [Default: 0.3]
  -p PROCESSES, --process PROCESSES   Specify the number of CPU cores that are used. [Default: 1]

  -a, --alignment-only                Only performs the alignment calculation, without prior clustering. 
                                        NOTE: This is not recommended for large datasets. [Default: False]
  -c, --cluster-only                  Only performs the clustering step of sequences, without the alignment. [Default: False]
  

Version:
  VeGETA v0.1 (alpha)
"""

import sys
import os
import logging

import numpy as np
from multiprocessing import Pool

from datetime import datetime

from colorlog import ColoredFormatter
from docopt import docopt
from Bio import Phylo

inputSequences = None
goi = None
outdir = None
alnOnly = None
clusterOnly = None
k = None
proc = None
cutoff = None

def warn(*args, **kwargs):
    pass
import warnings
warnings.warn = warn

from ClusterViruses import Clusterer
from AlignViruses import Aligner

def create_logger():
    """
    doc string.
    """

    logger = logging.getLogger()
    #logger.setLevel(logging.WARNING)

    handle = logging.StreamHandler()
    #handle.setLevel(logging.WARNING)

    formatter = ColoredFormatter("%(log_color)sVeGETA %(levelname)s -- %(asctime)s -- %(message)s", "%Y-%m-%d %H:%M:%S", 
                                    log_colors={
                                            'DEBUG':    'bold_cyan',
                                            'INFO':     'bold_white',
                                            'WARNING':  'bold_yellow',
                                            'ERROR':    'bold_red',
                                            'CRITICAL': 'bold_red'}
                                )

    #handle.setFormatter(logging.Formatter("ViMiFi %(levelname)s -- %(asctime)s -- %(message)s", "%Y-%m-%d %H:%M:%S"))
    handle.setFormatter(formatter)
    logger.addHandler(handle)
    return logger

def create_outdir(outdir):
    try:
      os.makedirs(outdir)
      logger.info(f"Creating output directory: {outdir}")
    except FileExistsError:
      logger.warning(f"The output directory exists. Files will be overwritten.")

def parse_arguments(d_args):
  """
  Parse all given arguments and check for error (e.g. file names).

  Arguments:
  d_args -- dict with input parameters and their values

  Returns:
  parsed and checked parameter list
  """

  if d_args['--version']:
    print("VeGETA version 0.1")
    exit(0)

  
  verbose = d_args['--verbose']
  if verbose:
    logger.setLevel(logging.INFO)

  inputSequences = d_args['<inputSequences>']
  if not os.path.isfile(inputSequences):
    logger.error("Couldn't find input sequences. Check your file")
    exit(1)

  goi = d_args['<genomeOfInterest>']
  if goi and not os.path.isfile(goi):
    logger.error("Couldn't find genome of interest. Check your file")
    exit(1)

  try:
    k = int(d_args['--kmer'])
  except ValueError:
    logger.error("Invalid parameter for k-mer size. Please input a number.")
    exit(2)
  
  try:
    proc = int(d_args['--process'])
  except ValueError:
    logger.error("Invalid number for CPU cores. Please input a number.")
    exit(2)

  try:
    cutoff = float(d_args['--cutoff'])
  except ValueError:
    logger.error("Invalid number for the cutoff threshold. Please input a number between 0 and 1.")
    exit(2)
  if not (0 <= cutoff <= 1):
    logger.error("Invalid number for the cutoff threshold. Please input a number between 0 and 1.")
    exit(2)
    
  output = d_args['--output']
  if output == 'pwd':
    output = os.getcwd()
  now = str(datetime.now()).split('.')[0].replace(' ','_').replace(':','-')
  #output = f"{output}/vegeta-{now}"
  output = f"{output}/vegeta"
  create_outdir(output)


  alnOnly = d_args['--alignment-only']
  clusterOnly = d_args['--cluster-only']


  return (inputSequences, goi, output, alnOnly, clusterOnly, k, proc, cutoff)

def perform_clustering():

  multiPool = Pool(processes=proc)
  virusClusterer = Clusterer(inputSequences, k, cutoff, proc)
  logger.info("Determining k-mer profiles for all sequences.")
  virusClusterer.determine_profile(multiPool)
  logger.info("Clustering with UMAP and HDBSCAN.")
  virusClusterer.apply_umap()
  clusterInfo = virusClusterer.allCluster
  logger.info(f"Summarized {virusClusterer.dim} sequences into {clusterInfo.max()+1} clusters. Filtered {np.count_nonzero(clusterInfo == -1)} sequences due to uncertainty.")
  logger.info("Extracting centroid sequences and writing results to file.\n")
  virusClusterer.get_centroids(outdir, multiPool)
  virusClusterer.split_centroids(outdir)

def perform_alignment(seq=None):
  if seq:
    clusteredSequences = seq
  else:
    clusteredSequences = f'{outdir}/representative_viruses.fa'
  logger.info("Starting the alignment step of VeGETA.\n")

  if goi:
    logger.info(f"Including your virus of interest:\n{goi}\n")
    with open(clusteredSequences, 'a') as outputStream:
      with open(goi, 'r') as inputStream:
        outputStream.write("".join(inputStream.readlines()))

  
  virusAligner = Aligner(logger, clusteredSequences, k, proc)
  logger.info("Calculating initial mafft alignment")
  virusAligner.calculate_pw_distances()
  virusAligner.get_tree_from_dist()
  treePath = f"{os.path.dirname(outdir)}/test_tree.nwk"
  Phylo.write(virusAligner.tree, treePath, 'newick')
  virusAligner.refine_pairwise_instances(virusAligner.tree, None)
  

if __name__ == "__main__":
  logger = create_logger()
  (inputSequences, goi, outdir, alnOnly, clusterOnly, k, proc, cutoff) = parse_arguments(docopt(__doc__))

  if alnOnly:
    logger.info("Skipping clustering and directly calculate the alignment.")
    perform_alignment(seq=inputSequences)
  elif clusterOnly:
    logger.info("Only clustering is performed. The alignment calculation will be skipped.")
    perform_clustering()
  else:
    logger.info("Doing both, the clustering and the alignment step.")
    perform_clustering()
    perform_alignment()



  """
  Clustering of input sequences done.
  """

    

    