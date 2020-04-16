#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

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
  vegeta.py [options] [-a|-c] <inputSequences> [<genomeOfInterest>]

Options:
  -h, --help                              Show this help message and exits.
  -v, --verbose                           Get some extra information from VeGETA during calculation. [Default: False]
  --version                               Prints the version of VeGETA and exits.
  -o DIR, --output DIR                    Specifies the output directory of VeGETA. [Default: pwd]

  -k KMER, --kmer KMER                    Length of the considered kmer. [Default: 7]
  --cutoff CUTOFF                         Cutoff threshold for the initial graph during clustering. The larger the value the more relationships are
                                          neglected for clustering, despite being closely related. [Default: 0.3]
  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]

  -a, --alignment-only                    Only performs the alignment calculation, without prior clustering. 
                                          NOTE: This is not recommended for large datasets. [Default: False]
  -c, --cluster-only                      Only performs the clustering step of sequences, without the alignment. [Default: False]

  --seedsize SEEDSIZE                     Specifies the length of a region that has to be conserved in order to serve as 
                                          a seed region in the sequence-based scaffold alignment. [Default: 10]
  --shannon SHANNON                       Cut-off value for a seed window based on its averaged shannon entropy.
                                          If none is set, VeGETA takes the best 10% windows as seeds. [Default: 0.1]

  -t THRESHOLD, --tbpp THRESHOLD          Basepairing-probability threshold for potential nucleotide interactions. 
                                          if bpp between two nucleotides is smaller than this value, 
                                          it isn't considered during ILP construction. [Default: 0.7]
  -w WINDOWSIZE, --windowsize WINDOWSIZE  Specifies the window length for the final structure calculation. [Default: 300]
  -s STEPSIZE, --stepsize STEPSIZE        Specifies the step size of the sliding window. [Default: 50]
  --allowLP                               If this is set, VeGETA will include lonely basepairs (isolated helices of length 1)
                                          into the final structure. [Default: False]
  

Version:
  VeGETA v0.1 (alpha)
"""

import sys
import os
import logging
import glob

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
from StructureViruses import StructCalculator

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

    handle.setFormatter(formatter)
    logger.addHandler(handle)
    return logger

def create_outdir(outdir):
    try:
      os.makedirs(outdir)
      os.makedirs(f"{outdir}/tmpSequences")
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
    tbpp = float(d_args['--tbpp'])
  except ValueError:
    logger.error("Invalid number for the basepair probability threshold. Please input a number between 0 and 1.")
    exit(2)
  if not (0 <= tbpp <= 1):
    logger.error("Invalid number for the basepair probability threshold. Please input a number between 0 and 1.")
    exit(2)
  
  try:
    windowSize = int(d_args['--windowsize'])
  except ValueError:
    logger.error("Invalid parameter for the window size. Please input a number.")
    sys.exit(2)
    
  try:
    stepSize = int(d_args['--stepsize'])
  except ValueError:
    logger.error("Invalid parameter for the sliding window step size. Please input a number.")
    sys.exit(2)

  try:
    seedSize = int(d_args['--seedsize'])
  except ValueError:
    logger.error("Invalid parameter for the seed size. Please input a number.")
    sys.exit(2)

  try:
    shannon = float(d_args['--shannon'])
  except ValueError:
    logger.error("Invalid number for the shannon entropy cutoff threshold. Please input a number higher than 0.0.")
    exit(2)

  output = d_args['--output']
  if output == 'pwd':
    output = os.getcwd()
  now = str(datetime.now()).split('.')[0].replace(' ','_').replace(':','-')
  #output = f"{output}/vegeta-{now}"
  output = f"{output}/vegeta/"
  create_outdir(output)


  alnOnly = d_args['--alignment-only']
  clusterOnly = d_args['--cluster-only']
  allowLP = d_args['--allowLP']


  return (inputSequences, goi, output, alnOnly, clusterOnly, k, proc, tbpp, seedSize, windowSize, stepSize, shannon, allowLP)

def perform_clustering():

  multiPool = Pool(processes=proc)
  virusClusterer = Clusterer(logger, inputSequences, k, proc)
  logger.info("Removing 100% identical sequences.")
  virusClusterer.remove_redundancy()
  logger.info("Determining k-mer profiles for all sequences.")
  virusClusterer.determine_profile(multiPool)
  logger.info("Clustering with UMAP and HDBSCAN.")
  virusClusterer.apply_umap(outdir)
  clusterInfo = virusClusterer.clusterlabel
  logger.info(f"Summarized {virusClusterer.dim} sequences into {clusterInfo.max()+1} clusters. Filtered {np.count_nonzero(clusterInfo == -1)} sequences due to uncertainty.")
  logger.info("Extracting centroid sequences and writing results to file.\n")
  virusClusterer.get_centroids(outdir, multiPool)
  virusClusterer.split_centroids(outdir)

  
  
  logger.info(f"Extracting representative sequences for each cluster.")
  sequences = virusClusterer.d_sequences
  distanceMatrix = virusClusterer.matrix
  profiles = virusClusterer.d_profiles
  #virusClusterer.create_subcluster()
  del virusClusterer
  for file in glob.glob(f"{outdir}/cluster*.fa"):
    if file == f"{outdir.rstrip('/')}/cluster-1.fa":
      continue
    virusSubClusterer = Clusterer(logger, file, k, proc, subCluster=True)
    virusSubClusterer.d_sequences = virusSubClusterer.read_sequences()
    code = virusSubClusterer.apply_umap(outdir)
    if code == 1:
      logger.warn(f"Too few sequences for clustering in {os.path.basename(file)}. Alignment will be calculated with all sequences of this cluster.")
      del virusSubClusterer
      continue
    virusSubClusterer.get_centroids(outdir, multiPool)
    virusSubClusterer.split_centroids(outdir)
    del virusSubClusterer

def perform_alignment(seq=None):

  #if seq:
  #  clusteredSequences = seq
  #else:
  #  clusteredSequences = f'{outdir}/representative_viruses.fa'
  logger.info("Starting the alignment step of VeGETA.\n")

  #if goi:
  #  logger.info(f"Including your virus of interest:\n{goi}\n")
  #  with open(clusteredSequences, 'a') as outputStream:
  #    with open(goi, 'r') as inputStream:
  #      outputStream.write("".join(inputStream.readlines()))

  if seq:
    files = [seq]
  else:
    files = glob.glob(f"{outdir}/*_repr.fa")
  
  for file in files:
    prefix = os.path.splitext(os.path.basename(file))[0]
    logger.info(f"Calculating Alignment for {prefix.split('_repr')[0]}")
    virusAligner = Aligner(logger, file, proc, outdir, seedSize, shannon, structureParameter, prefix)
    logger.info("Calculating initial mafft alignment")
    virusAligner.mafft_scaffold()
    logger.info("Finding conserved seeds in the alignment")
    virusAligner.find_seeds_in_scaffold()
    logger.info(f"Found {len(virusAligner.seeds)} seed regions in the alignment")
    logger.info("Extracting sequences between seeds")
    virusAligner.extract_non_seeds()
    logger.info("Applying LocARNA on fragments")
    virusAligner.refine_fragments(windowSize, stepSize)
    logger.info("Merging all fragments to a whole alignment")
    virusAligner.merge_fragments()
    logger.info("Refined alignment calculated.")
    structure = derive_structure(prefix)
    logger.info("Saving the final alignment in STOCKHOLM format")
    write_final_alignment(virusAligner.refinedAlignment, structure, prefix)

def derive_structure(prefix):
  struc = StructCalculator(f"{outdir}/{prefix}_refinedAlignment.aln", logger, outdir, windowSize, stepSize, proc, allowLP, tbpp, prefix)
  struc.apply_alifold()
  struc.calculate_avg_bpp()
  logger.info("Generating ILP based on all basepairing probabilities.")
  logger.info("Solving the ILP may take a while.")
  struc.generate_ilp()
  logger.info("Deriving structural elements from ILP solution.")
  struc.finalize_structure()
  logger.info("Testing individual structural elements for significance (dinucleotide shuffling).")
  return(struc.finalStructure)

def write_final_alignment(alignment, structure, prefix):
  with open(f"{outdir}/{prefix}_finalAlignment.stk",'w') as outputStream:
    outputStream.write("# STOCKHOLM 1.0\n")
    outputStream.write("#=GF AU  Kevin Lamkiewicz\n")
    outputStream.write("#=GF BM  VeGETA v. 0.1\n")
    outputStream.write(f"#=GF SQ  {len(alignment)}\n\n")
    
    for record in alignment:
    #for header, sequence in alignment.items():
      outputStream.write(f"{record.id}\t\t{record.seq}\n")
    outputStream.write(f"#=GC SS_cons\t\t{structure}\n")

  #virusAligner.calculate_pw_distances()
  #virusAligner.get_tree_from_dist()
  #treePath = f"{os.path.dirname(outdir)}/test_tree.nwk"
  #Phylo.write(virusAligner.tree, treePath, 'newick')
  #virusAligner.refine_pairwise_instances(virusAligner.tree, None)
  

if __name__ == "__main__":
  logger = create_logger()
  (inputSequences, goi, outdir, alnOnly, clusterOnly, k, proc, tbpp, seedSize, windowSize, stepSize, shannon, allowLP) = parse_arguments(docopt(__doc__))

  structureParameter = (logger, outdir, windowSize, stepSize, proc, allowLP, tbpp)

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
