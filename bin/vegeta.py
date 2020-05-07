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
  BioPython
  colorlog
  numpy
  scipy
  umap-learn
  hdbscan

Other Dependencies:
  ViennaRNA package
  cd-hit
  MAFFT
  LocARNA
  GLPK

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
  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]

  -a, --alignment-only                    Only performs the alignment calculation, without prior clustering. 
                                          NOTE: This is not recommended for large datasets. [Default: False]
  -c, --cluster-only                      Only performs the clustering step of sequences, without the alignment. [Default: False]


  --subcluster                            Additionally to the initial alignment, each cluster gets analyzed for 
                                          local structures and relations which results in an alignment for each initial 
                                          cluster. WARNING: This will increase the overall runtime drastically! [Default: False]
  --seedsize SEEDSIZE                     Specifies the length of a region that has to be conserved in order to serve as 
                                          a seed region in the sequence-based scaffold alignment. [Default: 10]
  --shannon SHANNON                       Cut-off value for a seed window based on its averaged shannon entropy.
                                          If none is set, VeGETA takes the best 10% windows as seeds. [Default: 0.1]

  -t THRESHOLD, --tbpp THRESHOLD          Basepairing-probability threshold for potential nucleotide interactions. 
                                          if bpp between two nucleotides is smaller than this value, 
                                          it isn't considered during ILP construction. [Default: 0.7]
  -w WINDOWSIZE, --windowsize WINDOWSIZE  Specifies the window length for the final structure calculation. [Default: 300]
  -s STEPSIZE, --stepsize STEPSIZE        Specifies the step size of the sliding window. [Default: 50]

  --allowLP                               NOT IMPLEMENTED COMPLETELY YET -- If this is set, VeGETA will include lonely basepairs (isolated helices of length 1)
                                          into the final structure. [Default: False]
  
  --shuffle SHUFFLE                       Number of sampled sequences created for each structural element; required for 
                                          significance testing (z-score analyses). [Default: 500]
  --pvalue PVALUE                         p-Value threshold whether a structure gets accepted or not in the final solution
                                          based on its z-score analyses. [Default: 0.05]

Version:
  VeGETA v0.3 (alpha)
"""

import sys
import os
import logging
import glob
import shutil

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
    print("VeGETA version 0.3")
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

  try:
    shuffle = int(d_args['--shuffle'])
  except ValueError:
    logger.error("Invalid number for the number of shuffle events. Please input a number higher than 0.")
    exit(2)

  try:
    pvalue = float(d_args['--pvalue'])
  except ValueError:
    logger.error("Invalid number for the p-value threshold. Please input a number between 0.0 and 1.0")
    exit(2)

  if not (0 <= pvalue <= 1):
    logger.error("Invalid number for the p-value threshold. Please input a number between 0.0 and 1.0")
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
  subcluster = d_args['--subcluster']

  return (inputSequences, goi, output, alnOnly, clusterOnly, k, proc, tbpp, subcluster, seedSize, windowSize, stepSize, shannon, allowLP, shuffle, pvalue)

def __abort_cluster(clusterObject, filename):
    logger.warn(f"Too few sequences for clustering in {os.path.basename(filename)}. Alignment will be calculated with all sequences of this cluster.")
    del clusterObject

def perform_clustering():

  multiPool = Pool(processes=proc)
  virusClusterer = Clusterer(logger, inputSequences, k, proc, outdir, goi=goi)

  logger.info("Removing 100% identical sequences.")
  code = virusClusterer.remove_redundancy()
  logger.info("Sequences are all parsed.")

  if code == 1:
    __abort_cluster(virusClusterer, inputSequences)
    return 0

  if goi:
    logger.info(f"Found {len(virusClusterer.goiHeader)} genome(s) of interest.")
  logger.info("Determining k-mer profiles for all sequences.")
  virusClusterer.determine_profile(multiPool)
  logger.info("Clustering with UMAP and HDBSCAN.")
  code = virusClusterer.apply_umap()
  if code == 1:
    __abort_cluster(virusClusterer, inputSequences)
    #logger.warning(f"All sequences fall into one cluster. Aligning this one without dividing the sequences anymore.")
    return 0
  clusterInfo = virusClusterer.clusterlabel
  logger.info(f"Summarized {virusClusterer.dim} sequences into {clusterInfo.max()+1} clusters. Filtered {np.count_nonzero(clusterInfo == -1)} sequences due to uncertainty.")

  goiCluster = virusClusterer.goi2Cluster
  if goiCluster:
    for header, cluster in goiCluster.items():
      logger.info(f"You find the genome {header} in cluster {cluster}.")

  logger.info("Extracting centroid sequences and writing results to file.\n")
  virusClusterer.get_centroids(multiPool)
  virusClusterer.split_centroids()
  
  logger.info(f"Extracting representative sequences for each cluster.")
  sequences = virusClusterer.d_sequences
  distanceMatrix = virusClusterer.matrix
  profiles = virusClusterer.d_profiles
  del virusClusterer

  if not subcluster:
    return 0

  for file in glob.glob(f"{outdir}/cluster*.fa"):
    if file == f"{outdir.rstrip('/')}/cluster-1.fa":
      continue
    virusSubClusterer = Clusterer(logger, file, k, proc, outdir, subCluster=True)
    code = virusSubClusterer.remove_redundancy()
    
    if code == 1:
      __abort_cluster(virusSubClusterer, file)
      #logger.warn(f"Too few sequences for clustering in {os.path.basename(file)}. Alignment will be calculated with all sequences of this cluster.")
      #del virusSubClusterer
      continue

    code = virusSubClusterer.apply_umap()
    
    if code == 1:
      __abort_cluster(virusSubClusterer, file)
      #logger.warn(f"Too few sequences for clustering in {os.path.basename(file)}. Alignment will be calculated with all sequences of this cluster.")
      #del virusSubClusterer
      continue

    virusSubClusterer.get_centroids(multiPool)
    virusSubClusterer.split_centroids()
    del virusSubClusterer

def perform_alignment(seq=None):

  #if seq:
  #  clusteredSequences = seq\
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
    
    try:
      os.makedirs(f"{outdir}/tmpSequences")
    except FileExistsError:
      pass # I always wanted to do this; except-pass == high-quality code

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
    shutil.rmtree(f"{outdir}/tmpSequences")

def derive_structure(prefix):
  struc = StructCalculator(f"{outdir}/{prefix}_refinedAlignment.aln", logger, outdir, windowSize, stepSize, proc, allowLP, tbpp, prefix, shuffle, pvalue)
  logger.info("Applying RNAalifold on alignment windows.")
  struc.apply_alifold()
  logger.info("Parsing basepairing probabilities out of windows.")
  struc.calculate_avg_bpp()
  logger.info("Generating ILP based on all basepairing probabilities.")
  logger.info("Solving the ILP may take a while.")
  struc.generate_ilp()
  logger.info("Deriving structural elements from ILP solution and testing individual structural elements for significance (nucleotide shuffling).")
  #logger.info("testing individual structural elements for significance (nucleotide shuffling).")
  struc.finalize_structure()
  return(struc.finalStructure)

def write_final_alignment(alignment, structure, prefix):
  longestID = max([len(x.id) for x in alignment])
  with open(f"{outdir}/{prefix}_finalAlignment.stk",'w') as outputStream:
    outputStream.write("# STOCKHOLM 1.0\n")
    outputStream.write("#=GF AU  Kevin Lamkiewicz\n")
    outputStream.write("#=GF BM  VeGETA v. 0.3\n")
    outputStream.write(f"#=GF SQ  {len(alignment)}\n\n")
    
    for record in alignment:
    #for header, sequence in alignment.items():
      spacesToFill = longestID - len(record.id) + 5
      outputStream.write(f"{record.id}{' '*spacesToFill}{str(record.seq).replace('T','U')}\n")
    spacesToFill = longestID - len('#=GC SS_cons') + 5
    outputStream.write(f"#=GC SS_cons{' '*spacesToFill}{structure}\n//\n")

  #virusAligner.calculate_pw_distances()
  #virusAligner.get_tree_from_dist()
  #treePath = f"{os.path.dirname(outdir)}/test_tree.nwk"
  #Phylo.write(virusAligner.tree, treePath, 'newick')
  #virusAligner.refine_pairwise_instances(virusAligner.tree, None)
  

if __name__ == "__main__":
  logger = create_logger()
  (inputSequences, goi, outdir, alnOnly, clusterOnly, k, proc, tbpp, subcluster, seedSize, windowSize, stepSize, shannon, allowLP, shuffle, pvalue) = parse_arguments(docopt(__doc__))

  structureParameter = (logger, outdir, windowSize, stepSize, proc, allowLP, tbpp, shuffle, pvalue)

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
