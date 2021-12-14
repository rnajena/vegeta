#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de

"""
VeGETA -- Viral GEnome sTructure Alignments

VeGETA calculates a multiple sequence alignment of 
the viral genomes. This alignment is progressively refined
using 'MAFFT' and 'LocARNA'. In the end, we provide a local-structure guided
MSA of the viruses.

Python Dependencies:
  docopt
  BioPython
  colorlog
  numpy
  scipy

Other Dependencies:
  ViennaRNA package
  MAFFT
  LocARNA
  GLPK

Contact:
  kevin.lamkiewicz@uni-jena.de

Citation:
  Lamkiewicz, K., et al. (20xx), "Structure-guided multiple sequence alignments of complete viral genomes.", Journal, Issue, Volume.

Usage:
  vegeta.py [options] <inputSequences>

Options:
  -h, --help                              Show this help message and exits.
  -v, --verbose                           Get some extra information from VeGETA during calculation. [Default: False]
  --version                               Prints the version of VeGETA and exits.
  -o DIR, --output DIR                    Specifies the output directory of VeGETA. [Default: pwd]

  -p PROCESSES, --process PROCESSES       Specify the number of CPU cores that are used. [Default: 1]

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
  VeGETA v0.4 (alpha)
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

from vegeta import Aligner, StructCalculator 

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
    print("VeGETA version 0.4")
    exit(0)

  
  verbose = d_args['--verbose']
  if verbose:
    logger.setLevel(logging.INFO)

  inputSequences = d_args['<inputSequences>']
  if not os.path.isfile(inputSequences):
    logger.error("Couldn't find input sequences. Check your file")
    exit(1)
  
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

  allowLP = d_args['--allowLP']

  return (inputSequences, output, proc, tbpp, seedSize, windowSize, stepSize, shannon, allowLP, shuffle, pvalue)

def perform_alignment(seq):
  logger.info("Starting the alignment step of VeGETA.\n")


  files = [seq] 
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
    logger.info("Merging all fragments to a whole alignment")
    virusAligner.naive_consensus_majority()
    exit(0)

    structure = derive_structure(prefix)
    logger.info("Saving the final alignment in STOCKHOLM format")
    write_final_alignment(virusAligner.refinedAlignment, structure, prefix)
    #shutil.rmtree(f"{outdir}/tmpSequences")

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
  longestID = max([len(x.id) for x in alignment]+[len("#=GC SS_cons")])
  with open(f"{outdir}/{prefix}_finalAlignment.stk",'w') as outputStream:
    outputStream.write("# STOCKHOLM 1.0\n")
    outputStream.write("#=GF AU  Kevin Lamkiewicz\n")
    outputStream.write("#=GF BM  VeGETA v. 0.4\n")
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
  (inputSequences, outdir, proc, tbpp, seedSize, windowSize, stepSize, shannon, allowLP, shuffle, pvalue) = parse_arguments(docopt(__doc__))

  structureParameter = (logger, outdir, windowSize, stepSize, proc, allowLP, tbpp, shuffle, pvalue)
  perform_alignment(inputSequences)
