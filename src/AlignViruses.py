#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import subprocess
from multiprocessing import Pool

class Aligner(object):
  """
  """

  def __init__(self, logger, inputFile, k, proc):
    """
    """

    self.logger = logger
    self.k = k
    self.proc = proc
    
    self.inputFile = inputFile
    self.sequences = self.read_sequences()
    #print(self.sequences.keys())

  def read_sequences(self):
    """
    """

    fastaContent = {}
    
    with open(self.inputFile, 'r') as inputStream:
      header = ''
      seq = ''

      for line in inputStream:
        if line.startswith(">"):
          if header:
              fastaContent[header] = seq
          header = line.rstrip("\n").replace(':','_').replace(' ','_').lstrip(">")
          seq = ''
        else:
          seq += line.rstrip("\n").upper().replace('U','T')

      fastaContent[header] = seq
    return fastaContent
  
  def perform_mafft(self):
    
    dirName = os.path.dirname(self.inputFile)
    with open(f"{dirName}/initial_mafft.fasta", 'w') as outputStream:
      cmd = f"mafft --quiet --thread {self.proc} {self.inputFile}"
      subprocess.run(cmd.split(), stdout=outputStream , check=True)
    
