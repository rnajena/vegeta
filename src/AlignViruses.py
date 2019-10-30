#!/usr/bin/env python3
# -*- coding: utf-8 -*-

class Aligner(object):
  """
  """

  def __init__(self, logger, inputFile, k, proc, windowSize, stepSize, outdir):
    """
    """

    self.logger = logger


  def __getstate__(self):
    self_dict = self.__dict__.copy()
    del self_dict['pool']
    return self_dict

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