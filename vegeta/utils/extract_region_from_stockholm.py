#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Author: Kevin Lamkiewicz
# Email: kevin.lamkiewicz@uni-jena.de


"""
Usage:
  extract_region_from_stockholm [--accession ID] <ALIGNMENT> <START> <STOP>

Option:
  --accession ID    Specifies whether the start and stop coordinates are based on genome coordinates or alignment columns. 
                    If none is provided, the alignment gets extracted at START and STOP column. If an ID is given
                    the alignment gets extracted at genomic coordinates. [Default: ]
"""


import sys
import os
from collections import Counter

from Bio import AlignIO
from docopt import docopt


def parse_arguments(d_args):
  """
  """

  alignment = AlignIO.read(d_args['<ALIGNMENT>'], 'stockholm')

  start = int(d_args['<START>'])
  stop = int(d_args['<STOP>'])
  accID = d_args['--accession']

  return alignment, start, stop, accID



def extract_fragment(alignment, start, stop, accID):
  """
  """

  if accID:
    for record in alignment:
      if str(record.id) == accID:
        gapsUntilStart = Counter(record.seq[:start])
        gapsUntilStop = Counter(record.seq[:stop])
        if not gapsUntilStart:
          startShift = 0
        else:
          startShift = gapsUntilStart['-']
        stopShift = gapsUntilStop['-'] 

        start = start+startShift
        stop = stop+stopShift 

  return alignment[:, start:stop]

if __name__ == '__main__':
  alignment, start, stop, accID = parse_arguments(docopt(__doc__))
  fragment = extract_fragment(alignment, start, stop, accID)
  output = f"fragment_{start}_{stop}.stk"
  
  if accID:
    output = f"{accID}_{output}"

  with open(output, 'w') as outputStream:
    AlignIO.write(fragment, outputStream, 'stockholm')
  exit(0)
