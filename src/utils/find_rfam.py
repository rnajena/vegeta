#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os

stkAln = sys.argv[1]

d_sequences = {}

with open(stkAln, 'r') as inputStream:
  for line in stkAln:
    if line.startswith('#=GC SS'):
      
