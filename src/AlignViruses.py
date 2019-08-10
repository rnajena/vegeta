#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
from multiprocessing import Pool
from collections import Counter
from itertools import groupby, combinations

from scipy.spatial.distance import squareform
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix

#######################################
# import sys
# import types
# #Difference between Python3 and 2
# import copyreg

# def _pickle_method(m):
#     class_self = m.im_class if m.im_self is None else m.im_self
#     return getattr, (class_self, m.im_func.func_name)

# copyreg.pickle(types.MethodType, _pickle_method)
#######################################

class Aligner(object):
  """
  """

  def __init__(self, logger, inputFile, k, proc):
    """
    """

    self.logger = logger
    self.k = k
    self.pool = Pool(processes=proc)
    self.proc = proc
    
    self.inputFile = inputFile
    self.sequences = self.read_sequences()
    self.pwDistances = {}
    self.tree = None


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
  
  def scan_kmer(self, sequence):
    """
    yet another doc string
    """
    return([sequence[x:x+self.k] for x in range(0,len(sequence)-self.k)])

  def get_pairwise_distances(self, pwComparison):
    """
    nice doc string.
    """
    kmers = self.scan_kmer(self.sequences[pwComparison[0]])
    secondKmers = self.scan_kmer(self.sequences[pwComparison[1]])
            
    count_one = Counter(kmers)
    count_two = Counter(secondKmers)

    sharedKmer = count_one.keys() & count_two.keys()
    uniqueOne = count_one.keys() - count_two.keys()
    uniqueTwo = count_two.keys() - count_one.keys()
            
    return (pwComparison,(sum(pow((count_one[x] - count_two[x]),2) for x in sharedKmer) + sum(pow(count_one[x],2) for x in uniqueOne) + sum(pow(count_two[x],2) for x in uniqueTwo)))

  def calculate_pw_distances(self): 
    for element in combinations(self.sequences,2):
      keys, pwDistance = self.get_pairwise_distances(element)
      self.pwDistances[keys] = pwDistance 
    

  def get_tree_from_dist(self):
    """
    nice doc string
    """
    distances = list(self.pwDistances.values())
    names = list(self.sequences.keys())
    dm = [[int(entry) for idx, entry in enumerate(sublist) if idx <= subID] for subID, sublist in enumerate(squareform(distances))]
    constuctor = DistanceTreeConstructor()
    self.tree = constuctor.nj(DistanceMatrix(names, dm)) # http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction.DistanceTreeConstructor-class.html
    self.tree.root_at_midpoint()
    


  def refine_pairwise_instances(self, tree, msa):
    """
    """

    if not tree or tree.total_branch_length() == 0:
        return msa
    

    sys.exit(0)
    copiedTree = copy.deepcopy(tree)
    
    resolvableNodes = [clade for clade in copiedTree.get_nonterminals() if clade.count_terminals() == 2]
    finalMSA = None
    for node in resolvableNodes:
        print(f"Resolving {node.name} now.")
        sequences = {leaf : singleSequences[leaf] for leaf in nodes2leaves[node.name]}
        #print(sequences.keys())
        
        if all([type(x) == str for x in sequences.values()]):
            with open(f"{tmpdir}/{node.name}.fasta", 'w') as outputStream:
                for entry in sequences.items():
                    outputStream.write(f">{entry[0]}\n{entry[1]}\n")
        
            with open(f"{tmpdir}/{node.name}_mafft.aln", 'w') as outputStream:
                cmd = f"mafft --thread {nCore} --quiet {tmpdir}/{node.name}.fasta"
                subprocess.run(cmd.split(), stdout=outputStream, check=True)
        else:
            with open(f"{tmpdir}/{node.name}_msaTable.txt", 'w') as msaMergeTable:
                with open(f"{tmpdir}/{node.name}.fasta", 'w') as outputStream:
                    seqCounter = 0
                    for name, record in sequences.items():
                        if type(record) != str:
                            for entry in record:
                                seqCounter += 1
                                
                                msaMergeTable.write(f"{seqCounter} ")
                                outputStream.write(f">{entry.id}\n{entry.seq}\n")
                            msaMergeTable.write("\n")
                        else:
                            seqCounter += 1
                            msaMergeTable.write(f"{seqCounter}\n")
                            outputStream.write(f">{name}\n{record}\n")
                        
                        
                
            
            with open(f"{tmpdir}/{node.name}_mafft.aln", 'w') as outputStream:
                cmd = f"mafft --thread {nCore} --quiet --merge {tmpdir}/{node.name}_msaTable.txt {tmpdir}/{node.name}.fasta"
                subprocess.run(cmd.split(), stdout=outputStream, check=True)
            
            
        currentMSA = AlignIO.read(f"{tmpdir}/{node.name}_mafft.aln", 'fasta')
        alnLength = len(currentMSA[0])
        #print("Length after mafft " + str(alnLength))
        create_locarna_alignments(currentMSA, alnLength, stepSize, windowSize, nCore, tmpdir)
        
        locarna_windows = read_locarna_alignments(tmpdir)
        windowMerger = WindowMerger(locarna_windows, stepSize, windowSize, sequences)
        
        mergedAlignment = windowMerger.merge_windows2()
        #mergedAlignment = windowMerger.mergedAlignment
        #structureAlignment = format_alignment(mergedAlignment, sequences)
        structureAlignment = mergedAlignment
        print([len(v) for k,v in structureAlignment.items()])
        alnLength = len(list(structureAlignment.values())[0])
        #print("Length after merged locarna " + str(alnLength))
        #print(structureAlignment.values())
        #exit(0)

        localStructures = fold_msa_windowed(list(structureAlignment.values()), windowSize, stepSize, nCore)
        finalStructure = derive_final_structure(localStructures, alnLength)
        finalMSA = MultipleSeqAlignment([SeqRecord(Seq(sequence, generic_rna), id=recordID) for recordID, sequence in structureAlignment.items()])
        #print(finalMSA)
        #exit(0)
        os.system(f"rm -r {tmpdir}/final* {tmpdir}/local*")
        finalMSA.column_annotations['secondary_structure'] = ''.join(finalStructure)
        
        singleSequences[node.name] = finalMSA
        
        for child in node.get_terminals():
            copiedTree.collapse(child)
        
        outputPath = f"{tmpdir}/{node}_viralign.fasta"
        AlignIO.write(finalMSA, outputPath, 'fasta')
        outputPath = f"{tmpdir}/{node}_viralign.stk"
        AlignIO.write(finalMSA, outputPath, 'stockholm')
        #exit(0)
    return refine_pairwise_instances(copiedTree, finalMSA)
  def perform_mafft(self):
    """
    """

    dirName = os.path.dirname(self.inputFile)
    with open(f"{dirName}/initial_mafft.fasta", 'w') as outputStream:
      cmd = f"mafft --quiet --thread {self.proc} {self.inputFile}"
      print(cmd)
      subprocess.run(cmd.split(), stdout=outputStream , check=True)
    
