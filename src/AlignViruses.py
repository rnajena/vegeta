#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import copy
from multiprocessing import Pool
import multiprocessing as mp
from collections import Counter, defaultdict
from itertools import groupby, combinations
import glob

from scipy.spatial.distance import squareform
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import SeqIO, AlignIO, Phylo
from Bio.Seq import Seq
from Bio.Alphabet import generic_rna
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

import RNA

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

  def __init__(self, logger, inputFile, k, proc, windowSize, stepSize, outdir):
    """
    """

    self.logger = logger
    self.k = k
    self.proc = proc
    self.windowSize = windowSize
    self.stepSize = stepSize


    self.inputFile = inputFile
    self.sequences = self.read_sequences()
    self.pwDistances = {}
    self.tree = None
    self.nodes2leaves = {}

    #self.dirName = os.path.dirname(self.inputFile)
    self.dirName = outdir


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
    self.nodes2leaves = {x.name : [y.name for y in x.get_terminals()] for x in self.tree.get_nonterminals()}



  def locarna(self, q, lock):
    while True:
        file = q.get()
        if file is None:
            break
        #with lock:
        #print(f"Running mlocarna on {file}")
        cmd = f"mlocarna -q --keep-sequence-order --stockholm --noLP {file}"
        subprocess.run(cmd.split())
        # os.system(f"mlocarna -q --keep-sequence-order --stockholm --noLP {file}")

  def create_locarna_alignments(self, msa, alnLength):
    """
    doc string.
    """
    #self.stepSize, self.windowSize, self.proc, self.dirName
    # alnWindows = [(start1, end1), (start2,end2), ... , (start_n,end_n)]
    tmpdir = f"{self.dirName}/tmpSequences/"
    alnWindows = [(0, x) for x in range(self.stepSize, self.windowSize, self.stepSize)] + \
        [(x, min(x+self.windowSize, alnLength)) for x in range(0, alnLength, self.stepSize)]

    for counter, (start, stop) in enumerate(alnWindows):
        current_window = [(x.id, str(x.seq)[start:stop]) for x in msa]
        gapThreshold = (stop - start)  * 0.2
        outpath = f"{tmpdir}/local_window_{counter}.fa"
        with open(outpath, 'w') as output_stream:
            gaps = [Counter(x[1])['-'] for x in current_window]
            if any([x > gapThreshold for x in gaps]):
                #print(current_window)
                output_stream.write("\n".join([f">{i}\n{seq}" for i, seq in current_window]))    
            else:
                gapless_sequences = [(x[0], x[1].replace('-', '')) for x in current_window]
                output_stream.write("\n".join([f">{i}\n{seq}" for i, seq in gapless_sequences]))
                #print(gapless_sequences)
    #sys.exit(0)

    q = mp.Queue(maxsize=self.proc)
    lock = mp.Lock()
    pool = Pool(self.proc, initializer=self.locarna, initargs=(q, lock))

    for file in glob.glob(f'{tmpdir}/*fa'):
        q.put(file)

    for _ in range(self.proc):
        q.put(None)  # tell workers, everything is done

    pool.close()
    pool.join()

    for file in glob.glob(f'{tmpdir}/*out/results/result.stk'):
        index = int(file.split('window_')[1].split('.')[0])
        os.system(f"cp {file} {tmpdir}/final_local_alignment_{index}.stk")
    os.system(f"rm -r {tmpdir}/*fa {tmpdir}/*out/")
    
  def format_window(self, msa):
    for record in msa:
        record.seq = str(record.seq).upper().replace('T', 'U')

  def read_locarna_alignments(self):
    """ 
    """
    tmpdir = f"{self.dirName}/tmpSequences/"
    locarna_windows = {}
    for file in glob.glob(f'{tmpdir}/final*.stk'):
        try:
            index = int(file.split('alignment_')[1].split('.')[0])
        except IndexError:
            print(file)
            sys.exit(1)
        locarna_windows[index] = AlignIO.read(
            file, 'stockholm', alphabet=generic_rna)
        self.format_window(locarna_windows[index])
    return locarna_windows

  def alifold(self, q, lock, structures):
    while True:
        # number, start, stop, fragment = q.get()
        item = q.get()
        if item is None:
            break
        
        number, start, stop, fragment = item
        
        fc = RNA.fold_compound(fragment)
        fc.pf()
        bpp = fc.bpp()
        structure, mfe = fc.mfe()
        basePairings = {}
        idx = []
        
        # Structure Hash
        for index, column in enumerate(structure):
            if column == '.':
                # self.basePairings[index+offset].append(-1)
                basePairings[index+start] = -1
                continue

            if column == '(':
                idx.append(index)

            if column == ')':
                index2 = idx.pop()
                # self.basePairings[index +
                #                   offset].append(index2+offset)
                basePairings[index + start] = index2+start
                # self.basePairings[index2 +
                #                   offset].append(index+offset)
                basePairings[index2 +
                                start] = index+start
        with lock:
            structures.append(StructureWindow(number, basePairings, bpp))

  def fold_msa_windowed(self, aln):
    """
    """
    RNA.cvar.ribo = 1
    RNA.cvar.noLP = 1

    
    alnLength = len(aln[list(aln.keys())[0]])
    windows = [(0, x) for x in range(self.stepSize, self.windowSize, self.stepSize)] + \
        [(x, min(x+self.windowSize, alnLength))
            for x in range(0, alnLength, self.stepSize)]

    #structures = []

    q = mp.Queue(maxsize=self.proc)
    structures = mp.Manager().list()
    lock = mp.Lock()
    pool = mp.Pool(self.proc, initializer=self.alifold, initargs=(q, lock, structures))

    for index, (start, stop) in enumerate(windows):
        alnFragment = [index, start, stop] + [[aln[x][start:stop] for x in aln]]
        if any([x for x in alnFragment[-1]]):    
            q.put(alnFragment)

    for _ in range(self.proc):
        q.put(None)
    pool.close()
    pool.join()
    #print(structures)
    
    return [x for x in structures]


  def derive_final_structure(self, localStructures, alnLength):
    mergedStructure = ['.' for _ in range(alnLength)]
    # mergedStructure = []
    currentStem = []
    usedNTs = set()
    for i in range(alnLength):
        if i in usedNTs:
            continue
        windows = [window for window in localStructures if i in window.get_range()]
        basePairs = [x.basePairs[i] for x in windows if x.basePairs[i] not in usedNTs]

        position = sorted(Counter(basePairs).most_common(), key=lambda x: (-x[1], -x[0]))
        
        mostOccuring = position[0][0]
        count = position[0][1]
        if count == 1:
            mostOccuring = -1
      
        if currentStem:
            if mostOccuring > currentStem[-1][1]:
                mostOccuring = -1
        
        if mostOccuring == -1:
            
            continue
        try:
            if i < mostOccuring:
                windowsInteractionPartner = [window for window in localStructures if mostOccuring in window.get_range()]
                basePairsInteraction = [x.basePairs[mostOccuring] for x in windowsInteractionPartner if x.basePairs[mostOccuring]]
                if not basePairsInteraction: #ToDo: Dafuq?
                    continue
                #print(basePairsInteraction)
                mostCommonPartner = sorted(Counter(basePairsInteraction).most_common(), key=lambda x: (-x[1], -x[0]))
                if mostCommonPartner[0][0] == i:
                #mergedStructure.append('(')
                    currentStem.append((i,mostOccuring))
                    mergedStructure[i] = '('
                    mergedStructure[mostOccuring] = ')'
                else:
                    continue
                #usedNTs.add(mostOccuring)
            elif (mostOccuring, i) in currentStem:
                currentStem.pop()
            else:
                continue
                #mergedStructure.append(')')
                #mergedStructure[i] = ')'
                #mergedStructure[mostOccuring] = '('
                #usedNTs.add(mostOccuring)
            #usedNTs.add(i)
        except IndexError:
            print("IndexError in merging!")
            print(i, mostOccuring, basePairsInteraction)
            print(basePairs)
            #for x in windowsInteractionPartner:
            #    print(x.basePairs)
            #print(windowsInteractionPartner)
            exit()
    
    # iterativeStep = ['.' for _ in range(alnLength)]

    done = False
    while not done:
        iterativeStep = self.discardLonelyPairs(mergedStructure)
        if iterativeStep == mergedStructure:
            done = True
        mergedStructure = iterativeStep
        

    return mergedStructure

  def discardLonelyPairs(self, structure):
    lps = []
    basePairs = []
    for i in range(1, len(structure)-1):
        if (structure[i] == '(' or structure[i] == ')') and structure[i-1] == '.' and structure[i+1] == '.':
            lps.append(i)

        if structure[i] == '(':
            basePairs.append(i)
        
        if structure[i] == ')':
            leftPartner = basePairs.pop()
            if leftPartner in lps or i in lps:
                structure[leftPartner] = '.'
                structure[i] = '.'
    return structure

  def refine_pairwise_instances(self, tree, msa):
    """
    """

    if not tree or tree.total_branch_length() == 0:
        return msa
    

    
    copiedTree = copy.deepcopy(tree)
    
    resolvableNodes = [clade for clade in copiedTree.get_nonterminals() if clade.count_terminals() == 2]
    finalMSA = None
    for node in resolvableNodes:
        print(f"Resolving {node.name} now.")
        sequences = {leaf : self.sequences[leaf] for leaf in self.nodes2leaves[node.name]}
        #print(sequences.keys())
        
        if all([type(x) == str for x in sequences.values()]):
            with open(f"{self.dirName}/{node.name}.fasta", 'w') as outputStream:
                for entry in sequences.items():
                    outputStream.write(f">{entry[0]}\n{entry[1]}\n")
        
            with open(f"{self.dirName}/{node.name}_mafft.aln", 'w') as outputStream:
                cmd = f"mafft --thread {self.proc} --quiet {self.dirName}/{node.name}.fasta"
                subprocess.run(cmd.split(), stdout=outputStream, check=True)
        else:
            with open(f"{self.dirName}/{node.name}_msaTable.txt", 'w') as msaMergeTable:
                with open(f"{self.dirName}/{node.name}.fasta", 'w') as outputStream:
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
                        
                        
                
            
            with open(f"{self.dirName}/{node.name}_mafft.aln", 'w') as outputStream:
                cmd = f"mafft --thread {self.proc} --quiet --merge {self.dirName}/{node.name}_msaTable.txt {self.dirName}/{node.name}.fasta"
                subprocess.run(cmd.split(), stdout=outputStream, check=True)

        currentMSA = AlignIO.read(f"{self.dirName}/{node.name}_mafft.aln", 'fasta')
        #self.sequences[node.name] = currentMSA

        #for child in node.get_terminals():
        #  copiedTree.collapse(child)
          
    #return self.refine_pairwise_instances(copiedTree, currentMSA)
        
        currentMSA = AlignIO.read(f"{self.dirName}/{node.name}_mafft.aln", 'fasta')
        alnLength = len(currentMSA[0])
        #print("Length after mafft " + str(alnLength))
        self.create_locarna_alignments(currentMSA, alnLength)
        
        locarna_windows = self.read_locarna_alignments()
        windowMerger = WindowMerger(locarna_windows, self.stepSize, self.windowSize, sequences)
        
        windowMerger.merge_windows2()
        mergedAlignment = windowMerger.mergedAlignment
        #structureAlignment = format_alignment(mergedAlignment, sequences)
        structureAlignment = mergedAlignment
        # print([len(v) for k,v in structureAlignment.items()])
        alnLength = len(list(structureAlignment.values())[0])
        #print("Length after merged locarna " + str(alnLength))
        
        #exit(0)
        print("Folding Stuff!")
        localStructures = self.fold_msa_windowed(structureAlignment)
        finalStructure = self.derive_final_structure(localStructures, alnLength)
        finalMSA = MultipleSeqAlignment([SeqRecord(Seq(sequence, generic_rna), id=recordID) for recordID, sequence in structureAlignment.items()])
        #print(finalMSA)
        #exit(0)
        os.system(f"rm -r {self.dirName}/tmpSequences/final* ")
        finalMSA.column_annotations['secondary_structure'] = ''.join(finalStructure)
        
        self.sequences[node.name] = finalMSA
        
        for child in node.get_terminals():
            copiedTree.collapse(child)
        
        outputPath = f"{self.dirName}/{node}_vegeta.fasta"
        AlignIO.write(finalMSA, outputPath, 'fasta')
        outputPath = f"{self.dirName}/{node}_vegeta.stk"
        AlignIO.write(finalMSA, outputPath, 'stockholm')
        #exit(0)
    #exit(0)
    return self.refine_pairwise_instances(copiedTree, finalMSA)

  def perform_mafft(self):
    """
    """

    dirName = os.path.dirname(self.inputFile)
    with open(f"{dirName}/initial_mafft.fasta", 'w') as outputStream:
        cmd = f"mafft --quiet --thread {self.proc} {self.inputFile}"
        print(cmd)
        subprocess.run(cmd.split(), stdout=outputStream , check=True)
    
################################################

class Window(object):
    """

    """

    # def __init__(self, number, basePairs, seqOrder, nuclColumns, columnOverview):
    def __init__(self, number):
        self.number = number
        # self.basePairs = basePairs
        # self.seqOrder = seqOrder
        # self.nuclColumns = nuclColumns
        # self.columnOverview = columnOverview
        
        #sorted_columns = sorted(list(self.basePairs.keys()))

        self.range = (0, 0)

    def get_range(self):
        """

        """
        
        return range(self.range[0], self.range[1]+1)

    # def update_info(self, basePairs, seqOrder, nuclColumns, columnOverview):
        # self.basePairs = basePairs
        # self.seqOrder = seqOrder
        # self.nuclColumns = nuclColumns
        # self.columnOverview = columnOverview

class SequenceWindow(Window):

    def __init__(self, number, seqOrder, nuclColumns, columnOverview):
        super().__init__(number)
        self.seqOrder = seqOrder
        self.nuclColumns = nuclColumns
        self.columnOverview = columnOverview

        sorted_columns = sorted(list(self.nuclColumns.keys()))
        self.range = (min(sorted_columns), max(sorted_columns))

class StructureWindow(Window):

    def __init__(self, number, basePairs, bpp):
        super().__init__(number)
        self.basePairs = basePairs
        self.bpp = bpp
        sorted_columns = sorted(list(self.basePairs.keys()))
        self.range = (min(sorted_columns), max(sorted_columns))

################################################

class WindowMerger(object):
    """
    ToDos for window merging:
    Easy:
    1) Majority Vote of columns. If two nucleotides are pairing in every window, it is a pretty strong hint.
    Medium - Hard:
    2) Energy shift: If a nucleotide has different bp partner in different windows, RNAsubopt might help out
        Trying to evaluate how much energy has to be put into the different structures. Minimize energy this way
    Hard:
    3) If nucleotides are different for the same column in different windows (lets make it easy and this column only
        has one bp partner). Then, the context of the sequence between the two bp columns and the compensatory score
        may decide which one to trust
    Very Hard:
    4) Nucleotides are different in a column that also has different bp partners in different columns.
        In that case, I will simply quit my PhD and become a carpenter!
    """

    def __init__(self, allWindows, stepSize, windowSize, sequences):
        """

        """

        self.mergedAlignment = defaultdict(self.__ddl)

        self.allWindows = allWindows
        self.stepSize = stepSize
        self.windowSize = windowSize
        self.sequences = sequences
        #print(self.sequences.keys())

        self.windowList = []

        #self.col2genomic = defaultdict(self.ddl)
        #self.basePairings = defaultdict(list)

        self.alnLength = 0
        self.convert_windows_to_hashes()
        
    def __ddl(self):
        return defaultdict(list)

    def __ddi(self):
        return defaultdict(int)

    def __convert_window_to_columnHash(self, window, struc, offset):
        """

        """
        windowHash = defaultdict(dict)
        if not window or not struc:
            return windowHash

        alnLength = window.get_alignment_length()
        seqOrder = [x.id for x in window]

        for i in range(alnLength):
            column = window[:, i]
            for seqID, nt in (zip(seqOrder, column)):
                windowHash[i+offset].update({seqID: nt})
            windowHash[i+offset].update({'struc': struc[i]})
            windowHash[i+offset].update({'column': column})
        return(windowHash)

    def __get_windows_with_col(self, col_idx):
        """

        """

        return [window for window in self.windowList if col_idx in window.get_range()]

    def __convert_order_info(self, colSeqOrder, colNucleotides):
        order = defaultdict(list)
        for x, y in zip(colSeqOrder, colNucleotides):
            for a, b in zip(x, y):
                order[a].append(b)
        return order

    def __get_windows_with_nucleotide(self, nt_idx, seqID):
        """

        """

        coveringWindows = []
        alnColumns = []


        for window in self.windowList:
            seqIDX = window.seqOrder.index(seqID)
           # print(window.nuclColumns)
            for colIDX, nucleotides in window.nuclColumns.items():
                if nt_idx == nucleotides[seqIDX]:
                    coveringWindows.append(window)
                    alnColumns.append(colIDX)
                    continue
            
        #return(coveringWindows)
        return(alnColumns)

    def convert_windows_to_hashes(self):

        gapInSequence = defaultdict(self.__ddi)
        gaps = 0
        for i in range(0, len(self.allWindows.keys())):
            basePairings = {}
            record_nts = defaultdict(list)
            currentWindow = self.allWindows[i]
            if i*self.stepSize >= self.windowSize:
                #offset = int(self.stepSize * (i+1-self.windowSize/self.stepSize)  - gaps)
                offset = (i+1)*self.stepSize - self.windowSize - gaps
            else:
                offset = 0
        
            secondary_structure = currentWindow.column_annotations['secondary_structure']
            idx = []

            # Structure Hash
            for index, column in enumerate(secondary_structure):
                if column == '.':
                    # self.basePairings[index+offset].append(-1)
                    basePairings[index+offset] = -1
                    continue

                if column == '(':
                    idx.append(index)

                if column == ')':
                    index2 = idx.pop()
                    # self.basePairings[index +
                    #                   offset].append(index2+offset)
                    basePairings[index +
                                    offset] = index2+offset
                    # self.basePairings[index2 +
                    #                   offset].append(index+offset)
                    basePairings[index2 +
                                    offset] = index+offset

            # Sequence Hash
            seqOrder = []
            for record in currentWindow:
                #record_nts = self.col2genomic[record.id]
                gapless = str(record.seq).replace('-', '').upper().replace('U','T')
                try:
                    startIndex = self.sequences[record.id].upper().find(gapless)
                    #print(startIndex)
                    #print(gapless)
                    #print(self.sequences[record.id].upper())
                    #print(startIndex)
                except TypeError:
                    print(record.id)
                    sys.exit(0)
                
                seqOrder.append(record.id)
                for col_idx, column in enumerate(record.seq):
                    if column == '-':
                        record_nts[col_idx+offset].append(-1)
                        gapInSequence[record.id][col_idx +
                                                    startIndex] = gapInSequence[record.id][col_idx+startIndex-1] + 1
                    else:
                        gapInSequence[record.id][col_idx +
                                                    startIndex] = gapInSequence[record.id][col_idx+startIndex-1]

                        record_nts[col_idx+offset].append(
                            col_idx - Counter(record.seq[:col_idx])['-'] + startIndex)
                    self.alnLength = col_idx+offset

            colOverview = self.__convert_window_to_columnHash(currentWindow, secondary_structure, offset)
            # self.windowList.append((SequenceWindow(i, basePairings, seqOrder, record_nts, colOverview)))
            #print(i, record_nts)
            #print(colOverview)
            #print(self.alnLength)
            #print(i, record_nts)
            self.windowList.append((SequenceWindow(i, seqOrder, record_nts, colOverview)))


    def merge_windows2(self):
        """

        """
        aln = {}
        for seqID, sequence in self.sequences.items():
            lastColumnFilled = -1
            alnRow = ['-'] * (self.alnLength+2)
            #print(len(alnRow))
            #print(seqID)
            usedNTs = []
            try:
                for idx, nt in enumerate(sequence):
                    windowsWithNT = self.__get_windows_with_nucleotide(idx, seqID)
                    #print(idx, windowsWithNT)
                    windowsWithNT = [x for x in windowsWithNT if x > lastColumnFilled and lastColumnFilled + 100 > x and x < len(alnRow)]
                    #print(windowsWithNT)

                    if not windowsWithNT:
                        alnCol = lastColumnFilled + 1
                        if alnCol >= len(alnRow):
                            raise ValueError(f"Alignment out of bounds\n{idx},{seqID}, {len(sequence)}:\n{self.__get_windows_with_nucleotide(idx, seqID)} -- {lastColumnFilled},{len(alnRow)}\n{windowsWithNT}")
                            print(seqID, idx, windowsWithNT, alnCol)
                            #print(alnRow)
                            exit(1)
                    else:
                        alnCol = Counter(windowsWithNT).most_common(1)[0][0]
                        if alnCol >= len(alnRow):
                            raise ValueError(f"Alignment out of bounds\n{idx},{seqID}, {len(sequence)}:\n{self.__get_windows_with_nucleotide(idx, seqID)} -- {lastColumnFilled},{len(alnRow)}\n{windowsWithNT}")
                            print(seqID, idx, windowsWithNT)
                            exit(1)

                    #print(idx, nt)
                    #print(self.__get_windows_with_nucleotide(idx, seqID))
                    #print(windowsWithNT, alnCol)    
                    #print()
                    try:
                        alnRow[alnCol] = nt
                        lastColumnFilled = alnCol
                        usedNTs.append(nt)
                    except IndexError:
                        print(alnCol, len(alnRow))
                    #for window in windowsWithNT:
                    #    print(window.nuclColumns)
            except ValueError as err:
                print(f"Value Error: {err}")
                exit(1)

            #print(seqID, len(usedNTs), len(set(usedNTs)), len(sequence))
            #print(set(usedNTs))
            assert(len(usedNTs) == len(sequence))
            assert(''.join(alnRow).replace('-','') == sequence)
            aln[seqID] = alnRow
        
        seqIDs = []
        alnRows = []
        for x, y in aln.items():
            seqIDs.append(x)
            alnRows.append(y)
        
        toRemove = []
        for idx in range(len(alnRows[0])):
            try:
                if all([row[idx] == '-' for row in alnRows]):
                    toRemove.append(idx)
            except IndexError:
                print(idx, [len(x) for x in alnRows])
        
        #print(toRemove)
        for index in sorted(toRemove, reverse=True):
            for row in alnRows:
                del row[index]
        
        aln = {k : ''.join(v) for k,v in zip(seqIDs, alnRows)}

        self.mergedAlignment = aln
        #return(aln)


    def merge_windows(self):
        """

        """
        
        usedPositions = {x: -1 for x in self.sequences.keys()}
        for col_idx in range(self.alnLength+1):
            windowInfo = self.__get_windows_with_col(col_idx)

            colNucleotides = [x.nuclColumns[col_idx] for x in windowInfo]
            colSeqOrder = [x.seqOrder for x in windowInfo]
            
            order = self.__convert_order_info(colSeqOrder, colNucleotides)
            for recordID, genomicPositions in order.items():
                countedPositions = Counter(genomicPositions).most_common()
                highestOccurence = countedPositions[0][1]

                allHighOccurence = [x[0] for x in countedPositions if x[1] == highestOccurence]

                majorityGap = -1 in allHighOccurence

                noGapsList = [x for x in genomicPositions if x != -1]
                if noGapsList:
                    lastChance = min(noGapsList) == usedPositions[recordID] + 1 and len(set(noGapsList)) > 1
                else:
                    lastChance = False
                #ntIsMostCommon = sorted(countedPositions, key = lambda x:(x[1], x[0]), reverse=True)
                #ntIsMostCommon = usedPositions[recordID] + 1 in allHighOccurence
                ntIsMostCommon = genomicPositions.count(usedPositions[recordID] + 1) > 1

                print()
                print(recordID)
                print(col_idx)
                print(usedPositions)
                print(genomicPositions)
                print(sorted(Counter(genomicPositions).most_common(), key = lambda x:(x[1], x[0]), reverse=True)[0][0])
                print(noGapsList)
                print(lastChance)

                if lastChance or (ntIsMostCommon and not majorityGap):
                    print("Found it")
                    position = usedPositions[recordID]+1
                    usedPositions[recordID] = usedPositions[recordID]+1
                else:
                    position = -1
                print()

                self.mergedAlignment[col_idx][recordID] = position
        print(len(self.mergedAlignment))

# class WindowMerger(object):
#     """
#     ToDos for window merging:
#     Easy:
#     1) Majority Vote of columns. If two nucleotides are pairing in every window, it is a pretty strong hint.
#     Medium - Hard:
#     2) Energy shift: If a nucleotide has different bp partner in different windows, RNAsubopt might help out
#         Trying to evaluate how much energy has to be put into the different structures. Minimize energy this way
#     Hard:
#     3) If nucleotides are different for the same column in different windows (lets make it easy and this column only
#         has one bp partner). Then, the context of the sequence between the two bp columns and the compensatory score
#         may decide which one to trust
#     Very Hard:
#     4) Nucleotides are different in a column that also has different bp partners in different columns.
#         In that case, I will simply quit my PhD and become a carpenter!
#     """

#     def __init__(self, allWindows, stepSize, windowSize, sequences):
#         """

#         """

#         self.mergedAlignment = defaultdict(self.__ddl)

#         self.allWindows = allWindows
#         self.stepSize = stepSize
#         self.windowSize = windowSize
#         self.sequences = sequences

#         self.windowList = []

#         #self.col2genomic = defaultdict(self.ddl)
#         #self.basePairings = defaultdict(list)

#         self.alnLength = 0
#         self.convert_windows_to_hashes()
        
#     def __ddl(self):
#         return defaultdict(list)

#     def __ddi(self):
#         return defaultdict(int)

#     def __convert_window_to_columnHash(self, window, struc, offset):
#         """

#         """
#         windowHash = defaultdict(dict)
#         if not window or not struc:
#             return windowHash

#         alnLength = window.get_alignment_length()
#         seqOrder = [x.id for x in window]

#         for i in range(alnLength):
#             column = window[:, i]
#             for seqID, nt in (zip(seqOrder, column)):
#                 windowHash[i+offset].update({seqID: nt})
#             windowHash[i+offset].update({'struc': struc[i]})
#             windowHash[i+offset].update({'column': column})
#         return(windowHash)

#     def __get_windows_with_col(self, col_idx):
#         """

#         """

#         return [window for window in self.windowList if col_idx in window.get_range()]

#     def __convert_order_info(self, colSeqOrder, colNucleotides):
#         order = defaultdict(list)
#         for x, y in zip(colSeqOrder, colNucleotides):
#             for a, b in zip(x, y):
#                 order[a].append(b)
#         return order

#     def convert_windows_to_hashes(self):

#         gapInSequence = defaultdict(self.__ddi)
#         gaps = 0
#         for i in range(0, len(self.allWindows.keys())):
#             basePairings = {}
#             record_nts = defaultdict(list)
#             currentWindow = self.allWindows[i]
#             if i*self.stepSize >= self.windowSize:
#                 offset = (i+1)*self.stepSize - self.windowSize - gaps
#             else:
#                 offset = 0
        
#             secondary_structure = currentWindow.column_annotations['secondary_structure']
#             idx = []

#             # Structure Hash
#             for index, column in enumerate(secondary_structure):
#                 if column == '.':
#                     basePairings[index+offset] = -1
#                     continue

#                 if column == '(':
#                     idx.append(index)

#                 if column == ')':
#                     index2 = idx.pop()

#                     basePairings[index +
#                                 offset] = index2+offset
#                     basePairings[index2 +
#                                 offset] = index+offset

#             # Sequence Hash
#             seqOrder = []
#             for record in currentWindow:
#                 gapless = str(record.seq).replace('-', '')
#                 try:
#                     startIndex = self.sequences[record.id].upper().find(gapless)

#                 except TypeError:
#                     print(record.id)
#                     sys.exit(0)
                
#                 seqOrder.append(record.id)
#                 for col_idx, column in enumerate(record.seq):
#                     if column == '-':
#                         record_nts[col_idx+offset].append(-1)
#                         gapInSequence[record.id][col_idx +
#                                                 startIndex] = gapInSequence[record.id][col_idx+startIndex-1] + 1
#                     else:
#                         gapInSequence[record.id][col_idx +
#                                                 startIndex] = gapInSequence[record.id][col_idx+startIndex-1]

#                         record_nts[col_idx+offset].append(
#                             col_idx - Counter(record.seq[:col_idx])['-'] + startIndex)
#                     self.alnLength = col_idx+offset

#             colOverview = self.__convert_window_to_columnHash(currentWindow, secondary_structure, offset)
#             self.windowList.append((SequenceWindow(i, seqOrder, record_nts, colOverview)))


#     def merge_windows(self):
#         """

#         """
        
#         usedPositions = {x: -1 for x in self.sequences.keys()}
#         for col_idx in range(self.alnLength+1):
#             windowInfo = self.__get_windows_with_col(col_idx)

#             colNucleotides = [x.nuclColumns[col_idx] for x in windowInfo]
#             colSeqOrder = [x.seqOrder for x in windowInfo]
            
#             order = self.__convert_order_info(colSeqOrder, colNucleotides)
#             for recordID, genomicPositions in order.items():
#                 if usedPositions[recordID] + 1 in genomicPositions:# and not majorityGap:
#                     position = usedPositions[recordID]+1
#                     usedPositions[recordID] = usedPositions[recordID]+1
#                 else:
#                     position = -1

#                 self.mergedAlignment[col_idx][recordID] = position