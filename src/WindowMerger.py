#!/usr/bin/env python3

import sys
import RNA
from collections import defaultdict
from collections import Counter
import operator

from Window import SequenceWindow


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
                gapless = str(record.seq).replace('-', '')
                try:
                    startIndex = self.sequences[record.id].upper().find(gapless)
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
            self.windowList.append((SequenceWindow(i, seqOrder, record_nts, colOverview)))


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
                if usedPositions[recordID] + 1 in genomicPositions:# and not majorityGap:
                    position = usedPositions[recordID]+1
                    usedPositions[recordID] = usedPositions[recordID]+1
                else:
                    position = -1

                self.mergedAlignment[col_idx][recordID] = position

    def minimize_gaps(self):
        """

        """

        pass

    def minimize_energy(self, e=2):
        """

        """

        pass
    pass
