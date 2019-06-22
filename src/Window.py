#!/usr/bin/env python3
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