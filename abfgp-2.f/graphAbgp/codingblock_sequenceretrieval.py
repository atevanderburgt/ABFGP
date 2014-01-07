"""
Functions used in CodingBlockGraphs and LowSimilarityCodingBlockGraphs for
retrieving the sequences of the orfs of the CBGs
"""

class CodingBlockGraphSequenceRetievalFunctions:

    def getorfproteinsequences(self,organism=None):
        """
        Get the Orf proteins sequences of each of the Organism identifiers

        @type  organism: *
        @param organism: organism identifier (or None)

        @rtype:  dictionary (or string if organism is specified)
        @return: dictionary with organisms (keys) and orf protein sequence
                 strings (values), or only a orf protein sequence string
                 if an organism identifier was specified
        """
        sequences = {}
        for org,orflist in self.get_orfs_of_graph().iteritems():
            if organism and org != organism: continue
            sequences[org]=orflist[0].protein_sequence
        if not organism:
            return sequences
        else:
            return sequences.values()[0]

    # end of function getorfproteinsequences


    def getomsrdnasequences(self,organism=None):
        """
        Get the OMSR DNA sequence of each of the Organism identifiers
        @type  organism: *
        @param organism: organism identifier (or None)
    
        @rtype:  dictionary (or string if organism is specified)
        @return: dictionary with organisms (keys) and *COORD* protein sequence
                 strings (values), or only a *COORD* protein sequence string
                 if an organism identifier was specified
        """
        coords = self.overall_minimal_spanning_range()
        sequences = {}
        for org,orflist in self.get_orfs_of_graph().iteritems():
            if organism and org != organism: continue
            node = self.node_by_organism(org)
            orf = orflist[0]
            absAaSta = min(coords[node])
            absAaEnd = max(coords[node]) + 1
            dnaseq = orf.getntseqslice(
                    abs_aa_startPY=absAaSta,
                    abs_aa_endPY=absAaEnd)
            sequences[org]=dnaseq
        if not organism:
            return sequences
        else:
            return sequences.values()[0]
        
    # end of function getomsrdnasequences


    def getomsrproteinsequences(self,organism=None):
        """
        Get the OMSR protein sequence of each of the Organism identifiers

        @attention: see _seqs_by_coords for argument documentation
        """
        coords = self.overall_minimal_spanning_range()
        return _seqs_by_coords(self,coords,organism=organism)

    # end of function getomsrproteinsequences


    def getmaxsrproteinsequences(self,organism=None):
        """
        Get the MAXSR protein sequences of each of the Organism identifiers

        @attention: see _seqs_by_coords for argument documentation
        """
        coords = self.maximal_spanning_range()
        return _seqs_by_coords(self,coords,organism=organism)

    # end of function getmaxsrproteinsequences


    def getomsr2orfendproteinsequences(self,organism=None):
        """
        Get the OMSR2ORFEND protein sequences of each of the Organism identifiers

        @attention: see _seqs_by_coords for argument documentation
        """
        coords = self.omsr2orfend()
        return _seqs_by_coords(self,coords,organism=organism)

    # end of function getomsr2orfendproteinsequences


# end of class CodingBlockGraphSequenceRetievalFunctions


def _seqs_by_coords(cbg,coords,organism=None):
    """
    Helper function for get***proteinsequences from Orfs in CBGs

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph

    @type  coords: dict
    @param coords: coordinate Set dict (omsr,maxsr,omsr2orfend,etc)

    @type  organism: *
    @param organism: organism identifier (or None)

    @rtype:  dictionary (or string if organism is specified)
    @return: dictionary with organisms (keys) and *COORD* protein sequence
             strings (values), or only a *COORD* protein sequence string
             if an organism identifier was specified
    """
    sequences = {}
    for org,orflist in cbg.get_orfs_of_graph().iteritems():
        if organism and org != organism: continue
        node = cbg.node_by_organism(org)
        absSta = min(coords[node]) - orflist[0].protein_startPY
        absEnd = max(coords[node]) + 1 - orflist[0].protein_startPY
        sequences[org]=orflist[0].protein_sequence[absSta:absEnd]
    if not organism:
        return sequences
    else:
        return sequences.values()[0]

# end of function _seqs_by_coords
