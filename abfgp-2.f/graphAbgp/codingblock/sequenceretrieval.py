"""
Functions for CodingBlockGraphs (and LowSimilarityCodingBlockGraphs) for
retrieving the sequences of the Orfs of the CBGs
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"


class CodingBlockGraphSequenceRetievalFunctions:
    """
    Subclass with seqeuence retrieval functions for CodingBlockGraphs
    
    @attention: not functional on its own, required to be inherited from
    """
    def getorfproteinsequences(self,organism=None,node=None):
        """
        Get the orf sequence of each of the organisms of this graph

        @type  organism: *
        @param organism: Organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)

        @rtype:  dictionary (or string if organism is specified)
        @return: dictionary with organisms (keys) and orf protein sequence
                 strings (values), or only a orf protein sequence string if an
                 Organism identifier was specified
        """
        sequences = {}
        for org,orflist in self.get_orfs_of_graph().iteritems():
            if organism and org != organism: continue
            thenode = self.node_by_organism(org)
            if node and thenode != node: continue
            sequences[org]=orflist[0].protein_sequence
        if not organism:
            return sequences
        else:
            return sequences.values()[0]

    # end of function getorfproteinsequences


    def getomsrdnasequences(self,organism=None,node=None):
        """
        Get the OMSR DNA sequence of each of the Organism identifiers

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: Node identifier (or None)
    
        @rtype:  dictionary (or string if organism is specified)
        @return: dictionary with organisms (keys) and *COORD* protein sequence
                 strings (values), or only a *COORD* protein sequence string
                 if an organism identifier was specified
        """
        coords = self.overall_minimal_spanning_range()
        sequences = {}
        for org,orflist in self.get_orfs_of_graph().iteritems():
            if organism and org != organism: continue
            thenode = self.node_by_organism(org)
            if node and thenode != node: continue
            orf = orflist[0]
            absAaSta = min(coords[thenode])
            absAaEnd = max(coords[thenode]) + 1
            dnaseq = orf.getntseqslice(
                    abs_aa_startPY=absAaSta,
                    abs_aa_endPY=absAaEnd)
            sequences[org]=dnaseq
        if not organism:
            return sequences
        else:
            return sequences.values()[0]
        
    # end of function getomsrdnasequences


    def getomsrproteinsequences(self,**kwargs):
        """
        Get the OMSR protein sequence of each of the Organisms/Genes/Nodes

        @attention: see self._get_sequences_by_coords() for args and kwargs
        """
        coords = self.overall_minimal_spanning_range()
        return self._get_sequences_by_coords(coords,**kwargs)

    # end of function getomsrproteinsequences


    def getmaxsrproteinsequences(self,**kwargs):
        """
        Get the MAXSR protein sequences of each of the Organisms/Genes/Nodes

        @attention: see self._get_sequences_by_coords() for args and kwargs
        """
        coords = self.maximal_spanning_range()
        return self._get_sequences_by_coords(coords,**kwargs)

    # end of function getmaxsrproteinsequences


    def getminsrproteinsequences(self,**kwargs):
        """
        Get the MINSR protein sequences of each of the Organisms/Genes/Nodes

        @attention: see self._get_sequences_by_coords() for args and kwargs
        """
        coords = self.minimal_spanning_range()
        return self._get_sequences_by_coords(coords,**kwargs)

    # end of function getminsrproteinsequences


    def getomsr2orfendproteinsequences(self,**kwargs):
        """
        Get the OMSR2ORFEND AA sequences of each of the Organisms/Genes/Nodes

        @attention: see self._get_sequences_by_coords() for args and kwargs
        """
        coords = self.omsr2orfend()
        return self._get_sequences_by_coords(coords,**kwargs)

    # end of function getomsr2orfendproteinsequences


    def get_maxsr_proteinsequences_and_coords(self,**kwargs):
        """
        Get the MAXSR protein sequences and coordinate dict of each of the Nodes

        @attention: see self._get_sequences_by_coords() for args and kwargs
        """
        coords = self.maximal_spanning_range()
        return self._get_sequences_by_coords(coords,**kwargs), coords

    # end of function get_maxsr_proteinsequences_and_coords


    def get_omsr_proteinsequences_and_coords(self,**kwargs):
        """
        Get the OMSR protein sequences and coordinate dict of each of the Nodes

        @attention: see self._get_sequences_by_coords() for args and kwargs
        """
        coords = self.overall_minimal_spanning_range()
        return self._get_sequences_by_coords(coords,**kwargs), coords

    # end of function get_omsr_proteinsequences_and_coords


    def get_msr_proteinsequences_and_coords(self,**kwargs):
        """
        Get the MSR protein sequences and coordinate dict of each of the Nodes

        @attention: see self._get_sequences_by_coords() for args and kwargs
        """
        coords = self.minimal_spanning_range()
        return self._get_sequences_by_coords(coords,**kwargs), coords

    # end of function get_msr_proteinsequences_and_coords


    def get_left_sprdif_proteinsequences(self,organism=None,node=None,**kwargs):
        """ """
        coords = self.left_spanningrange_difference(**kwargs)
        if organism:
            node = self.node_by_organism(organism)
        if node and not coords.has_key(node):
            return {}
        elif node:
            coords = { node: coords[node] }
        else:
            pass
        return self._get_sequences_by_coords(coords,organism=organism,node=node)

    # end of function get_left_sprdif_proteinsequences


    def get_left_sprdif_proteinsequences_and_coords(self,organism=None,node=None,**kwargs):
        """ """
        coords = self.left_spanningrange_difference(**kwargs)
        if organism:
            node = self.node_by_organism(organism)
        if node and not coords.has_key(node):
            return {}
        elif node:
            coords = { node: coords[node] }
        else:
            pass
        return self._get_sequences_by_coords(coords,organism=organism,node=node), coords

    # end of function get_left_sprdif_proteinsequences_and_coords



    def _get_sequences_by_coords(self,coords,organism=None,node=None):
        """
        @type  organism: *
        @param organism: organism identifier (or None)

        @type  node: *
        @param node: node identifier (or None)

        @type  coords: dict
        @param coords: dict with Nodes as keys, ranges or (min,max) as values

        @rtype:  dictionary (or string if organism is specified)
        @return: dictionary with Node (keys) and orf protein sequence
                 substrings (values), or only a orf protein sequence string
                 if an Organism or Node identifier was specified
        """
        sequences = {}
        for org,orflist in self.get_orfs_of_graph().iteritems():
            if organism and org != organism: continue
            thenode = self.node_by_organism(org)
            if not coords.has_key(thenode): continue
            if node and node != thenode: continue
            absSta = min(coords[thenode]) - orflist[0].protein_startPY
            absEnd = max(coords[thenode]) + 1 - orflist[0].protein_startPY
            sequences[thenode]=orflist[0].protein_sequence[absSta:absEnd]
        if not organism:
            return sequences
        else:
            return sequences.values()[0]

    # end of function _get_sequences_by_coords

# end of class CodingBlockGraphSequenceRetievalFunctions