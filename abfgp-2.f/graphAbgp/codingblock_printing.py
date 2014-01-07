"""
Class with functions for printing and making multiple alignments
of CodingBlockGraphs used in Aligment Based Gene Predictions
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from sets import Set


class CodingBlockGraphPrintingFunctions:
    """
    """
    def printmultiplealignment(self,organism=None,_linesize=150,asfasta=False):
        """
        Print the multiple alignment of the FC-CBG, guided on a organism

        @attention: See function multiplealignment for documentation
        @attention: no arguments means default (CLUSTALW-like) output, not a multifasta representation
        @attention: when no organism is applied, strongest connected node is taken as organism argument
        """
        print self.multiplealignment(
                organism=organism,
                _linesize=_linesize,
                asfasta=asfasta
                )

    # end of function printmultiplealignment


    def multiplealignment(self,organism=None,_linesize=100,asfasta=True):
        """
        Print the multiple alignment of the FC-CBG, guided on a organism

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  _linesize: positive integer (default 100)
        @param _linesize: length of a single line of alignment

        @type  asfasta: Boolean
        @param asfasta: True returns an aligned multi-fasta, False default CLUSTALW-like output

        @rtype:  string
        @return: multiline string to be printed as multiplealignment

        @attention: no arguments means a multifasta representation, not the default (CLUSTALW-like) output
        @attention: when no organism is applied, strongest connected node is taken as organism argument
        """
        if not self.has_overall_minimal_spanning_range():
            print "NO overall_minimal_spanning_range present, no printing...."
            print self
            return False

        # if no organism specified, take strongest connected node
        if not organism:
            thenode  = self.strongest_connected_node()
            organism = self._organism_from_node( thenode )
        else:
            thenode  = self.node_by_organism(organism)

        # get omsr of this codingblock
        allomsr = self.overall_minimal_spanning_range()
        omsr = list( allomsr[thenode] )
        omsr.sort()
        aoa = []
        # get pacbporfs of this organism
        pacbporfs = self.get_pacbps_by_organism(organism)

        for pos in omsr:
            aoa.append( [] )
            is_first = True
            gaps = []
            gap_seen = False
            similarities = []
            for pacbporf in pacbporfs:
                algPos = pacbporf.alignmentposition_by_query_pos(qposPY=pos)
                algPosObj = pacbporf._positions[algPos]
                if len(pacbporf._positions) > algPos+1:
                    # get next position to check if there are gaps that follow
                    # on current position
                    nextPosObj = pacbporf._positions[algPos+1]
                    if nextPosObj.query == '-':
                        # a gap in the query sequence itself!
                        gap_seen = True
                        gaps.append( [ nextPosObj.sbjct ] )
                        gappos = algPos+2
                        while pacbporf._positions[gappos].query == '-':
                            gaps[-1].append( pacbporf._positions[gappos].sbjct ) 
                            gappos+=1
                    else:
                        # no gaps for this sbjct sequence
                        gaps.append( [] )
                else:
                   # last position of _positions list reached; no nextPosObj so no gaps as well
                   pass
                if is_first:
                    aoa[-1].append( algPosObj.query )
                aoa[-1].append( algPosObj.sbjct )
                if algPosObj.isa in ['similarity','identity']:
                    similarities.append(1)
                is_first=False

            if not asfasta:
                # set alignment overview symbol (*,+ or space)
                if len(Set(aoa[-1]))==1:
                    aoa[-1].append('*')
                elif len(similarities) == len(pacbporfs):
                    aoa[-1].append('+')
                else:
                    aoa[-1].append(' ')

            # create gaps!
            if gap_seen:
                gap_len = max([len(g) for g in gaps])
                # make all the gaps uniform in length
                # TODO: confirm if the AAs in the sbjct sequences
                # in the gap can be aligned 
                for g in gaps:
                    while len(g) < gap_len:
                        g.append('-')
                # insert gaps for the query
                gaps.insert(0,['-']*gap_len)
                for i in range(0,len(gaps[0])):
                    aoa.append( [] )
                    for g in gaps:
                        aoa[-1].append( g[i] )
                    # append a space (no alignment) to the alignment line
                    if not asfasta: aoa[-1].append(' ')

        # correct for leading & tailing gaps in the multiple alignment.
        # Why? -> example:
        #     TLTMFYFFSNSKIDTIPVNFPGGPGSPSNQNFPGQQAPVAEEGDTALTKDDAPAQKPA-------
        #     TLTMFYFFSNSKIDTIPVNFPGGPGSPSNQNFPGQQAPVAEEGDTALTKDDAPAQKPAGEGAGTA
        #     VLAVFYFVSHSSYDVRPDQYPGGSG--DIQDTPKVDTP---------KLDDAPKKDQT-------
        #     VIAVLYFISTPSSVVTPLKFQQSTGNSK-----GSTAGSTTGGTTEVAHDEKTPAKSN-------
        #      + + ** *       *  +    *                        *+              
        # This can be the case now. Such a multiple alignment can arise when
        # the OMSR is confined to an area that ends with a gap -> obvious!
        # This must be corrected in the OMSR!!!
        # TODO: FIX THIS in the OMSR!

        # make the multiple alignment lines uniform in length
        # EXAMPLE: this will crash in a HMM search.
        #     RRFPLLFATAMHFLPFTLLSLAA--SILS-
        #     RAFNPLAAAMRCLRPGPILGLASSIPLLS-
        #     RSFLHASAHAAMRLPFHLLSLVA--TIITN
        #     EHRHPIMRIFGSLGTAFTLGLAT---IAS. 
        #                       * *     + +.
        # In this example, the dots (.) are empty / no space
        # This function adds spaces where it is needed.
        maxlen = max([ len(item) for item in aoa ])
        for pos in range(0,len(aoa)):
            while len(aoa[pos]) < maxlen:
                aoa[pos].insert(-1,'-')

        # print multiple aligment
        if not asfasta:
            return pacbporfs[0]._print_aoa(aoa,_linesize=_linesize)
        else:
            # add empty lines for each organism in between the alignment
            for i in range(0,len(aoa)):
                for linepos in range(0,len(pacbporfs)*2+2,2):
                    aoa[i].insert(linepos,'')
            # add the main organism identifier
            aoa[0][0] = ">%s_orf_%s_aapos_%s_%s" % ( thenode[0], thenode[1], min(omsr), max(omsr)+1 )

            # add identifiers of the other organisms
            othernodes = self.get_ordered_nodes()
            othernodes.remove(thenode)
            for headerpos in range(2,len(pacbporfs)*2+2,2):
                thisnode = othernodes[ (headerpos-2)/2 ]
                thisomsr = list( allomsr[thisnode] )
                aoa[0][headerpos] = ">%s_orf_%s_aapos_%s_%s" % ( thisnode[0], thisnode[1], min(thisomsr), max(thisomsr)+1 )
            # and return the multiple alignment; ignore _linesize!!
            return pacbporfs[0]._print_aoa(aoa,_linesize=len(aoa))


    # end of function multiplealignment


    def printdetailedmultiplealignment(self,organism=None,_linesize=40):
        """
        Print the detailed multiple alignment of the FC-CBG, guided on a organism

        @type  organism: *
        @param organism: organism identifier (or None)

        @type  _linesize: positive integer (default 100)
        @param _linesize: length of a single line of alignment
        """
        if not self.has_overall_minimal_spanning_range():
            print "NO overall_minimal_spanning_range present, no printing...."
            return False
        # if no organism specified, take strongest connected node
        if not organism:
            organism = self._organism_from_node( self.strongest_connected_node() )
        # get omsr of this codingblock
        omsr = self.overall_minimal_spanning_range(organism=organism)
        omsr = list(omsr)
        omsr.sort()
        aoa = []
        # get pacbporfs of this organism
        pacbporfs = self.get_pacbps_by_organism(organism)

        for pos in omsr:
            aoa.append( [] )
            is_first = True
            gaps = []
            similarities = []
            tcodescores, cnt = [], 0
            for pacbporf in pacbporfs:
                algPos = pacbporf.alignmentposition_by_query_pos(qposPY=pos)
                algPosObj = pacbporf._positions[algPos]
                nextPosObj = pacbporf._positions[algPos+1]
                if nextPosObj.sbjct_pos > algPosObj.sbjct_pos+1:
                    gaps.append( nextPosObj.sbjct_pos - algPosObj.sbjct_pos -1 )
                if is_first:
                    aoa[-1].append( algPosObj.query_dna_seq )
                    aoa[-1].append( " %s " % algPosObj.query )
                aoa[-1].append( " %s " % algPosObj.sbjct )
                if algPosObj.isa in ['similarity','identity']:
                    similarities.append(1)
                is_first=False
                # do tcode score calculation
                if algPosObj.sbjct_dna_start:
                    tcodescores.append( pacbporf.orfS.aa_pos_tcode_score(algPosObj.sbjct_pos) )
                    cnt+=1

            # set alignment overview symbol (*,+ or space)
            if len(Set(aoa[-1][1:]))==1:
                aoa[-1].append('***')
            elif len(similarities) == len(pacbporfs):
                aoa[-1].append('+++')
            else:
                aoa[-1].append('   ')

            # now get the EMBOSS-TCODE data about this position
            if not (algPosObj.isa_gap and algPosObj.query_dna_start==0):
                score = pacbporfs[0].orfQ.tcode_score(start=pos,end=pos+1)
                tcodescores.append(score)
                cnt+=1
                score = int(round(score*10))
                if score < 10:  aoa[-1].append(' %s ' % score)
                else:           aoa[-1].append('%s ' % score)
            else:
                aoa[-1].append('   ')

            # now calculate total tscore on this position
            summedtcodescore = sum(tcodescores) / float(cnt)
            score = int(round(summedtcodescore*10))
            if score < 10:  aoa[-1].append(' %s ' % score)
            else:           aoa[-1].append('%s ' % score)
            aoa[-1].append(' %s ' % cnt)

            # create gaps!
            if gaps:
                for i in range(0,max(gaps)):
                    aoa.append( [] )
                    aoa[-1].append( ['   '] )
                    aoa[-1].append( ['-'] )
                    aoa[-1].extend( [' ? '*len(pacbporfs)] )
                    aoa[-1].extend( ['   ','   '] )


        # print multiple aligment
        print self
        print pacbporfs[0]._print_aoa(aoa,_linesize=_linesize)

    # end of function printdetailedmultiplealignment

# end of class CodingBlockGraphPrintingFunctions

