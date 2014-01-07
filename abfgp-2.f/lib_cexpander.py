"""
Cexpander functions for Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# import cexpander items
from parsers.cexpander import *

# ABGP imports
from graphAbgp.exceptions import (
    NoOverallMinimalSpanningRange,
    OrganismNotPresentInGraph
    )

# CodingBlockGraph Imports
from graphAbgp.codingblock_operations import (
    _update_cbg_with_pacbporf_replacements
    )

# Other Imports
from lib_fasta import writeMultiFasta
from lib_stopwatch import StopWatch
import pacb

# Python Imports
from sets import Set
import re

# Global Variable Import
from settings.codingblockgraph import (
    CBG_CEXPANDER_OMSRBORDERGAPS_MAX_BITSCORERATIO_THRESHOLD,
    CBG_CEXPANDER_OMSRBORDERGAPS_NONUNIFORM_MAX_AA_LENGTH,
    CBG_CEXPANDER_OMSRBORDERGAPS_NONUNIFORM_AA_OFFSET,
    CBG_CEXPANDER_OMSRBORDERGAPS_GAP_SIZE,
    )


def cexpanderanalyses(cbg,min_cols=0,projected_on=":::",
    output='binary',cbgregion='omsr',verbose=False):
    """
    Run cexpander and get the CexpanderOutput object of this CBG
    
    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph

    @type  min_cols: DEPRECATED integer
    @param min_cols: DEPRECATED default 0, only change when expert user!

    @type  projected_on: DEPRECATED string
    @param projceted_on: DEPRECATED default ':::', only change when expert user!

    @type  output: string
    @param output: one of 'binary', 'float' (default 'binary')

    @type  cbgregion: string
    @param cbgregion: one of 'omsr', 'maxsr', 'omsr2orfend' (default 'omsr')

    @type  verbose: Boolean
    @param verbose: print intermediate info to STDOUT for debugging purposes

    @rtype  cxpdrOutput: CexpanderOutput
    @return cxpdrOutput: CexpanderOutput object of this CBG

    """
    # create (~unique) basefname
    nodestringlist = []
    for node in cbg.get_ordered_nodes():
        nodestringlist.append( "%s%s" % (node[0],node[1]) )
    fname_fasta = "cbg_" + "".join(nodestringlist) + ".fa"

    if cbg.node_count() > 2:
        # write multi fasta of OMSR sequences
        if cbgregion == 'maxsr':
            fastaseqs= cbg.getmaxsrproteinsequences()
        elif cbgregion == 'omsr2orfend':
            fastaseqs= cbg.getomsr2orfendproteinsequences()
            # append dummy W amino acid that mimics the
            # alignment of STOP codons
            for k,v in fastaseqs.iteritems():
                fastaseqs[k] = v+"W"
        else:
            # omsr or Non-exstsing keyword...
            fastaseqs= cbg.getomsrproteinsequences()
        writeMultiFasta( fastaseqs , fname_fasta )
        if projected_on == ":::":
            pass
        elif not projected_on:
            strongestN   = cbg.strongest_connected_node()
            strongestO   = cbg.organism_by_node(strongestN)
            projected_on = strongestO
        elif projected_on in cbg.organism_set():
            pass
        else:
            raise OrganismNotPresentInGraph
    
        # get cxpdrOutput object; file-cleanup is taken care for
        cxpdrOutput = runcexpander(fname_fasta,
                cbalignp_commandline = " -y", output=output)

        # correct hard-added 'W' residue in cexpander OMSR2ORFEND 
        if cbgregion == 'omsr2orfend':
            IS_FIRST = True
            for trf in cxpdrOutput._transferblocks: 
                trf.positions-=1
                if trf.binarystring[-1] == "1":
                    trf.score-=1
                trf.binarystring = trf.binarystring[0:-1]
                trf.ratio        = trf._binarystring2matchratio(trf.binarystring)
                if IS_FIRST:
                    cxpdrOutput.set_transferblock(trf.header)
                    IS_FIRST = False
            for k,seq in cxpdrOutput.sequences.iteritems():
                cxpdrOutput.sequences[k] = seq[0:-1]
        # EOF correct hard-added 'W' residue in cexpander OMSR2ORFEND

    else:
        # weird case of CBG with only 2 nodes;
        # can happen when ``weakest organism`` is removed
        # from the GSG based on GTG analyses
        # Fake cexpander output here
        pacbporf = cbg.pacbps.values()[0]
        bstring  = "1" *( pacbporf._original_alignment_pos_end -\
                          pacbporf._original_alignment_pos_start )
        cxpdrOutput = CexpanderOutput()
        cxpdrOutput.binarystring = bstring
        cxpdrOutput.header = cbg.organism_by_node(
                cbg.get_ordered_nodes()[0])
        cxpdrOutput.positions = len(bstring)
        cxpdrOutput.score = len(bstring)

    # (2) return the output object
    return cxpdrOutput

# end of function cexpanderanalyses


def cexpanderanalyses_omsr(cbg,**kwargs):
    """
    Default cexpanderanalyses using the OMSR of a CBG

    @attention: see cexpanderanalyses() for argument documentation
    """
    return cexpanderanalyses(cbg,cbgregion='omsr',**kwargs)

# end of function cexpanderanalyses_omsr


def cexpanderanalyses_maxsr(cbg,**kwargs):
    """
    Default cexpanderanalyses using the OMSR of a CBG

    @attention: see cexpanderanalyses() for argument documentation
    """
    return cexpanderanalyses(cbg,cbgregion='maxsr',**kwargs)

# end of function cexpanderanalyses_maxsr


def cexpanderanalyses_omsr2orfend(cbg,**kwargs):
    """
    Default cexpanderanalyses using the OMSR of a CBG

    @attention: see cexpanderanalyses() for argument documentation
    """
    return cexpanderanalyses(cbg,cbgregion='omsr2orfend',**kwargs)

# end of function cexpanderanalyses_omsr2orfend


def get_cexpander_string_of_cbg(cbg,verbose=False):
    """
    Run cexpander and get the cexpander binary string for this CBG

    @attention: recommended not to use, use more detailed cexpanderanalyses()

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph

    @type  verbose: Boolean
    @param verbose: print intermediate info to STDOUT for debugging purposes
    """
    # (0) create (~unique) basefname
    basefname = "cbg_"+ "".join([ "%s%s" % (node[0],node[1]) for node in cbg.get_ordered_nodes() ])
    fname_fasta     = basefname+".fa"
    fname_allvsall  = basefname+".allvsall"
    fname_aligned   = basefname+".aligned"
    fname_cexpander = basefname+".cexpander"

    # (1) write multi fasta of OMSR sequences
    omsrseqs= cbg.getomsrproteinsequences()
    writeMultiFasta( omsrseqs , fname_fasta )
    strongestN  = cbg.strongest_connected_node()
    strongestO  = cbg.organism_by_node(strongestN)

    # create complete .fa -> cexpanderstring command
    command = """
        python %s %s %s;
        %s -i %s %s > %s;
        %s %s ">%s" > %s;
        cat %s | grep "\$start_values" -A 1000 | grep "\$end_values" -m 1 -B 1000 | sed 's/^.*_values$//' | sed '/^$/d' | tr -d "\\t" | tr -d "\\n"
        """ % (
        EXECUTABLE_CEXPANDER_ALLVSALL, fname_fasta, fname_allvsall,
        EXECUTABLE_CEXPANDER_CBALIGNP, fname_allvsall, "-y", fname_aligned,
        EXECUTABLE_CEXPANDER_CEXPANDER, fname_aligned, strongestO, fname_cexpander,
        fname_cexpander )

    ci,co,ce = osPopen3(command)
    ci.close()
    # output of EXECUTABLE_CEXPANDER_ALLVSALL is cast to STDOUT as well!
    cexpanderstring = co.readlines()[-1]
    co.close()
    error = ce.read()
    ce.close()

    # (6) cleanup files
    osSystem("rm -f %s.*" % basefname )

    ############################################################
    if verbose:
        linesize=100
        print cbg
        for offset in range(0,len(cexpanderstring),linesize):
            print omsrseqs[strongestO][offset:offset+linesize]
            print cexpanderstring[offset:offset+linesize]
    ############################################################

    # (7) return the cexpander string
    return ( cexpanderstring, strongestO )

# end of function get_cexpander_string_of_cbg



def cexpander2multiplealignment(cxpdr,verbose=False):
    """
    This function and its application are still under development. In future
    version, this cexpander obtained data will replace the (deprecated) PAOC
    and PASC VISTA-like tracks which were far to computationally expensive to
    obtain.
    """
    ########################################################################
    if verbose:
        stw = StopWatch(name="cxpdr2multiplealignment")
        stw.start()
    ########################################################################


    # for each of the _transferblocks (1 for each organism/gene), the
    # binarystring **should** contain an identical number of 1's
    # in freak-accident cases (1 in hundreds of thousand of cases),
    # it is observed that this not the case. Catch this exception here
    # before it hard-crashes with a raise somewhere later in this function
    if len(Set([ trf.binarystring.count("1") for trf in cxpdr._transferblocks ])) > 1:
        print "WARNING: unequal Cexpander.transferblocks.binarystring 1's count:",
        print Set([ trf.binarystring.count("1") for trf in cxpdr._transferblocks ]) 
        return False

    # split the cexpander binarystrings on character changes 0->1 and 1->0
    substrings = {}
    orgs   = [ trf.header for trf in cxpdr._transferblocks ]
    for ipos in range(0,len(orgs)):
        org = orgs[ipos]
        trf = cxpdr._transferblocks[ipos]
        substrings[org] = [ x.group() for x in re.finditer("(1+|0+)",trf.binarystring) ]

    # maximum number of blocks in the cexpander output
    # WARNING TODO THIS IS STILL NOT 100% SAFE!!
    try:
        maxblocks = max(Set([ len(substrings[org]) for org in substrings.keys() ]))
    except:
        print "ERROR in cexpander2multiplealignment"
        print substrings.keys()
        print "inputseqs:", len(cxpdr.sequences)
        for k,v in substrings.iteritems():
            print k, len(v)
            print v
        # now raise the error...
        maxblocks = max(Set([ len(substrings[org]) for org in substrings.keys() ]))

    curblock = 0
    ########################################################################
    if verbose:
        print "maxblocks:", maxblocks,
        print [ len(substrings[org]) for org in substrings.keys() ]
        if len(Set([ len(substrings[org]) for org in substrings.keys() ])) > 1:
            for ipos in range(0,len(orgs)):
                org = orgs[ipos]
                print org,
                print [ Set(substrings[org][block]) for block in range(0,len(substrings[org])) ] 
                trf = cxpdr._transferblocks[ipos]
                print trf.binarystring, len(trf.binarystring),
                print trf.binarystring.count("1"), trf.binarystring.count("0")
    ########################################################################
    while curblock < maxblocks:
        try:
            # create curblocktypeset
            curblocktypeset = Set("".join([ substrings[org][curblock] for org in substrings.keys()]))
        except IndexError:
            # substrings[org][curblock](s) IndexError
            # can happen on EOF blocks if some have zeros, others have nothing
            # append empty block; this will be dealth with in the 
            # curblocktypeset Set("0")
            for org in substrings.keys():
                if len(substrings[org]) == curblock:
                    substrings[org].append("")
            # recreate curblocktypeset in 2th instance
            curblocktypeset = Set("".join([ substrings[org][curblock] for org in substrings.keys()]))

        ########################################################################
        if verbose:
            print "curiter::", curblock, maxblocks,
            print [ len(substrings[org][curblock]) for org in substrings.keys()]
        ########################################################################

        if curblocktypeset == Set("1"):
            # block of just ones; settle this block by limiting on minimal length
            # of all organisms of 111-string.
            curblocklengths = Set([ len(substrings[org][curblock]) for org in substrings.keys() ])
            if len(curblocklengths) == 1:
                pass # all normal...
            else:
                minlength = min(curblocklengths)
                for org in substrings.keys():
                    if len(substrings[org][curblock]) > minlength:
                        blocklen = len(substrings[org][curblock])
                        substrings[org][curblock] = substrings[org][curblock][0:minlength]
                        substrings[org].insert(curblock+1,"1"*(blocklen-minlength))
                        substrings[org].insert(curblock+1,"")
                # increase maxblocks counter
                maxblocks = max(Set([ len(substrings[org]) for org in substrings.keys() ]))
                ####################################################################
                if verbose:
                    print "TRBLOCKS CHANGED!, curblock, maxblocks:", curblock, maxblocks,
                    print [ len(substrings[org]) for org in substrings.keys() ]
                    for ipos in range(0,len(orgs)):
                        org = orgs[ipos]
                        print org,
                        print [ Set(substrings[org][block]) for block in range(0,len(substrings[org])) ]
                ####################################################################
        elif curblocktypeset == Set("0"):
            # check lengths of the blocks
            lengths = [ len(substrings[org][curblock]) for org in substrings.keys()]
            for org in substrings.keys():
                if len(substrings[org][curblock]) != max(lengths):
                    substrings[org][curblock] += "."*(max(lengths)-len(substrings[org][curblock]))
        elif curblocktypeset == Set(["0","1"]):
            # situation where frontal or intermediate zeros complicate the multiplealignment
            for org in substrings.keys():
                if Set(substrings[org][curblock]) == Set(['1']):
                    substrings[org].insert(curblock,"")
            # next, do as if cublocktypeset == Set("0") (which it is now!
            # check lengths of the blocks
            lengths = [ len(substrings[org][curblock]) for org in substrings.keys()]
            for org in substrings.keys():
                if len(substrings[org][curblock]) != max(lengths):
                    substrings[org][curblock] += "."*(max(lengths)-len(substrings[org][curblock]))
        else:
            print "MIXED!!", curblocktypeset, "curblock:", curblock, "maxblocks:", maxblocks
            print "ERROR WILL LIKELY OCCUR QUICKLY AFTER HERE..."
            pass
            import sys
            sys.exit()

        # increase the blocks counter
        curblock+=1
    ########################################################################
    if verbose:
        for org in substrings.keys():
            # print the sequence itself
            for block in range(0,maxblocks):
                offset   = sum([ substrings[org][i].count("1")+substrings[org][i].count("0") for i in range(0,block) ])
                blocklen = len(substrings[org][block])
                if Set(substrings[org][block]) == Set("1"):
                    print cxpdr.sequences[org][offset:offset+blocklen].upper(),
                elif Set(substrings[org][block]) == Set("0"):
                    print cxpdr.sequences[org][offset:offset+blocklen].lower(),
                else:
                    gaps    = substrings[org][block].count(".")
                    nongaps = blocklen-gaps
                    print cxpdr.sequences[org][offset:offset+nongaps].lower()+"-"*gaps,
            print org
            for block in range(0,maxblocks):
                print substrings[org][block],
            print org
    ########################################################################
    if verbose:
        for block in range(0,maxblocks):
            if substrings[substrings.keys()[0]][block].count("1") > 0: continue
            for org in substrings.keys():
                offset = sum([ substrings[org][i].count("1")+substrings[org][i].count("0") for i in range(0,block) ])
                blocklen = len(substrings[org][block])
                if Set(substrings[org][block]) == Set("1"):
                    print cxpdr.sequences[org][offset:offset+blocklen].upper(),
                elif Set(substrings[org][block]) == Set("0"):
                    print cxpdr.sequences[org][offset:offset+blocklen].lower(),
                else:
                    gaps    = substrings[org][block].count(".")
                    nongaps = blocklen-gaps
                    print cxpdr.sequences[org][offset:offset+nongaps].lower()+"-"*gaps,
                print substrings[org][block],
                print org
    ########################################################################
    if verbose:
        for org in substrings.keys():
            print org, "\t", 
            for block in range(0,maxblocks):
                print len(substrings[org][block]),
                if substrings[org][block].count("1") == 0:
                    print "(%s,%s)" % ( substrings[org][block].count('0'),substrings[org][block].count('.') ),
            print "\t\t", sum([ len(substrings[org][block])  for block in  range(0,maxblocks) ])
        print stw.lap()
    ########################################################################
    return substrings

# end of function cexpander2multiplealignment


def checkCBGs4omsrbordergaps(gsg,
    ignore_lsrcbg_boundaries = True,
    ignore_ksminx_cbgs=False, ignore_ks_cbgs=False,
    nonuniform_max_aa_length =\
            CBG_CEXPANDER_OMSRBORDERGAPS_NONUNIFORM_MAX_AA_LENGTH,
    nonuniform_aa_offset = CBG_CEXPANDER_OMSRBORDERGAPS_NONUNIFORM_AA_OFFSET,
    gap_size = CBG_CEXPANDER_OMSRBORDERGAPS_GAP_SIZE,
    max_bitscoreratio_threshold =\
            CBG_CEXPANDER_OMSRBORDERGAPS_MAX_BITSCORERATIO_THRESHOLD,
    verbose = False ):
    """
    @type  gsg: GeneStructureOfCodingBlockGraphs
    @param gsg: GeneStructureOfCodingBlockGraphs instance

    @type  ignore_lsrcbg_boundaries: Boolean
    @param ignore_lsrcbg_boundaries: what to do with a CBG that is neighboured
                                     by an lsrCBG?

    @type  ignore_ks_cbgs: Boolean
    @param ignore_ks_cbgs: if True, K(s) CBGs are skipped

    @type  ignore_ksminx_cbgs: Boolean
    @param ignore_ksminx_cbgs: if True, K(s-x) CBGs are skipped

    @type  nonuniform_max_aa_length: integer
    @param nonuniform_max_aa_length: maximum alowed length of the nonuniform
                                     stretch to process

    @attention: long nonuniform stretches often occur in low id% protein alignments!

    @type  nonuniform_aa_offset: integer
    @param nonuniform_aa_offset: area of the nonuniform stretch to check for gaps

    @type  gap_size: integer
    @param gap_size: continuous gap length to occur in the nonuniform_aa_offset
                     in order to shorten the CBG

    @type  max_bitscoreratio_threshold: float
    @param max_bitscoreratio_threshold: maximal bitscore ratio of Q vs. S slice
                                        to enforce a CBG shortening

    @type  verbose: Boolean 
    @param verbose: print debugging/intermediate information to STDOUT 

    @rtype:  Boolean ( or NoOverallMinimalSpanningRange exception )
    @return: status weather or not the CBG was shortened
    """
    # counters for number of changed/deleted CBGs
    cbgs_changed_cnt, cbgs_deleted_cnt = 0, 0

    # loop REVERSED because elements might be added/removed during this for loop
    for pos in range(len(gsg)-1,-1,-1):
        this = gsg.codingblockgraphs[pos]
        # ignore IGNORED and lsrCBGs graphs
        if this.IS_IGNORED:
            continue
        if ignore_lsrcbg_boundaries and this.__class__.__name__ ==\
        'LowSimilarityRegionCodingBlockGraph':
            continue
        if ignore_ks_cbgs and this.node_count() == self.EXACT_SG_NODE_COUNT: 
            continue
        if ignore_ksminx_cbgs and this.node_count() < self.EXACT_SG_NODE_COUNT:
            continue

        # check if this CBG has nonuniform cexpander tail(s)
        try:
            hasconsistency = this._cexpander.binarystring.count("1") >= 1
            has5Pomsrflaw  = this._cexpander.binarystring[0] == "0"
            has3Pomsrflaw  = this._cexpander.binarystring[-1] == "0"
        except:
            print "NO CORRECT CEXPANDER BINARYSTRING!!!!"
            print this
            this.printmultiplealignment()
            print "'%s'" % this._cexpander.binarystring
            print "'%s'" % this._cexpander.header
            print "'%s'" % this._cexpander
            continue

        if hasconsistency == False:
            # only zeros; this problem should have been solved
            # somewhere earlier, but not here!
            continue

        # check the 5p inconsistency for length & its previous CBG
        prev, omit5pside = None, False
        if has5Pomsrflaw:
            if this._cexpander.binarystring.find("1") >=\
            nonuniform_max_aa_length:
                omit5pside = True
            # if the prev cbg is a lsrCBG -> omit3pside = True
            elif ignore_lsrcbg_boundaries and pos > 0:
                prev = gsg.codingblockgraphs[pos-1]
                if prev.__class__.__name__ ==\
                'LowSimilarityRegionCodingBlockGraph':
                    omit5pside = True
            else:
                pass

        # check the 3p inconsistency for length & the next CBG
        next, omit3pside = None, False
        if has3Pomsrflaw:
            if len(this._cexpander.binarystring) -\
            this._cexpander.binarystring.rfind("1") >=\
            nonuniform_max_aa_length:
                omit3pside = True
            # if the next cbg is a lsrCBG -> omit3pside = True
            elif ignore_lsrcbg_boundaries and pos < len(gsg)-1:
                next = gsg.codingblockgraphs[pos+1]
                if next.__class__.__name__ ==\
                'LowSimilarityRegionCodingBlockGraph':
                    omit3pside = True
            else:
                pass


        # now call the cexpander_checkCBG4omsrbordergaps function
        try:
            status = cexpander_checkCBG4omsrbordergaps(
                    this, verbose = verbose,
                    omit5pside=omit5pside,
                    omit3pside=omit3pside,
                    nonuniform_aa_offset = nonuniform_aa_offset,
                    max_bitscoreratio_threshold = max_bitscoreratio_threshold,
                    gap_size = gap_size )
            # if status (==True), this CBG was succesfully shortened
            if status: cbgs_changed_cnt+=1
        except NoOverallMinimalSpanningRange:
            # NoOverallMinimalSpanningRange Exception;
            # that is the signal for deleting this CBG
            this.IS_IGNORED = True
        except ZeroUniformlyAlignedPositions:
            # due to optimization, the multiple alignment collapsed
            # that is the signal for deleting this CBG
            this.IS_IGNORED = True
        except:
            # unexpected exception -> raise!
            raise "UnExpectedException in checkCBGs4omsrbordergaps"

    # Done. Cleanup by actually removing the IS_IGNORED cbgs
    _dsGSG, _usGSG, etcGSG = gsg.separate_ds_and_us_gsg()
    cbgs_deleted_cnt = len(_dsGSG) + len(_usGSG) + len(etcGSG)

    if not ignore_lsrcbg_boundaries and (cbgs_changed_cnt or cbgs_deleted_cnt):
        LSRCBG_CORRECTED = gsg.correct_lsrcbgs_after_optimization(verbose=verbose)
        ################################################################
        if verbose: print "LSRCBG_CORRECTED:", LSRCBG_CORRECTED
        ################################################################

    # return the counters
    return ( cbgs_changed_cnt, cbgs_deleted_cnt )

# end of function checkCBGs4omsrbordergaps



def cexpander_checkCBG4omsrbordergaps(cbg,
    omit5pside = False, omit3pside = False,
    max_bitscoreratio_threshold =\
            CBG_CEXPANDER_OMSRBORDERGAPS_MAX_BITSCORERATIO_THRESHOLD,
    nonuniform_aa_offset = CBG_CEXPANDER_OMSRBORDERGAPS_NONUNIFORM_AA_OFFSET,
    gap_size = CBG_CEXPANDER_OMSRBORDERGAPS_GAP_SIZE,
    verbose = False):
    """
    Check the area directly around the OMSR of a CBG for non-uniform alignments

    @type  cbg: CodingBlockGraph
    @param cbg: CodingBlockGraph instance to optimize 

    @type  omit5pside: Boolean
    @param omit5pside: Do not process the 5' side (left) of the CBG

    @type  omit3pside: Boolean
    @param omit3pside: Do not process the 3' side (rigth) of the CBG

    @type  nonuniform_aa_offset: integer
    @param nonuniform_aa_offset: area of the nonuniform stretch to check for gaps

    @type  gap_size: integer
    @param gap_size: continuous gap length to occur in the nonuniform_aa_offset in order to shorten the CBG

    @type  max_bitscoreratio_threshold: float
    @param max_bitscoreratio_threshold: maximal bitscore ratio of Q vs. S slice to enforce a CBG shortening

    @type  verbose: Boolean 
    @param verbose: print debugging/intermediate information to STDOUT 

    @rtype:  Boolean ( or NoOverallMinimalSpanningRange or ZeroUniformlyAlignedPositions exception )
    @return: status weather or not the CBG was shortened
    """
    hasconsistency = cbg._cexpander.binarystring.count("1") >= 1
    has5Pomsrflaw  = cbg._cexpander.binarystring[0] == "0"
    PACBPS_CORRECTED = 0

    if not hasconsistency:
        # a priori error. CBGs must have at least a single Uniformly Aligned AA position
        raise ZeroUniformlyAlignedPositions

    if not omit5pside and (hasconsistency,has5Pomsrflaw) == (True,True):
        # start correction on the 5' side of the OMSR
        omsr = cbg.overall_minimal_spanning_range()
        replacements = {}
        ########################################################################
        if verbose:
            print "STARTING cexpander_checkCBG4omsrbordergaps 5p side"
            print cbg
            print "cexp::", cbg._cexpander.binarystring, cbg._cexpander.header
        ########################################################################
        for (currentkey,nodeQ,nodeS),pacbporf in cbg.pacbps.iteritems():
            # get slice of the pacbporf around the max(OMSR) query value 
            orgQ     = cbg.organism_by_node(nodeQ) 
            cexpQstr = cbg._cexpander.get_transferblock(orgQ).binarystring 
            endQpos  = min(omsr[nodeQ]) + cexpQstr.find("1")
            staQpos  = endQpos - nonuniform_aa_offset

            # get slice of the pacbporf around the max(OMSR) sbjct value 
            orgS     = cbg.organism_by_node(nodeS) 
            cexpSstr = cbg._cexpander.get_transferblock(orgS).binarystring 
            endSpos  = min(omsr[nodeS]) + cexpSstr.find("1")
            staSpos  = endSpos - nonuniform_aa_offset

            # correct staQpos if < pacbporf.orfQ.protein_startPY
            staQpos  = max([ pacbporf.orfQ.protein_startPY, staQpos ])
            editedQ  = staQpos != endQpos - nonuniform_aa_offset

            # correct staSpos if < pacbporf.orfS.protein_startPY
            staQpos  = max([ pacbporf.orfS.protein_startPY, staSpos ])
            editedS  = staSpos != endSpos - nonuniform_aa_offset

            if editedQ and editedS: 
                if ( endSpos - (staSpos + nonuniform_aa_offset) ) >=\
                (endQpos - (staQpos + nonuniform_aa_offset) ): 
                    # editing on Sbjct is gte as on Query -> take Sbjct 
                    (q,m,s,coords) = pacbporf.alignmentpart_by_sbjct( staSpos, endSpos ) 
                else: 
                    # other way around -> take Query 
                    (q,m,s,coords) = pacbporf.alignmentpart_by_query( staQpos, endQpos ) 
            elif editedS: 
                # take by sbjct coords 
                (q,m,s,coords) = pacbporf.alignmentpart_by_sbjct( staSpos, endSpos ) 
            else: 
                # unedited or edited Query -> take by query coords 
                (q,m,s,coords) = pacbporf.alignmentpart_by_query( staQpos, endQpos ) 

            # check minval of coords; CBGs at the far 5' end of the
            # input DNA sequence can get negative coords for their
            # non-existing Orf frontal STOPcodon (up to -3)
            if min(coords) < 0: continue

            # get bitscore-ratio of this Query/Sbjct slice
            (qS,qE,sS,sE)  = coords
            bitscoreratio = pacb.calculate_bitscoreratio(q,s)

            # get bitscore-ratio of this Query/Sbjct slice
            bitscoreratio = pacb.calculate_bitscoreratio(q,s)

            # if more gaps in this alignment slice then expected -> a pacbp split will follow
            # the slice is of size (2*omsr_offset)+1
            if q.find('-'*gap_size) >= 0:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                pos = q.find('-'*gap_size)
                while pos+gap_size < len(q) and q[pos+gap_size] == "-":
                    pos+=1
                splitpos = qS + pos + gap_size
                # correct splitpos by pacbp.query_start
                splitpos = splitpos - pacbp.query_start

            elif s.find('-'*gap_size) >= 0:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                pos = s.find('-'*gap_size)
                while pos+gap_size < len(s) and s[pos+gap_size] == "-":
                    pos+=1
                splitpos = sS + pos + gap_size
                # correct splitpos by pacbp.sbjct_start
                splitpos = splitpos - pacbp.sbjct_start

            elif bitscoreratio <= max_bitscoreratio_threshold:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                # correct for matches on the rigth of the match string
                splitpos = qE - ( len(m) - m.rfind(" ") - 1)
                # correct splitpos by pacbp.query_start
                splitpos = splitpos - pacbp.query_start

            else:
                ################################################################
                if verbose:
                    print nodeQ, nodeS, "'%s' '%s' '%s'" % (q,m,s), coords,
                    print "settings:", (nonuniform_aa_offset, gap_size),
                    print "bitsratio: %1.3f" % bitscoreratio
                ################################################################
                # not passing the cut-off for splitting this pacbp
                continue

            ####################################################################
            if verbose:
                print "5p,", nodeQ,nodeS, (q,m,s,coords),
                print "bitsratio: %1.3f (thr:%1.3f)" % (
                    bitscoreratio,max_bitscoreratio_threshold)
                print pacbp, "relative splitpos:", splitpos
                pacbp.print_protein(_linesize=120)
            ####################################################################

            # now split the pacbp on this position and recreate the pacbporf
            pacbpR = pacb.splitting.split_pacb_on_coordinates(pacbp,(
                        splitpos,splitpos),returnside='rigth')

            if pacbpR:
                newpacbporf = pacb.conversion.pacbp2pacbporf(
                                pacbpR,pacbporf.orfQ,pacbporf.orfS)
                newpacbporf.extend_pacbporf_after_stops()
                # store to replacements dict
                replacements[(currentkey,nodeQ,nodeS)] = newpacbporf
                ################################################################
                if verbose:
                    print pacbpR
                    pacbpR.print_protein(_linesize=120)
                    print newpacbporf
                ################################################################
                # increase counter for how much pacbps are corrected
                PACBPS_CORRECTED+=1

        # do the replacements of 5' PacbP corrections
        status = _update_cbg_with_pacbporf_replacements(cbg,replacements)
        if status == True:
            pass    # cbg succesfully updated; still an OMSR
        elif status == False:
            # raise a NoOverallMinimalSpanningRange Exception
            print "WARNING: NoOverallMinimalSpanningRange", cbg
            raise NoOverallMinimalSpanningRange, str(cbg)
        else:
            pass 

    # check (again!) if the is any consistency and if there is a 3' inconsistency
    hasconsistency = cbg._cexpander.binarystring.count("1") >= 1
    has3Pomsrflaw  = cbg._cexpander.binarystring[-1] == "0"

    if not hasconsistency:
        # due to 5' optimization, the complete CBG alignment collapsed!
        raise ZeroUniformlyAlignedPositions


    if not omit3pside and (hasconsistency,has3Pomsrflaw) == (True,True):
        # start correction on the 3' side of the OMSR
        omsr = cbg.overall_minimal_spanning_range()
        replacements = {}
        ########################################################################
        if verbose:
            print "STARTING cexpander_checkCBG4omsrbordergaps 3p side"
            print cbg, "\ncexp::", cbg._cexpander.binarystring,
            print cbg._cexpander.header
        ########################################################################

        for (currentkey,nodeQ,nodeS),pacbporf in cbg.pacbps.iteritems():
            # get slice of the pacbporf around the max(OMSR) query value
            orgQ     = cbg.organism_by_node(nodeQ)
            cexpQstr = cbg._cexpander.get_transferblock(orgQ).binarystring
            staQpos  = max(omsr[nodeQ]) - ( len(cexpQstr) - cexpQstr.rfind("1") )
            endQpos  = staQpos + nonuniform_aa_offset

            # get slice of the pacbporf around the max(OMSR) sbjct value
            orgS     = cbg.organism_by_node(nodeS)
            cexpSstr = cbg._cexpander.get_transferblock(orgS).binarystring
            staSpos  = max(omsr[nodeS]) - ( len(cexpSstr) - cexpSstr.rfind("1") )
            endSpos  = staSpos + nonuniform_aa_offset
            
            # correct endQpos if > pacbporf.orfQ.protein_endPY
            endQpos = min([ pacbporf.orfQ.protein_endPY, endQpos ])
            editedQ = endQpos != staQpos + nonuniform_aa_offset

            # correct endSpos if > pacbporf.orfQ.protein_endPY
            endSpos = min([ pacbporf.orfS.protein_endPY, endSpos ])
            editedS = endSpos != staSpos + nonuniform_aa_offset
            
            if editedQ and editedS:
                if ( endSpos - (staSpos + nonuniform_aa_offset) ) >=\
                (endQpos - (staQpos + nonuniform_aa_offset) ):
                    # editing on Sbjct is gte as on Query -> take Sbjct
                    (q,m,s,coords) = pacbporf.alignmentpart_by_sbjct( staSpos, endSpos )
                else:
                    # other way around -> take Query
                    (q,m,s,coords) = pacbporf.alignmentpart_by_query( staQpos, endQpos )
            elif editedS:
                # take by sbjct coords
                (q,m,s,coords) = pacbporf.alignmentpart_by_sbjct( staSpos, endSpos )
            else:
                # unedited or edited Query -> take by query coords
                (q,m,s,coords) = pacbporf.alignmentpart_by_query( staQpos, endQpos )

            # get bitscore-ratio of this Query/Sbjct slice
            (qS,qE,sS,sE)  = coords
            bitscoreratio = pacb.calculate_bitscoreratio(q,s)

            # if more gaps in this alignment slice then expected -> a pacbp split will follow
            # the slice is of size (2*omsr_offset)+1
            if q.find('-'*gap_size) >= 0:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                splitpos = qS + q.find('-'*gap_size)
                # correct splitpos by pacbp.query_start
                splitpos = splitpos - pacbp.query_start

            elif s.find('-'*gap_size) >= 0:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                splitpos = sS + s.find('-'*gap_size)
                # correct splitpos by pacbp.sbjct_start
                splitpos = splitpos - pacbp.sbjct_start

            elif bitscoreratio <= max_bitscoreratio_threshold:
                # convert (back) to pacbp and obtain position where to split
                pacbp = pacb.conversion.pacbporf2pacbp(pacbporf)
                # correct for matches on the left of the match string
                splitpos = qS + m.find(" ")
                # correct splitpos by pacbp.query_start
                splitpos = splitpos - pacbp.query_start

            else:
                ################################################################
                if verbose:
                    print nodeQ, nodeS, "'%s' '%s' '%s'" % (q,m,s), coords,
                    print "settings:", (nonuniform_aa_offset, gap_size), 
                    print "bitsratio: %1.3f" % bitscoreratio
                ################################################################
                # not passing the cut-off for splitting this pacbp
                continue

            ####################################################################
            if verbose:
                print "3p,", nodeQ,nodeS, (q,m,s,coords),
                print "bitsratio: %1.3f (thr:%1.3f)" % (
                            bitscoreratio,max_bitscoreratio_threshold)
                print pacbp, "relative splitpos:", splitpos
                pacbp.print_protein(_linesize=120)
            ####################################################################

            # now split the pacbp on this position and recreate the pacbporf
            pacbpL = pacb.splitting.split_pacb_on_coordinates(
                            pacbp,(splitpos,splitpos),returnside='left')

            if pacbpL:
                newpacbporf = pacb.conversion.pacbp2pacbporf(
                                pacbpL,pacbporf.orfQ,pacbporf.orfS)
                newpacbporf.extend_pacbporf_after_stops()
                # store to replacements dict
                replacements[(currentkey,nodeQ,nodeS)] = newpacbporf
                ################################################################
                if verbose:
                    print pacbpL
                    pacbpL.print_protein(_linesize=120)
                    print newpacbporf
                ################################################################
                # increase counter for how much pacbps are corrected
                PACBPS_CORRECTED+=1


        # do the replacements of 3' PacbP corrections
        status = _update_cbg_with_pacbporf_replacements(cbg,replacements)
        if status == True:
            pass    # cbg succesfully updated; still an OMSR
        elif status == False:
            # raise a NoOverallMinimalSpanningRange Exception
            raise NoOverallMinimalSpanningRange, str(cbg)
        elif status == None:
            pass    # no updates done at all
        else:
            # NOT POSSIBLE -> status isa NoneBoolean
            pass

    ####################################################################
    if verbose and PACBPS_CORRECTED:
        print "REPLACEMENTS DONE:", PACBPS_CORRECTED, "omit5pside:",
        print omit5pside, "omit3pside", omit3pside
        print cbg
        cbg.printmultiplealignment()
    ####################################################################


    # return if there is something improved
    if PACBPS_CORRECTED: return True
    else:                return False

# end of function cexpander_checkCBG4omsrbordergaps


def cexpander_checkCBGs4inframeintrons(gsg):
    """
    Check the CBGs in the GSG for non-uniform regions that represent inframe intron(s)
    """
    pass

# end of function cexpander_checkCBGs4inframeintrons

def cexpander_checklsrCBGs4inframeintrons(gsg):
    """
    Check the lsrCBGs in the GSG for non-uniform regions that represent inframe intron(s)

    aligned_intron_min_aa_length=ALIGNED_INTRON_MIN_AA_LENGTH,verb

        return codingblock_splitting.potentially_contains_aligned_intron(
                self,gap_size=gap_size,length_discrepancy=length_discrepancy,
                intron_window_min_similarity_score=intron_window_min_similarity_score,
                dif_potentially_observed_vs_expected=dif_potentially_observed_vs_expected,
                MIN_TOTAL_PSSM_INFRAME_INTRON=MIN_TOTAL_PSSM_INFRAME_INTRON
                )
        # in codingblock_splitting.potentially_contains_aligned_intron,
        # an actual search for introns is performed

    """

    pass

# end of function cexpander2lsrCBGs
