"""
Auxilary functions that must be placed somewhere....
"""

# Python Imports
from os.path import join as osPathJoin
from os.path import exists as osPathExists
from os.path import isdir as osPathIsdir
from os import remove as osRemove
from os import system as osSystem
from os import listdir as osListdir
from sys import exit as sysExit
from copy import deepcopy
from sets import Set

# Global Variable Imports
from settings.genestructure import ORF_IS_UNIGENE_LABEL

def _abgp_safe_identifier(name):
    """ """
    safename = deepcopy(name)
    for char in list( ":" ):
        safename = safename.replace(char,"x")
    for char in list( "$%^&*!@~<>+" ):
        safename = safename.replace(char,"")
    for char in list( "|#)(][}{" ):
        if safename.find(char) >= 0:
            safename = safename.split(char)[-1]
    for char in list( " _" ):
        if safename.find(char) >= 0:
            safename = safename.split(char)[-2]
    return safename

# end of function _abgp_safe_identifier


def debug(crossdata,input):
    """
    Debugging function feeded with `crossdata` and `input`.
    After each step, it can be checked what is in what-where-how.
    If nothing needed, return True immediately!
    """
    return True
    # below here, type functionality for debugging
    for k,v in crossdata[('fgsg', 'mgg')].iteritems():
        if type(v) != type({}): continue
        for (a,b,c,d),pacbp in v.iteritems():
            if c == 85 and d == 1:
                print k, (a,b,c,d), pacbp

# end of function debug


def _file_cleanup(fnamelist,include_directories=False):
    """ """
    for fname in fnamelist:
        if osPathExists(str(fname)):
            try:
                osRemove(str(fname))
            except OSError:
                if osPathIsdir(str(fname)) and include_directories:
                    osSystem("rm -rf %s"  % fname)
            except:
                # failed !? on cluster computing, this
                # could be the case when an identical filename
                # is created/owned by another user.
                # I suspect this might happen for formatdb.log ;-)
                pass

# end of function _file_cleanup


def _blastdb_cleanup(input):
    """ """
    fnamelist = ['formatdb.log']
    for organism in input.keys():
        if not input[organism]['blastdb']: continue
        fnamelist.append( input[organism]['blastdb'] )
        for extention in ['.phr','.pin','.psq']:
            fnamelist.append( input[organism]['blastdb'] + extention )
    # and perform generic file cleanup
    _file_cleanup(fnamelist) 

# end of function _blastdb_cleanup 


def recreate_global_informantdata(fromthis,target):
    """ """
    if fromthis.__class__.__name__ == 'dict':
        # make from input data dict structure
        ORGANISMS = fromthis.keys()
        ORGANISMS.sort()
        GENECOMBIS = [ (target,org) for org in ORGANISMS ]
        GENECOMBIS.remove((target,target))
        GENE_INFORMANT_SET    = Set()
        GENE_IDENTIFIER_SET   = Set()
        UNIGENE_INFORMANT_SET = Set()
        for informant in ORGANISMS:
            if fromthis[informant].has_key('is_unigene') and\
            fromthis[informant]['is_unigene']:
                UNIGENE_INFORMANT_SET.add(informant)
            else:
                GENE_IDENTIFIER_SET.add(informant)
                if informant != target:
                    GENE_INFORMANT_SET.add(informant)

    if fromthis.__class__.__name__ == 'PacbpCollectionGraph':
        ## make from PacbpCollectionGraph
        ## check if all nodes have PacbP(ORF)s
        ##for informant in fromthis.organism_set():
        ##    if informant == target: continue
        ##    nodes = fromthis.get_organism_nodes(informant)
        ##    for node in nodes:
        ##        if not fromthis.nodes[node]:
        ##            for n in nodes: print n, fromthis.nodes[n]
        ##            fromthis.del_node(node)
        ##            print "_delete_node:", node
        ##            continue
        ##        tnodes = list(fromthis.nodes[node])
        ##        for tnode in tnodes:
        ##            try:
        ##                pacbporfs = fromthis.get_pacbps_by_nodes(tnode,node)
        ##            except:
        ##                # no pacbporf(s) for this node -> delete it
        ##                fromthis.del_edge(tnode,node)
        ##                fromthis.del_node(node)
        ##                print "_delete_edge:", tnode,node

        # recreate lists
        ORGANISMS = list( fromthis.organism_set() )
        ORGANISMS.sort()
        GENECOMBIS = [ (target,org) for org in ORGANISMS ]
        # this should always be the case; in weird exceptions -> target is lost !?
        if (target,target) in GENECOMBIS: GENECOMBIS.remove((target,target))
        GENE_INFORMANT_SET    = Set()
        GENE_IDENTIFIER_SET   = Set()
        UNIGENE_INFORMANT_SET = Set()
        omitted_informants = []
        for informant in ORGANISMS:
            if informant == target:
                GENE_IDENTIFIER_SET.add(informant)
            else:
                try:
                    anypacbporf = fromthis.get_pacbps_by_organisms(target,informant)[0]
                    if hasattr(anypacbporf.orfS,ORF_IS_UNIGENE_LABEL):
                        UNIGENE_INFORMANT_SET.add(informant)
                    else:
                        GENE_IDENTIFIER_SET.add(informant)
                        GENE_INFORMANT_SET.add(informant)
                except IndexError:
                    # apparantly informant is removed from the PCG
                    omitted_informants.append(informant)
        # check if omitted_informants found
        for informant in omitted_informants:
            for node in fromthis.get_organism_nodes(informant):
                fromthis.del_node(node)
            # cleanup data structures
            ORGANISMS.remove(informant)
            try: GENE_INFORMANT_SET.remove(informant)
            except: pass
            try: GENE_IDENTIFIER_SET.remove(informant)
            except: pass
            try: UNIGENE_INFORMANT_SET.remove(informant)
            except: pass
            GENECOMBIS.remove((target,informant))

    # return data structures
    return ( ORGANISMS, GENECOMBIS,
             GENE_INFORMANT_SET, GENE_IDENTIFIER_SET, UNIGENE_INFORMANT_SET )

# end of function recreate_global_informantdata

def abgpsysexit(input,OPTIONS,message=""):
    """ """
    if input.has_key(OPTIONS.target):
        key = OPTIONS.target
    elif OPTIONS.target == None:
        key = None    
    else:
        # find by protein fref
        for k in input.keys():
            if input[k]['gldobj'].protein_fref() == OPTIONS.target:
                key = k
                break
        else:
            key = None

    # get filename to write to
    if key:
        fname = "%s.bailout.log" % input[key]['gldobj'].protein_fref()
    elif key == None and OPTIONS.filewithloci:
        hdr = OPTIONS.filewithloci.split("/")[-1].split(".")[0]
        fname = "%s.bailout.log" % hdr 
    elif key == None and OPTIONS.dirwithloci:
        hdr = OPTIONS.dirwithloci.split("/")[-1].split(".")[0]
        fname = "%s.bailout.log" % hdr
    else:
        fname = "%s.bailout.log" % ("Notargetapplied") 

    # clean the complete directory from all files
    _file_cleanup( osListdir(OPTIONS.outdir), include_directories=True )
        
    fname = osPathJoin(OPTIONS.outdir,fname)
    fh = open(fname,'w')
    fh.write(message+"\n")
    fh.close()
    # safely break out of the ABGP algorithm
    sysExit()

# end of function abgpsysexit
