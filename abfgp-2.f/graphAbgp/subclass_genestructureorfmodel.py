"""
Functions joined in a class for Orfmodel creation used in GeneStructureOfCodingBlocks class for Aligment Based Gene Predictions
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# Imports
from sets import Set

class GeneStructureOrfModelFunctions:
    """
    """
    def create_genestructure_orfmodel(self,printit=False,ignore_ignored=True):
        """
        Create orf-id genestructure orfmodel from genestructure_graphs

        @type  printit: Boolean
        @param printit: True or False; True prints a nicely layouted version of the output

        Example of its printed output:

        ## genestructure model
        fveg    49     49      33      121     121     117     117
        mgg     39     39      76      76      103     103     103
        ncu     14     93      68      68      68      41      111
        foxg    -     115      109     49      49      58      58
        fgsg    13     26      120     120     120     67      67
        """
        orfs2org = {}
        printheader = []
        for pos in range(0,len(self)):
            sg = self.codingblockgraphs[pos]
            if ignore_ignored and sg.IS_IGNORED: continue  # skip IGNORED codingblocks
            printheader.append( sg.organism_set_size() )
            for (org_gene_id,orf_id) in sg.get_nodes():
                if not orfs2org.has_key(org_gene_id):
                    orfs2org[org_gene_id] = []
                    if pos > 0:
                        for i in range(0,pos):
                            # again, omit the IGNORED CBGs
                            if self.codingblockgraphs[i].IS_IGNORED: continue
                            # append adash for a missing organism/orf in this CBG
                            orfs2org[org_gene_id].append('-')
                # append this orf id to the organism key list
                orfs2org[org_gene_id].append( orf_id )
            # Check if all possible organisms/genes are dealth with
            # In case we have partial codingblocks, we have to insert
            # a dash at this position
            if len(orfs2org) > sg.node_count():
                #print True, sg.get_nodes(), sg.organism_set(), orfs2org.keys()
                for org_gene_id in orfs2org.keys():
                    if org_gene_id not in sg.organism_set():
                        orfs2org[org_gene_id].append('-')

        # check if we need to print it
        if printit:
            orgkeys = orfs2org.keys()
            orgkeys.sort()
            print "CBGR\t", "\t".join([ "(%s)" % cnt for cnt in printheader ])
            for k in orgkeys:
                print k,"\t", "\t".join([ str(id) for id in orfs2org[k] ])
        # return the structure
        return orfs2org

    # end of function create_genestructure_orfmodel


    def condensed_genestructure_orfmodel(self,orfs2org={},printit=False):
        """
        Create condensed orf-id genestructure orfmodel from genestructure_graphs
    
        @type  orfs2org: dictionary
        @param orfs2org: dictionary with organism/gene as key and a list of codingblock-orfids as value  
    
        @type  printit: Boolean
        @param printit: True or False; True prints a nicely layouted version of the output
    
        Example of (printed) input:
    
        ## genestructure model
        fveg    49      33      121     121     117     117
        mgg     39      76      76      103     103     103
        ncu     93      68      68      68      41      111
        foxg    115     109     49      49      58      58
        fgsg    26      120     120     120     67      67
    
        Translates this into (printed) output
    
        ## genestructure model
        fveg    49      33      121     117
        mgg     39      76      103
        ncu     93      68      41      111
        foxg    115     109     49      58
        fgsg    26      120     67
    
        """
        # if not given, get the orfs2org data structure
        if orfs2org == {}: orfs2org = self.create_genestructure_orfmodel()

        condensed = {}
        for k,vlist in orfs2org.iteritems():
            condensed[k] = [ vlist[0] ]
            for item in vlist[1:]:
                if item != condensed[k][-1]:
                    condensed[k].append( item )
                else:
                    pass
            # remove all the '-' signs (lacking org/orf in a CBG)
            while '-' in condensed[k]: condensed[k].remove('-')
        if printit:
            orgkeys = condensed.keys()
            orgkeys.sort()
            for k in orgkeys:
                print k,"\t", "\t".join([ str(id) for id in condensed[k] ])
        return condensed
    
    # end of function condensed_genestructure_orfmodel


    def create_scaffold_genestructure(self,organism=None,printit=False):
        """
        Build a ``genestructure scaffold``

        The genestructure scaffold is the path with
            a) the smallest number of ORFs (from any organism) that
            b) avoids introns when possible by
            c) tiling codingblocks on ORFs that lack an intron
        through this gene structure graph (from ATG to TGA).

        When there is a choice >1 orf, the strongest connected node is taken.

        Now by example to explain:

        Given the (assumed correct!) `create_genestructure_orfmodel()` output:
            fveg    33      33      121     121     117     117
            mgg     39      76      76      103     103     103
            ncu     93      68      68      68      41      111
            foxg    109     109     49      49      58      58
            fgsg    26      120     120     120     67      67
        This will result in any of below 2 options:
            fveg    33       33      -       -       -       -
            mgg      -       -       -      103     103     103
            ncu      -       -       -       -       -       -
            foxg    109     109      -       -       -       -
            fgsg     -      120     120     120      -       -
        Depending on which of both ORFs (fveg-33) or foxg-109) is the
        strongest connected node in these first 2 CBGs, the scaffold is obtained.
        Assuming it is fveg, the genestructre scaffold will become
                    33       33      -      103     103     103
                     -      120     120     120      -       -
        This scaffold represents a protein sequence that lacks introns.
        This scaffold can be used to align the sequence from individual
        organisms of CBGs against. In such way, discrepansies can be discovered. 

        When there is an intron junction that is not bridged by a single ORF
        in any organism, the scaffold structure will look like:
                     33     120     120     120      -       -
                      -      -       -      103     103     103

        The output will be list a 2-element tuples, with as elements in
        the tuple the nodes that represent the ORF. By example:
            [   ( ('fveg',33),      None    ),
                ( ('fgsg',120),     None    ),
                ( ('fgsg',120),     None    ),
                ( ('fgsg',120), ('mgg',103) ),
                (     None,     ('mgg',103) ),
                (     None,     ('mgg',103) ),
            ]
        """

        # get the genestructure of the orfmodel
        org2orflist = self.create_genestructure_orfmodel()

        datastruct        = [] 
        printheader       = []
        pos2posinself     = {}
        pos               = 0
        for posinself in range(0,len(self)):
            cbg = self.codingblockgraphs[posinself]
            if not cbg.IS_IGNORED:
                datastruct.append([])
                printheader.append( cbg.organism_set_size() )
                pos2posinself[pos] = posinself
                pos+=1
        boundary_required = [1]*len(datastruct)
        true_required     = [1]*len(datastruct)
        forfilled         = [0]*len(datastruct)

        for length in range(len(datastruct),0,-1):
            if length == 1:
                for pos in range(0,len(datastruct)):
                    if forfilled[pos] == 0:
                        true_required[pos] = 1
                        if pos >= 1:
                            true_required[pos-1] = 1
                        if pos < len(datastruct)-1:
                            true_required[pos+1] = 1

            for org, orflist in org2orflist.iteritems():
                ####while "-" in orflist: orflist.remove("-")
                while "-" in orflist:
                    # Replace dashes (absent Orf) by minus integer values
                    # This maintains length of the list, but each negative
                    # value is unique -> no stretches spanning multiple CBGs.
                    # And, negative id's are recognized lateron
                    dashpos = orflist.index("-")
                    orflist[dashpos] = -dashpos

                for theorfid in Set(orflist):
                    if theorfid < 0:
                        # this wash a dash -> ignore
                        continue
                    if orflist.count(theorfid) == length:
                        inscaffold = True
                        therange = []
                        for pos in range(0,len(orflist)):
                            if orflist[pos] == theorfid:
                                therange.append(pos)
                                if length == 1 and forfilled[pos] > 0:
                                    inscaffold = False
                                    break
                                if forfilled[pos] >= boundary_required[pos]:
                                    inscaffold = False
                                    break
                        if inscaffold:
                            for pos in therange:
                                datastruct[pos].append((org,theorfid))
                            if therange[0] != 0 and length > 1:
                                boundary_required[ therange[0] ] += 1
                                true_required[ therange[0] ] = 2
                            if therange[-1] != len(datastruct)-1 and length > 1:
                                boundary_required[ therange[-1] ] += 1
                                true_required[ therange[-1] ] = 2
            # now update required/forfilled
            for pos in range(0,len(datastruct)):
                forfilled[pos] = len(datastruct[pos])


        # correct false doublets caused by excess of options
        for pos in range(0,len(datastruct)):
            if forfilled[pos] < true_required[pos]:
                # correct `true_required`
                true_required[pos] = 1
            # ignore positions for which a single Orf suffices
            if true_required[pos] == 1: continue
            # check positions that need >1 Orf and are not first or second
            if pos > 0 and pos < len(datastruct)-1:
                prev = datastruct[pos-1]
                this = Set(datastruct[pos])
                next = datastruct[pos+1]
                if len(this.intersection(prev))>0 and len(this.intersection(next))>0:
                    pass
                else:
                    # wrong! correct this one
                    true_required[pos] = 1

        # list that gathers all nodes that are UNambigious
        # and, gather headers needed for printing (if requested for)
        unambigious = []
        for pos in range(0,len(datastruct)):
            if forfilled[pos] == true_required[pos]:
                unambigious.extend( datastruct[pos] )


        ### TODO DEBUG preting to be deleted...
        ###for d in datastruct:
        ###    print "DD:", d
        ###print "REQ-B", boundary_required
        ###print "REQ-T", true_required
        ###print "FFLld", forfilled
        ###print "UNAMBIGIOUS:", unambigious

        # now loop over all ambigious positions and gather
        # node weighted connectivity values
        allnodes2scores = {}
        # `posinself` is the position in self.codingblockgraphs
        # `pos` is the position in datastruct
        # both numbers can be different in case of IS_IGNORED graphs
        pos = 0
        for pos in range(0,len(datastruct)):
            posinself = pos2posinself[pos]
            sg = self.codingblockgraphs[posinself]
            if forfilled[pos] > true_required[pos]:
                for node in sg.get_nodes():
                    if node in datastruct[pos] and node not in unambigious:
                        # gather weighted connectivity scoring of this node
                        if allnodes2scores.has_key(node):
                            allnodes2scores[node].append( sg.get_node_weighted_connectivity(node) )
                        else:
                            allnodes2scores[node] = [ sg.get_node_weighted_connectivity(node) ]
            # if here, increase pos
            pos += 1

        # solve the ambigiouty by taking the strongest connected node
        # for ambigious positions
        for pos in range(0,len(datastruct)):
            this = datastruct[pos]
            prev = []
            next = []
            if pos > 0:
                prev = datastruct[pos-1]
            if pos < len(datastruct)-1:
                next = datastruct[pos+1]
            if forfilled[pos] > true_required[pos]:
                intersectPrev = Set(this).intersection(prev)
                intersectNext = Set(this).intersection(next)
                # in case there is no intersection at all (first or last nodes...)
                # reset the first intersect to current node list
                if not intersectPrev and not intersectNext:
                    intersectPrev = Set(this)
                # new combination of nodes on this position
                new_combi = []
                posinself = pos2posinself[pos]
                sg = self.codingblockgraphs[posinself]
                # append unambigious nodes first
                for node in datastruct[pos]:
                    if node in unambigious:
                        new_combi.append(node)
                for nodeCombi in [intersectPrev,intersectNext]:
                    if not nodeCombi: continue
                    evalscores = []
                    for node in sg.get_nodes():
                        if node in datastruct[pos] and node not in unambigious and node in nodeCombi:
                            evalscores.append( ( sum( allnodes2scores[node] ), node ) )
                    evalscores.sort()
                    # check if there is a best node selected (this should be the case)
                    if evalscores:
                        bestscore, bestnode = evalscores[-1]
                        if bestnode not in new_combi:
                            new_combi.append(bestnode)
                # and update this position in the datastruct
                datastruct[pos] = new_combi

        # add None values on empty positions
        for d in datastruct:
            if len(d) == 1: d.append(None)
            d = tuple(d)

        # print scaffold data if requested for
        if printit:
            orgkeys = org2orflist.keys()
            orgkeys.sort()
            print "CBGR\t", "\t".join([ "(%s)" % cnt for cnt in printheader ])
            for org in orgkeys:
                tmp = []
                for items in datastruct:
                    if len(items)==2:
                        (node1,node2) = items
                        if node1 and org in node1:   tmp.append( node1[1] )
                        elif node2 and org in node2: tmp.append( node2[1] )
                        else:                        tmp.append('-')
                    elif len(items)==0:
                        tmp.append('-')
                    else:
                        tmp.append('?')
                print org,"\t", "\t".join([ str(id) for id in tmp ])

        # and return this scaffold data structure
        return datastruct

    # end of function create_scaffold_genestructure

# end of class GeneStructureOrfModelFunctions
