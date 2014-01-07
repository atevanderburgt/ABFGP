from sets import Set


def order_graphlist_by_total_weight(gralist):
    """
    """
    # sort the subgraphs based on their total_weight
    tmp = []
    for sg in gralist:
        tmp.append( ( sg.total_weight(), sg ) )
    tmp.sort()
    tmp.reverse()
    return [ sg for (tw,sg) in tmp ]

# end of function order_graphlist_by_total_weight


def make_consensus_genestructure_from_compatible_pacb_graphs(subgraphs,orgsimgra,INPUT_SEQUENCE_LENGTH=None):
    """
    Loop over all pacb graphs and order them, resulting a (a kind of) ordered list of coding blocks.
    This ordered list roughly represents the genstructure, except for the modeling of splice sites etc.
    This is performed by starting with the graph with highest evidence (graph's total weight).
    A graph is checked for being 'correct': is there an ``overall minimal spanning range`` of
    a minimal required size.
    When a graph with lower evidence is conflicting with an already placed graph, it is abandoned!


    subgraphs   list of codingblock graphs to be ordered
    orgsimgra   OrganismSimilarityGraph (of highest scoring graph)
    """

    # update edge weights to minimal spanning range
    for sg in subgraphs:
        sg.update_edge_weights_by_minimal_spanning_range()

    # order the graphs by their total weight
    subgraphs = order_graphlist_by_total_weight(subgraphs)

    genestructure = [ subgraphs[0] ]
    rejected_graphs = []

    ###orgsimgra.organism_set_size()

    for sg in subgraphs[1:]:
        for i in range(0,len(genestructure)):
            acb = genestructure[i]
            tba1, tba2, tba3 = 0, 0, 0
            tbb1, tbb2, tbb3 = 0, 0, 0
            print_check_buffer = []
            for org in sg.organism_set():
                acb_range = acb.minimal_spanning_range(organism=org)
                sg_range  =  sg.minimal_spanning_range(organism=org)
                (a1,a2,a3), (b1,b2,b3) = rel_pos_sets(acb_range,sg_range)
                print_check_buffer.append( ( (a1,a2,a3), (b1,b2,b3), org ) )
                (Ba1,Ba2,Ba3), (Bb1,Bb2,Bb3) = relpos2binaryrelpos( (a1,a2,a3), (b1,b2,b3) )
                # cumulative binary totals
                tba1+=Ba1
                tba2+=Ba2
                tba3+=Ba3
                tbb1+=Bb1
                tbb2+=Bb2
                tbb3+=Bb3

            # set all cumulative binary totals to binary
            (ba1,ba2,ba3), (bb1,bb2,bb3) = relpos2binaryrelpos( (tba1,tba2,tba3), (tbb1,tbb2,tbb3) )

            if ( (ba1,ba2,ba3), (bb1,bb2,bb3) ) == ( (1,0,0), (0,0,1) ):
                # this sg is BEHIND of current AcceptedCodingBlock (acb)
                # pass here, untill we are at the correct place in the genestructure
                pass
            elif ( (ba1,ba2,ba3), (bb1,bb2,bb3) ) == ( (0,0,1), (1,0,0) ):
                # this sg is in FRONT of current AcceptedCodingBlock (acb)
                genestructure.insert(i,sg)
                # break in order to get into next iteration
                break
            elif ( (tba1,tba3), (tbb1,tbb3) ) == ( (0,5), (5,0) ):
                ###print "case new XX", "IGNORED", sg.get_nodes(), len(genestructure), (tba1,tba2,tba3), (tbb1,tbb2,tbb3)
                for ( (ca1,ca2,ca3), (cb1,cb2,cb3), org ) in print_check_buffer:
                    if ca3-ca2 <= 0 or cb1-cb2 <= 0:
                        print "case new XX BUT IGNORED", "IGNORED", sg.get_nodes(), len(genestructure), (tba1,tba2,tba3), (tbb1,tbb2,tbb3)
                        rejected_graphs.append( sg )
                        "\n".join([ str(item) for item in print_check_buffer])
                        break
                else:
                    # okay, all organisms as we expected!
                    genestructure.insert(i,sg)
                    ### print "INSERTED new XX", sg.get_nodes(), len(genestructure)
                    break
            elif ( (tba1,tba3), (tbb1,tbb3) ) == ( (5,0), (0,5) ):
                for ( (ca1,ca2,ca3), (cb1,cb2,cb3), org ) in print_check_buffer:
                    if ca1-ca2 <= 0 or cb3-cb2 <= 0:
                        print "case new YY BUT IGNORED", "IGNORED", sg.get_nodes(), len(genestructure), (tba1,tba2,tba3), (tbb1,tbb2,tbb3)
                        rejected_graphs.append( sg )
                        "\n".join([ str(item) for item in print_check_buffer])
                        break
                else:
                    # okay, all organisms as we expected!
                    pass
            elif ( (tba3,), (tbb1,tbb3) ) == ( (0,), (0,5) ):
                # this sg is in FRONT of current AcceptedCodingBlock (acb)
                ### print "to be APPENDED NEWEST CASE", sg.get_nodes(), len(genestructure),  len(genestructure), (tba1,tba2,tba3), (tbb1,tbb2,tbb3)
                pass
            elif ( (tba1,tba3), (tbb3,) ) == ( (0,5), (0,) ):
                # this sg is BEHIND of current AcceptedCodingBlock (acb)
                ### print "NEWEST CASE BEHIND", ( (tba1,tba3), (tbb1,tbb3) ), sg.get_nodes(), len(genestructure)
                pass
            elif ( (tba1,tba2+tba3), (tbb1,tbb3) ) == ( (0,5), (5,0) ):
                for ( (ca1,ca2,ca3), (cb1,cb2,cb3), org ) in print_check_buffer:
                    if cb1-cb2 <= 0:
                        print "case new OVERLAPPING BEFORE a XX BUT IGNORED", "IGNORED", sg.get_nodes(), len(genestructure), (tba1,tba2,tba3), (tbb1,tbb2,tbb3)
                        rejected_graphs.append( sg )
                        ### print "\n".join([ str(item) for item in print_check_buffer])
                        break
                else:
                    # okay, all organisms as we expected!
                    genestructure.insert(i,sg)
                    ### print "INSERTED new XX", sg.get_nodes(), len(genestructure)
                    break
            elif ( (tba1,tba3), (tbb1+tbb2,tbb3) ) == ( (0,5), (5,0) ):
                for ( (ca1,ca2,ca3), (cb1,cb2,cb3), org ) in print_check_buffer:
                    if ca3-ca2 <= 0:
                        print "case new OVERLAPPING BEFORE b XX BUT IGNORED", "IGNORED", sg.get_nodes(), len(genestructure), (tba1,tba2,tba3), (tbb1,tbb2,tbb3)
                        rejected_graphs.append( sg )
                        ### print "\n".join([ str(item) for item in print_check_buffer])
                        break
                else:
                    # okay, all organisms as we expected!
                    genestructure.insert(i,sg)
                    ### print "INSERTED new XX", sg.get_nodes(), len(genestructure)
                    break


            elif ( (tba1,tba3), (tbb1,tbb2+tbb3) ) == ( (5,0), (0,5) ):
                # this sg is BEHIND of current AcceptedCodingBlock (acb)
                ### print "NEWEST CASE a OVERLAPPING BEHIND", ( (tba1,tba3), (tbb1,tbb3) )
                pass
            elif ( (tba1+tba2,tba3), (tbb1,tbb3) ) == ( (5,0), (0,5) ):
                # this sg is BEHIND of current AcceptedCodingBlock (acb)
                ### print "NEWEST CASE b OVERLAPPING BEHIND", ( (tba1,tba3), (tbb1,tbb3) )
                pass

            else:
                # a more complex case
                # IGNORE FOR NOW...
                print "IGNORED", sg.get_nodes(), len(genestructure), (tba1,tba2,tba3), (tbb1,tbb2,tbb3)
                rejected_graphs.append( sg )
                ### print "\n".join([ str(item) for item in print_check_buffer])
                break
        else:
            # okay, just place behind it!
            genestructure.append(sg)
            ### print "APPENDED", sg.get_nodes(), len(genestructure)


    # get rid of all that are to near the edges...
    # base it on flanking sequence, known protein sequence, etc...
    delthese = []
    for i in range(0,len(genestructure)):
        sg = genestructure[i]
        if max(sg.maximal_spanning_range(organism='mgg')) < 350:
            delthese.append(i)
        if INPUT_SEQUENCE_LENGTH:
            threshold = INPUT_SEQUENCE_LENGTH/3 - 350
            if min(sg.maximal_spanning_range(organism='mgg')) > threshold:
                delthese.append(i)

    delthese.reverse()
    for delthis in delthese:
        print "#$#$#$# out of known sequence length range"
        poppedsg = genestructure.pop(delthis)
        rejected_graphs.append( poppedsg )


    # set the order to the objects
    for i in range(0,len(genestructure)):
        genestructure[i].order = ( i+1, len(genestructure) )

    return genestructure, rejected_graphs 

# end of function make_consensus_genestructure_from_compatible_pacb_graphs


def rel_pos_sets(rA,rB):
    """
    How are 2 continious ranges of coordinates overlapping with each other?
    First argument is supposed to be LEFT / IN FRONT OF / BEFORE second argument.
    rA  Set     continious Set of coordinates; Set(range(A1,A2))
    rB  Set     continious Set of coordinates; Set(range(B1,B2))
    (a1,a2,a3)  tuple with tree absolute integer values >= 0
    (b1,b2,b3)  tuple with tree absolute integer values >= 0
    x1          outsticking coordinates on the LEFT  side of coordinate range
    x2          overlapping coordinates in both coordinate ranges
    x3          outsticking coordinates on the RIGTH side of coordinate range
    Example 1:
        rA = Set(list(range(10,31))
        rB = Set(list(range(28,31))
    Outcome 1:
        (a1,a2,a3) = (18,3,0)
        (b1,b2,b3) = (0,3,0)
    Example 2:
        rA = Set(list(range(10,31))
        rB = Set(list(range(6,28))
    Outcome 2:
        (a1,a2,a3) = (0,18,3)
        (b1,b2,b3) = (4,18,0)
    """
    (a1,a2,a3) = (0,0,0)
    (b1,b2,b3) = (0,0,0)
    if max(rA) < min(rB):
        a1,b3 = len(rA), len(rB)
    elif min(rA) > max(rB):
        a3,b1 = len(rA), len(rB)
    elif rA.symmetric_difference(rB) == Set():
        b2 = len(rA.intersection(rB))
        a2 = b2
    else:
        b2 = len(rA.intersection(rB))
        a2 = b2
        remA = rA.difference(rB)
        remB = rB.difference(rA)
        if not remA and remB:
            if len(remB) == max(remB) - min(remB) + 1:
                # a continious range
                if min(remB) < min(rA):
                    b1 = len(remB)
                else:
                    b3 = len(remB)
            else:
                # hmm.... a non-continious range ...
                intersect = rA.intersection(rB)
                b1 = min(intersect) - min(remB)
                b3 = max(remB) - max(intersect)
        elif not remB and remA:
            if len(remA) == max(remA) - min(remA) + 1:
                # a continious range
                if min(remA) < min(rB):
                    a1 = len(remA)
                else:
                    a3 = len(remA)
            else:
                # hmm.... a non-continious range ...
                intersect = rA.intersection(rB)
                a1 = min(intersect) - min(remA)
                a3 = max(remA) - max(intersect)
        else:
            if max(remA) < min(remB):
                a1,b3 = len(remA), len(remB)
            else:
                a3,b1 = len(remA), len(remB)
        
    # and return the relative positioning tuple
    return (a1,a2,a3), (b1,b2,b3)

# end of function rel_pos_sets


def relpos2binaryrelpos( (a1,a2,a3), (b1,b2,b3), overlap_offset = 3, percentual_overlap_offset = 1.5 ):
    """
    """
    aret, bret = [], []
    zerorange = range(0,overlap_offset+1)
    for a in (a1,a2,a3):
        if a > overlap_offset:      aret.append(1)
        elif a in zerorange:        aret.append(0)
        else:                       aret.append(-1)
    for b in (b1,b2,b3):
        if b > overlap_offset:      bret.append(1)
        elif b in zerorange:        bret.append(0)
        else:                       bret.append(-1)
    # check for the percentual overlap offset
    if not a3 and a1 and 100.0*(float(a2)/float(a1)) <= percentual_overlap_offset: aret[1] = 0
    if not a1 and a3 and 100.0*(float(a2)/float(a3)) <= percentual_overlap_offset: aret[1] = 0
    if not b3 and b1 and 100.0*(float(b2)/float(b1)) <= percentual_overlap_offset: bret[1] = 0
    if not b1 and b3 and 100.0*(float(b2)/float(b3)) <= percentual_overlap_offset: bret[1] = 0

    # and return the binary positioning tuple
    return ( tuple(aret), tuple(bret) )

# end of function relpos2binaryrelpos


