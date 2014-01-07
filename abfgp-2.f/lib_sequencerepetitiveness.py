from pythonlibs.sequencerepetitiveness import *
from abgp_warnings import RepetitiveProteinSequenceWarning


def annotatedproteinsequencerepetitivenesscheck(input,verbose=True,threshold=0.175):
    """
    @type  input: dict
    @param input: input dictionary data structure

    @type  verbose: Boolean
    @param verbose: print discrepancies to STDOUT (True) or not (False, default)

    @rtype:  Boolean 
    @return: Is any of the annotated/unigene ORFs recognize as repetitive? 
    """
    IS_REPETITIVE = False
    for org,data in input.iteritems():
        # check if 'orfid-unigenestructure' key exists -> if no unigeneconfirmation()
        # function call has been done, key might be absent
        if input[org].has_key('orfid-unigenestructure') and input[org]['orfid-unigenestructure']:
            for id in input[org]['orfid-unigenestructure']:
                orf = input[org]['orfs'].get_orf_by_id(id)
                score = proteinsequencelowcomplexityscore(orf.protein_sequence)
                if score >= threshold:
                    IS_REPETITIVE = True
                    if verbose:
                        print RepetitiveProteinSequenceWarning("%1.2f" % score, "of unigene ORF ('%s',%s)" % (org,orf.id))


        if input[org]['orfid-genestructure']:
            for id in input[org]['orfid-genestructure']:
                orf = input[org]['orfs'].get_orf_by_id(id)
                score = proteinsequencelowcomplexityscore(orf.protein_sequence)
                if score >= threshold:
                    IS_REPETITIVE = True
                    if verbose:
                        print RepetitiveProteinSequenceWarning("%1.2f" % score, "of annotated ORF ('%s',%s)" % (org,orf.id))

    # return the IS_REPETITIVE status
    return IS_REPETITIVE

# end of function annotatedproteinsequencerepetitivenesscheck

