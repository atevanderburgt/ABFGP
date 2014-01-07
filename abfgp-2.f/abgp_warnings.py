"""
Warnings used in Aligment Based Gene Prediction

call a Warning like this:
print AbgpWarning(*args)

"""

class AbgpWarning:
    """ Basal ABGP Warning class; all Warnings can directly inherit from this one """
    def __init__(self,*args):
        self.message = " ".join([str(arg) for arg in args])
    def __str__(self):
        return "WARNING::%s(%s)" % (self.__class__.__name__,self.message)


class UnexpectedEventWarning(AbgpWarning):
    """ Generic, unexpected event. Please provide data/evidence in the arguments """
    pass

class RepetitiveProteinSequenceWarning(AbgpWarning):
    """ A (annotated gene's) protein sequence is recognized as repetitive/lowcomplexity """
    pass

class CodingExonNotFoundAsOrfWarning(AbgpWarning):
    """ A (annotated) gene's or unigene's coding exon is not mappable on a single ORF """
    pass

class PotentialSequenceErrorWarning(AbgpWarning):
    """ A potential sequence error; please provide data/evidence in the arguments """
    pass

class SmallAnnotatedExonWarning(AbgpWarning):
    """ A (annotated) gene or unigene contains a small coding exon """
    pass

class SmallAnnotatedFirstExonWarning(SmallAnnotatedExonWarning):
    """ A (annotated) gene or unigene contains a small FIRST coding exon """
    pass

class SmallAnnotatedFinalExonWarning(SmallAnnotatedExonWarning):
    """ A (annotated) gene or unigene contains a small FINAL coding exon """
    pass

class SmallAnnotatedIntronWarning(AbgpWarning):
    """ A (annotated) gene or unigenes contains a very small intron """

class DnaSequenceContainsNsymbolsWarning(AbgpWarning):
    """ A DNA sequence (on which the gene might be encoded) contains N nucleotides """
    pass

class NonCanonicalSpliceSiteWarning(AbgpWarning):
    """ A (annotated) gene or unigene contains non-canonical splice site(s) """
    pass

class ImplausibleSpliceSiteWarning(AbgpWarning):
    """ A (annotated) gene or unigene contains implausible splice site(s) (non GT/GC/AG)"""
    pass

class UniGeneStructureIsNotMappableOnOrfsWarning(AbgpWarning):
    """ The ORF of an unigene structure is not mappable on the orfs OR contains implausible splice site(s) (non GT/GC/AG)"""
    pass

class GeneStructureIsNotMappableOnOrfsWarning(AbgpWarning):
    """ The (annotated) gene structure is not mappable on the orfs OR contains implausible splice site(s) (non GT/GC/AG)"""
    pass

class IncompleteGeneStructureWarning(AbgpWarning):
    """ The (annotated) gene structure is incomplete (no CDS tracks and/or no start/stop codons"""
    pass

