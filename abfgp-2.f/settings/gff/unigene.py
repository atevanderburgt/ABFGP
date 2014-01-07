################################################################################
# GFF Feature names for an unigene belonging to this gene locus
################################################################################

# main data for unigene group
GFF_UNIGENE_FSOURCE          = 'GenomeThreader'
GFF_UNIGENE_FMETHOD          = 'UGexon'
GFF_UNIGENE_GCLASS           = 'UniGene'
GFF_UNIGENE_GNAME            = ''                   # auto assigned
GFF_UNIGENE_OUTPUT           = True

# track for coding exons of a unigene/est
GFF_UGEXON_FSOURCE          = GFF_UNIGENE_FSOURCE   # use GFF_UNIGENE_FSOURCE
GFF_UGEXON_FMETHOD          = 'UGexon'              # free naming
GFF_UGEXON_GCLASS           = GFF_UNIGENE_GCLASS    # use GFF_UNIGENE_GCLASS
GFF_UGEXON_GNAME            = ''                    # auto assigned
GFF_UGEXON_OUTPUT           = GFF_UNIGENE_OUTPUT    # use GFF_UNIGENE_OUTPUT

# track for 5utr exon part of a unigene/est
GFF_UG5UTREXON_FSOURCE      = GFF_UNIGENE_FSOURCE   # use GFF_UNIGENE_FSOURCE
GFF_UG5UTREXON_FMETHOD      = 'UTR5UGexon'          # free naming
GFF_UG5UTREXON_GCLASS       = GFF_UNIGENE_GCLASS    # use GFF_UNIGENE_GCLASS
GFF_UG5UTREXON_GNAME        = ''                    # auto assigned
GFF_UG5UTREXON_OUTPUT       = GFF_UNIGENE_OUTPUT    # use GFF_UNIGENE_OUTPUT

# track for 3utr exon part of a unigene/est
GFF_UG3UTREXON_FSOURCE      = GFF_UNIGENE_FSOURCE   # use GFF_UNIGENE_FSOURCE
GFF_UG3UTREXON_FMETHOD      = 'UTR3UGexon'          # free naming
GFF_UG3UTREXON_GCLASS       = GFF_UNIGENE_GCLASS    # use GFF_UNIGENE_GCLASS
GFF_UG3UTREXON_GNAME        = ''                    # auto assigned
GFF_UG3UTREXON_OUTPUT       = GFF_UNIGENE_OUTPUT    # use GFF_UNIGENE_OUTPUT

