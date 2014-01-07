################################################################################
# GFF Feature names for the known/annotated gene structure of the protein
################################################################################

from settings.abgp import ABGP_VERSION

# IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT
# *_FMETHOD strings are  *literaly* searched for in the gff of the annotated gene.
# And, the current gene's annotation is confirmed based on these tracks!
# IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT IMPORTANT

# main data for gene group
GFF_GENE_FSOURCE            = 'BroadGeneCall'
GFF_GENE_FMETHOD            = ''                    # defined in subclasses
GFF_GENE_GCLASS             = 'gene_id'
GFF_GENE_GNAME              = ''                    # auto assigned
GFF_GENE_OUTPUT             = True

# track for locus around the currently annotated protein sequence
# RECOMMENDED: include ~1500nt on both sides of the cds sequence
# this should enable to capture putatively missed exons
GFF_LOCUS_FSOURCE           = GFF_GENE_FSOURCE      # use GFF_GENE_FSOURCE
GFF_LOCUS_FMETHOD           = 'locus'               # free naming
GFF_LOCUS_GCLASS            = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_LOCUS_GNAME             = ''                    # auto assigned
GFF_LOCUS_OUTPUT            = GFF_GENE_OUTPUT       # use GFF_GENE_OUTPUT

# Track for CDS-`exons` of currently annotated protein sequence
# It expects a CDS-type of track. So, tracks that do NOT include 5' and 3' UTRs.
# in abgp_geneconfirmation.geneconfirmation() it is checked and corrected if
# the final exon track (accidentially) includes the stop codon itself (not all
# provided annotations are consistent...!), but when tracks of type EXON are
# used as input, it will barf on the 5' and 3' UTRs ...

# track for cds of currently annotated protein sequence
GFF_CDS_FSOURCE             = GFF_GENE_FSOURCE      # use GFF_GENE_FSOURCE
GFF_CDS_FMETHOD             = 'CDS'                 # free naming
GFF_CDS_GCLASS              = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_CDS_GNAME               = ''                    # auto assigned
GFF_CDS_OUTPUT              = GFF_GENE_OUTPUT       # use GFF_GENE_OUTPUT

# tracks for the exons - INCLUDING stop codons and potential UTRs
GFF_EXON_FSOURCE            = GFF_GENE_FSOURCE      # use GFF_GENE_FSOURCE
GFF_EXON_FMETHOD            = 'exon'                # free naming, BUT USE CONSISTENTLY
GFF_EXON_GCLASS             = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_EXON_GNAME              = ''                    # auto assigned
GFF_EXON_OUTPUT             = GFF_GENE_OUTPUT       # use GFF_GENE_OUTPUT

# track for intron of currently annotated protein sequence
GFF_INTRON_FSOURCE          = GFF_GENE_FSOURCE      # use GFF_GENE_FSOURCE
GFF_INTRON_FMETHOD          = 'intron'              # free naming
GFF_INTRON_GCLASS           = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_INTRON_GNAME            = ''                    # auto assigned
GFF_INTRON_OUTPUT           = False                 # introns not needed...

# track for start codon of currently annotated protein sequence
GFF_GENESTART_FSOURCE       = ''                    # auto assigned
GFF_GENESTART_FMETHOD       = 'start_codon'         # free naming
GFF_GENESTART_GCLASS        = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_GENESTART_GNAME         = ''                    # auto assigned
GFF_GENESTART_OUTPUT        = False
GFF_GENESTART_FMETHOD_LIST  = ['start_codon', 'startcodon', 'start']

# track for stop codon of currently annotated protein sequence
GFF_GENESTOP_FSOURCE        = ''                    # auto assigned
GFF_GENESTOP_FMETHOD        = 'stop_codon'          # free naming
GFF_GENESTOP_GCLASS         = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_GENESTOP_GNAME          = ''                    # auto assigned
GFF_GENESTOP_OUTPUT         = False
GFF_GENESTOP_FMETHOD_LIST   = ['stop_codon', 'stopcodon', 'stop']

# tracks for 5' and 3' UTR exons in the annotated gene
GFF_UTR5_FSOURCE            = ABGP_VERSION          # on-the-fly created
GFF_UTR5_FMETHOD            = 'utr5exon'            # free naming
GFF_UTR5_GCLASS             = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_UTR5_GNAME              = ''                    # auto assigned
GFF_UTR5_OUTPUT             = False

GFF_UTR3_FSOURCE            = ABGP_VERSION          # on-the-fly created
GFF_UTR3_FMETHOD            = 'utr3exon'            # free naming
GFF_UTR3_GCLASS             = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_UTR3_GNAME              = ''                    # auto assigned
GFF_UTR3_OUTPUT             = False

# track for transcript data => rarely available
GFF_TRANSCRIPT_FSOURCE      = GFF_GENE_FSOURCE      # use GFF_GENE_FSOURCE
GFF_TRANSCRIPT_FMETHOD      = 'transcript'          # free naming
GFF_TRANSCRIPT_GCLASS       = GFF_GENE_GCLASS       # use GFF_GENE_GCLASS
GFF_TRANSCRIPT_GNAME        = ''                    # auto assigned
GFF_TRANSCRIPT_OUTPUT       = False
GFF_TRANSCRIPT_GCLASS_LIST  = ['transcript_id', 'transcriptId']


# list which contains the usual annotated gene elements
GFF_GENE_ELEMENT_FMETHOD_LIST = []
GFF_GENE_ELEMENT_FMETHOD_LIST.append(GFF_EXON_FMETHOD)
GFF_GENE_ELEMENT_FMETHOD_LIST.append(GFF_CDS_FMETHOD)
GFF_GENE_ELEMENT_FMETHOD_LIST.extend(GFF_GENESTART_FMETHOD_LIST)
GFF_GENE_ELEMENT_FMETHOD_LIST.extend(GFF_GENESTOP_FMETHOD_LIST)

