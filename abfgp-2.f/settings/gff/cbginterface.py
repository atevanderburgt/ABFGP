################################################################################
# GFF Feature names for CodingBlockGraphInterface objects,
# ConsideredSpliceRanges and SpliceSites 
################################################################################
from settings.abgp import ABGP_PROGRAM_NAME, ABGP_VERSION 

# track for CodingBlockGraphInterface
GFF_CBGIF_FSOURCE           = ABGP_VERSION          # use ABGP_VERSION
GFF_CBGIF_FMETHOD           = 'CbgIF'
GFF_CBGIF_GCLASS            = 'CodingBlockGraphInterface'
GFF_CBGIF_GNAME             = ''                    # auto assigned
GFF_CBGIF_OUTPUT            = True


# track for ConsideredSpliceDonorRange
GFF_DONORRANGE_FSOURCE      = ABGP_VERSION          # use ABGP_VERSION
GFF_DONORRANGE_FMETHOD      = 'ConsideredSpliceDonorRange'
GFF_DONORRANGE_GCLASS       = 'SpliceDonorRange'
GFF_DONORRANGE_GNAME        = ''                    # auto assigned
GFF_DONORRANGE_OUTPUT       = True

# track for ConsideredSpliceAcceptorRange
GFF_ACCEPTORRANGE_FSOURCE   = ABGP_VERSION          # use ABGP_VERSION
GFF_ACCEPTORRANGE_FMETHOD   = 'ConsideredSpliceAcceptorRange'
GFF_ACCEPTORRANGE_GCLASS    = 'SpliceAcceptorRange'
GFF_ACCEPTORRANGE_GNAME     = ''                    # auto assigned
GFF_ACCEPTORRANGE_OUTPUT    = True

# track for SpliceDonor
GFF_DONOR_FSOURCE           = 'PSSM4FUNGI'
GFF_DONOR_FMETHOD           = 'SpliceDonor'
GFF_DONOR_GCLASS            = 'SpliceDoror'
GFF_DONOR_GNAME             = ''                    # auto assigned
GFF_DONOR_OUTPUT            = True

# track for SpliceAcceptor
GFF_ACCEPTOR_FSOURCE        = 'PSSM4FUNGI'
GFF_ACCEPTOR_FMETHOD        = 'SpliceAcceptor'
GFF_ACCEPTOR_GCLASS         = 'SpliceAcceptor'
GFF_ACCEPTOR_GNAME          = ''                    # auto assigned
GFF_ACCEPTOR_OUTPUT         = True

# track for AlignedDonorSite
GFF_ALIGNED_DONOR_FSOURCE   = 'PSSM4FUNGI-PACBP'
GFF_ALIGNED_DONOR_FMETHOD   = 'aligned_donor'
GFF_ALIGNED_DONOR_GCLASS    = 'AlignedDonorSite'
GFF_ALIGNED_DONOR_GNAME     = ''                    # auto assigned
GFF_ALIGNED_DONOR_OUTPUT    = True
GFF_ALIGNED_DONOR_REPORTBEST= 3

# track for AlignedAcceptorSite
GFF_ALIGNED_ACCEPTOR_FSOURCE= 'PSSM4FUNGI-PACBP'
GFF_ALIGNED_ACCEPTOR_FMETHOD= 'aligned_acceptor'
GFF_ALIGNED_ACCEPTOR_GCLASS = 'AlignedAcceptorSite'
GFF_ALIGNED_ACCEPTOR_GNAME  = ''                    # auto assigned
GFF_ALIGNED_ACCEPTOR_OUTPUT = True
GFF_ALIGNED_ACCEPTOR_REPORTBEST= 3


