"""
Class with functions accessing/checking CodingBlockGraphInterfaces of CodingBlockGraphs
in a GeneStructureOfCodingBlocks used in Alignment Based Gene Prediction
"""

# Module metadata
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

# ABGP Imports
from codingblockgraphinterface import CodingBlockGraphInterface

class CodingBlockGraphInterfaceFunctions:
    """ """ 
    def create_cbginterfaces(self,ignore_optimal=True,ignore_compatible=True,
        allow_phase_shift=False,
        allow_non_canonical=False,
        optimizetinyexoninterface=False,
        verbose=False):
        """
        (Re)create CBGInterface objects in between CBGs in this GSG

        @type  ignore_optimal: Boolean
        @param ignore_optimal: Once a CBGInterface is optimal, do not recreate it

        @type  ignore_compatible: Boolean
        @param ignore_compatible: Once a CBGInterface is compatible, do not recreate it

        @type  allow_phase_shift: Boolean
        @param allow_phase_shift: (re)create CBGInterfaces allowing a phase shift of splice sites

        @type  allow_non_canonical: Boolean
        @param allow_non_canonical: (re)create CBGInterfaces allowing non-canonical (donor) splice sites

        @type  optimizetinyexoninterface: Boolean
        @param optimizetinyexoninterface: do a quick optimization of ths cbgIF
            (non-canonical, short suitable splice site range etc

        @type  verbose: Boolean
        @param verbose: print status messages to STDOUT

        @rtype:  Integer
        @return: number of CBGInterfaceobjects that is (re)created
        """
        RECREATED_CNT = 0
        for pos in range(1,len(self)):
            cbgD, cbgA = self[pos-1], self[pos]
            CREATE_INTERFACE       = False
            has_interface_donor    = self.has_donor_cbginterface(cbgD)
            has_interface_acceptor = self.has_acceptor_cbginterface(cbgA)
            if has_interface_donor and has_interface_acceptor:
                # interface objects already exist; only (re)create when not ignore_optimal
                pass
                #if self.cbginterface_is_optimal_donor(cbgD) and self.cbginterface_is_optimal_acceptor(cbgA):
                #    if not ignore_optimal: CREATE_INTERFACE = True
                #elif self.cbginterface_is_compatible_donor(cbgD) and self.cbginterface_is_compatible_acceptor(cbgA):
                #    if not ignore_compatible: CREATE_INTERFACE = True
                #else:
                #    CREATE_INTERFACE = True
            elif has_interface_donor:
                CREATE_INTERFACE = True
            elif has_interface_acceptor:
                CREATE_INTERFACE = True
            else:
                CREATE_INTERFACE = True

            if CREATE_INTERFACE:
                cbgIF = CodingBlockGraphInterface(cbgD,cbgA)
                cbgIF.harvest_splice_sites(allow_phase_shift=allow_phase_shift,allow_non_canonical=allow_non_canonical)
                cbgIF.find_conserved_splice_sites()
                if optimizetinyexoninterface:
                    cbgIF.optimizetinyexoninterface()
                # and set the interface objects to the CBGs in GSG
                cbgD._CBGinterface3p = cbgIF
                cbgA._CBGinterface5p = cbgIF
                RECREATED_CNT+=1
                if verbose: print cbgIF
            else:
                if verbose: print cbgD._CBGinterface3p, "EXISTING"

        # set current first & last CBG as IS_FIRST and IS_LAST
        if len(self):
            self.codingblockgraphs[0].IS_FIRST = True
            self.codingblockgraphs[-1].IS_LAST = True

        # return counter for how much CBGInterfaces are recreated
        return RECREATED_CNT

    # end of function create_cbginterfaces


    def clear_central_cbginterfaces(self):
        """
        Whipe out all the central CbgInterface objects
        """
        for pos in self.cbgpositerator()[1:]:
            self.codingblockgraphs[pos]._CBGinterface5p = None
        for pos in self.cbgpositerator()[0:-1]:
            self.codingblockgraphs[pos]._CBGinterface3p = None

    # end of function clear_central_cbginterfaces


    def has_donor_cbginterface(self,cbg):
        """
        Has this (lsr)CBG already a donor CBGInterface object assigned?

        @type  cbg: CodingBlockGraph
        @param cbg: CodingBlockGraph instance

        @rtype:  Boolean
        @return: True or False weather donor CBGInterface exists
        """
        if hasattr(cbg,"_CBGinterface3p") and cbg._CBGinterface3p:
            return True
        else:
            return False

    # end of function has_donor_cbginterface


    def has_acceptor_cbginterface(self,cbg):
        """
        Has this (lsr)CBG already an acceptor CBGInterface object assigned?

        @type  cbg: CodingBlockGraph
        @param cbg: CodingBlockGraph instance

        @rtype:  Boolean
        @return: True or False weather acceptor CBGInterface exists
        """
        if hasattr(cbg,"_CBGinterface5p") and cbg._CBGinterface5p:
            return True
        else:
            return False

    # end of function has_acceptor_cbginterface


    def cbginterface_is_optimal_donor(self,cbg):
        """
        Is the donor site of this CBGInterface optimal?

        @type  cbg: CodingBlockGraph
        @param cbg: CodingBlockGraph instance

        @rtype:  Boolean 
        @return: True or False weather CBGInterface's donor of this CBG is optimal
        """
        if self.has_donor_cbginterface(cbg) and cbg._CBGinterface3p.is_optimal_donor():
            return True
        else:
            return False

    # end of function cbginterface_is_optimal_donor


    def cbginterface_is_optimal_acceptor(self,cbg):
        """
        Is the acceptor site of this CBGInterface optimal?

        @type  cbg: CodingBlockGraph
        @param cbg: CodingBlockGraph instance

        @rtype:  Boolean 
        @return: True or False weather the CBGInterface's acceptor of this CBG is optimal
        """
        if self.has_acceptor_cbginterface(cbg) and cbg._CBGinterface5p.is_optimal_acceptor():
            return True
        else:
            return False

    # end of function cbginterface_is_optimal_acceptor


    def cbginterface_is_compatible_donor(self,cbg):
        """
        Is the donor site of this CBGInterface compatible (to the acceptor of the next CBG)?

        @type  cbg: CodingBlockGraph
        @param cbg: CodingBlockGraph instance

        @rtype:  Boolean
        @return: True or False weather CBGInterface's donor of this CBG is compatible 
        """
        if self.has_donor_cbginterface(cbg) and cbg._CBGinterface3p.is_compatible():
            return True
        else:
            return False

    # end of function cbginterface_is_compatible_donor


    def cbginterface_is_compatible_acceptor(self,cbg):
        """
        Is the acceptor site of this CBGInterface compatible (to the donor site of the previous CBG)?

        @type  cbg: CodingBlockGraph
        @param cbg: CodingBlockGraph instance

        @rtype:  Boolean
        @return: True or False weather the CBGInterface's acceptor of this CBG is compatible 
        """
        if self.has_acceptor_cbginterface(cbg) and cbg._CBGinterface5p.is_compatible():
            return True
        else:
            return False

    # end of function cbginterface_is_compatible_acceptor


    def improve_interfaces_by_shifting_splicesites(self,verbose=False):
        """
        @type  verbose: Boolean
        @param verbose: print status messages to STDOUT

        @rtype:  Integer
        @return: number of donor/acceptor CBGInterfaceobjects that is altered
        """
        acceptors_altered, donors_altered = [], []
        for cbg in self.codingblockgraphs[1:]:
            acceptors_altered.append(
                cbg._CBGinterface5p.optimize_aligned_acceptorgraph(verbose=verbose)
                )
            donors_altered.append(
                cbg._CBGinterface5p.optimize_aligned_donorgraph(verbose=verbose)
                )

        # return number of altered aligned sites (count the Trues in the list)
        return acceptors_altered.count(True) + donors_altered.count(True)

    # end of function improve_interfaces_by_shifting_splicesites

# end of class CodingBlockGraphInterfaceFunctions
