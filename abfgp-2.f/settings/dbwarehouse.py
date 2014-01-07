DBWAREHOUSE_PATH                    = '/home/avdb/storage/abfgp/dbwarehouse21122009'
DBWAREHOUSE_MAXIMAL_NUM_LOCI        = 50 
DBWAREHOUSE_MINIMAL_NUM_LOCI        = 2 
DBWAREHOUSE_DEFAULT_MAX_NUM_LOCI    = 15
DBWAREHOUSE_MINIMAL_OVERLAP_RATIO   = 0.700 # 50% mentioned in presentation!
DBWAREHOUSE_MINIMAL_BITSCORE_RATIO  = 0.100 # 0.08 mentioned in presentation!

DBWAREHOUSE_SEARCH_METHOD           = 'SAFEORTHOLOGS'
                                             # on of BDBH SAFEORTHOLOGS SIMILARITY 
                                             # SAFEORTHOLOGS sets DBWAREHOUSE_ALLOW_PARALOGS to False
DBWAREHOUSE_SEARCH_METHOD_LIST      = ['BDBH', 'SAFEORTHOLOGS', 'SIMILARITY', 'HOMOLOGS']
DBWAREHOUSE_SAFEORTHOLOGS_RATIO     = 0.90   # ratio between best and 2th best bitscore ratio
DBWAREHOUSE_ALLOW_PARALOGS          = False  # Boolean
DBWAREHOUSE_MAXIMAL_LENGTH_RATIO    = 1.25   # length ratio between two proteins


DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MINVAL = 0.03
DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MAXVAL = 0.95
DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MINVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MINVAL )
DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MAXVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_MINIMAL_BITSCORE_RATIO_MAXVAL )


DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MINVAL = 0.20 #0.50
DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MAXVAL = 0.95
DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MINVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MINVAL )
DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MAXVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_MINIMAL_OVERLAP_RATIO_MAXVAL )


DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MINVAL = 0.75
DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MAXVAL = 0.98
DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MINVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MINVAL )
DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MAXVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_SAFEORTHOLOGS_RATIO_MAXVAL )


DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MINVAL = 1.05
DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MAXVAL = 5.00 # protA 5x as long protein protB!
DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MINVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MINVAL )
DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MAXVAL_STRREPR = "%1.2f" % (
    DBWAREHOUSE_MAXIMAL_LENGTH_RATIO_MAXVAL )


# dictionary that enables gene/protein identifier's to their
# short genome name
ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING = {
    'MGG':    'mgg',
    'FVEG':   'fveg',
    'FOXG':   'foxg',
    'FGSG':   'fgsg',
    'NCU':    'ncu',
    'ANID':   'anid',
    'CNAG':   'cnag',
    'CAWG':   'cawg',
    'MYCGR':  'mycgr',
    'AFL2G':  'afl',
    'ANIG':   'anig',
    'SS1G':   'ss',
    'BC1G':   'bc',
    'PTRG':   'ptr',
    'PGTG':   'pgt',
    'SNOG':   'sno',
    'AO090':  'aory',
    'BC1T':   'boci',
    'CFU':    'cfu',
    'CHEC5':  'cheC5',
    'CIMT':   'cim',
    'CIMG':   'cim',
    'CPSG':   'cps',
    'HCAG':   'hca',
    'MYCFI':  'mycfi',
    'VDAG':   'vda',
    'VDBG':   'vdb',
}

ABGP_ORGANISM2FULLSPECIESNAME_MAPPING = {
    'afla':     'Aspergillus flavus',
    'anid':	'Aspergillus nidulans', 				
    'anig':	'Aspergillus niger',
    'aory':	'Aspergillus oryzae',
    'bc':       'Botrytis cinerea', # OLD ALIAS NAME!!!
    'boci':	'Botrytis cinerea',
    'cawg': 	'Candida albicans',				
    'cfu':      'Cladosporium fulvum',
    'cheC5':	'Cochliobolus heterostrophus C5',
    'cim':	'Coccidioides immitis RS',	
    'cnag': 	'Cryptococcus neoformans grubii h99', 	
    'cps':	'Coccidioides posadasii  str. silveira',		
    'crypa':    'Cryphonectria parasitica',
    'dotse':    'Dothistroma septosporum',
    'fgsg': 	'Fusarium graminiarum',
    'foxg':	'Fusarium oxysporum',
    'fveg': 	'Fusarium verticillioides',
    'hca':	'Histoplasma capsulatum nam1',
    'hyspu':    'Hysterium pulicare',
    'lepmu':    'Leptosphaeria maculans',
    'mgg':	'Magnaporthe grisea', 				
    'mycfi':	'Mycosphaerella fijiensis',		
    'mycgr':	'Mycosphaerella graminicola', 		
    'ncu':	'Neurospora crassa', 				
    'necha':    'Nectria haematococca',
    'pgt':	'Puccinia graminis f.sp. tritici',
    'ptr':	'Pyrenophora tritici-repentis',
    'rhyru':    'Rhytidhysteron rufulum',
    'sepmu':    'Septoria musiva SO2202',
    'sno':	'Stagonospora nodorum',			
    'ss':       'Sclerotinia sclerotiorum', # OLD ALIAS NAME!!
    'ss1':      'Sclerotinia sclerotiorum',
    'triat':    'Trichoderma atroviride',
    'vda': 	'Verticillium dahliae VdLs.17',		
    'vdb':	'Verticillium albo-atrum VaMs.102',		
    }


def _get_organism_full_name(organism,truncate=False):
    """ """
    if ABGP_ORGANISM2FULLSPECIESNAME_MAPPING.has_key(organism):
        # single homologous/orthologous protein of this species
        fullname = ABGP_ORGANISM2FULLSPECIESNAME_MAPPING[organism]
    elif ABGP_ORGANISM2FULLSPECIESNAME_MAPPING.has_key(organism[0:-1]):
        # paralogous proteins present!
        fullname = ABGP_ORGANISM2FULLSPECIESNAME_MAPPING[organism[0:-1]]
    else:
        # not found in mapping dict; set truncate to False
        fullname = organism
        truncate = False

    # if truncation of organism name is asked for, do so
    if truncate:
        fullname    = fullname.split(" ")
        fullname[0] = fullname[0][0]+"."
        fullname    = " ".join(fullname)

    # return organism/gene's full organism name
    return fullname

# end of function _get_organism_full_name


