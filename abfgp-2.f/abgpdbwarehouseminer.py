#!/usr/bin/env python2.4
from settings.dbwarehouse import *
from lib_stopwatch import StopWatch
import os, re

"""

# Example usage code

from abgpdbwarehouseminer import *
miner = AbgpDbwarehouseMiner()
genomedir = miner.identifier2genomedir('MGG_0001')
print genomedir
genomedir = miner.identifier2genomedir('MGG_00016')
print genomedir
locusdir = miner.identifier2locusdir('MGG_0016')
print locusdir
locusdir = miner.identifier2locusdir('MGG_00005')
print locusdir
locusdir = miner.identifier2locusdir('MGG_00005',genomedir='/home/avdb/data/fungalgenomes/mgg/')
print locusdir
miner.mine('MGG_00005')

from abgpdbwarehouseminer import *
miner = AbgpDbwarehouseMiner()
data = miner.mine('MGG_12531')
miner.ALLOW_PARALOGS = True
data = miner.mine('MGG_12531')
miner.SEARCH_METHOD  = "SIMILARITY"
data = miner.mine('MGG_12531')
miner.SEARCH_METHOD  = "SAFEORTHOLOGS"
data = miner.mine('MGG_12531')
miner.SAFEORTHOLOGS_RATIO = 0.999
data = miner.mine('MGG_12531')

"""

class AbgpDbwarehouseMiner:
    """
    AbgpDbwarehouseMiner class; functionality for mining orthologous/similar proteins from Dbwarehouse
    """
    def __init__(self,dbwarehouse_path=DBWAREHOUSE_PATH,
        minimal_overlap_ratio=DBWAREHOUSE_MINIMAL_OVERLAP_RATIO,
        minimal_bitscore_ratio=DBWAREHOUSE_MINIMAL_BITSCORE_RATIO,
        search_method=DBWAREHOUSE_SEARCH_METHOD,
        safeorthologs_ratio=DBWAREHOUSE_SAFEORTHOLOGS_RATIO,
        allow_paralogs=DBWAREHOUSE_ALLOW_PARALOGS,
        maximal_length_ratio=DBWAREHOUSE_MAXIMAL_LENGTH_RATIO,
        maximal_num_loci=DBWAREHOUSE_MAXIMAL_NUM_LOCI,
        minimal_num_loci=DBWAREHOUSE_MINIMAL_NUM_LOCI,
        optparseoptions=None,
        verbose=False,
        ):
        # set keyword arguments to DbwarehouseMiner object
        self.dbwarehouse_path=dbwarehouse_path
        self.genomedirs = self.get_genomedirectories()
        self.MINIMAL_OVERLAP_RATIO=minimal_overlap_ratio
        self.MINIMAL_BITSCORE_RATIO=minimal_bitscore_ratio
        self.SEARCH_METHOD=search_method
        self.SAFEORTHOLOGS_RATIO=safeorthologs_ratio
        self.ALLOW_PARALOGS=allow_paralogs
        self.MAXIMAL_LENGTH_RATIO=maximal_length_ratio
        self.maximal_num_loci=maximal_num_loci
        self.minimal_num_loci=minimal_num_loci
        self.genometags_to_use    = []
        self.genometags_to_ignore = []
        self.verbose=verbose

        # if optparseoptions object applied -> get from options!
        if optparseoptions:
            self.verbose                = optparseoptions.verbose
            self.dbwarehouse_path       = optparseoptions.dbwarehouse
            self.MINIMAL_OVERLAP_RATIO  = optparseoptions.minimal_overlap_ratio
            self.MINIMAL_BITSCORE_RATIO = optparseoptions.minimal_bitscore_ratio
            self.SEARCH_METHOD          = optparseoptions.search_mode
            self.SAFEORTHOLOGS_RATIO    = optparseoptions.safeorthologs_ratio
            self.ALLOW_PARALOGS         = optparseoptions.allow_paralogs
            self.MAXIMAL_LENGTH_RATIO   = optparseoptions.maximal_length_ratio
            self.maximal_num_loci       = optparseoptions.maximal_num_loci
            self.minimal_num_loci       = optparseoptions.minimal_num_loci
            if optparseoptions.genomedirs_to_use:
                self.genomedirs = []
                for fulldirname in optparseoptions.genomedirs_to_use:
                    if fulldirname[-1] == "/": fulldirname = fulldirname[0:-1]
                    self.genometags_to_use.append(os.path.basename(fulldirname))
                    self.genomedirs.append(fulldirname+"/")
            if optparseoptions.genomedirs_to_ignore:
                for fulldirname in optparseoptions.genomedirs_to_ignore: 
                    try:
                        if fulldirname[-1] != "/": fulldirname = fulldirname + "/"
                        self.genomedirs.remove(fulldirname)
                        fulldirname = fulldirname[0:-1]
                        self.genometags_to_ignore.append(os.path.basename(fulldirname))
                    except:
                        # not found -> ignore, no error message or warning yet.
                        pass

        # result list attributes of mined loci
        self._loci = []
        self._data = []

    # end of function __init__


    def get_genomedirectories(self):
        """ """
        return [ os.path.join(self.dbwarehouse_path,dirname.strip()) for dirname in os.popen('ls -F %s | grep "/" | grep -v "^_"' % self.dbwarehouse_path).readlines() ]

    # end of function get_genomedirectories


    def printsettings(self):
        """ Printer-friendly and human-readable representation of search settings """
        print "--search_mode:            %s"    % self.SEARCH_METHOD
        if self.SEARCH_METHOD == 'SAFEORTHOLOGS':
            print "--safeorthologs_ratio:    %1.2f"    % self.SAFEORTHOLOGS_RATIO
        else:
            print "--allow_paralogs:         %s"    % self.ALLOW_PARALOGS
        print "--minimal_overlap_ratio:  %1.2f" % self.MINIMAL_OVERLAP_RATIO
        print "--minimal_bitscore_ratio: %1.2f" % self.MINIMAL_BITSCORE_RATIO
        print "--maximal_length_ratio:   %1.2f" % self.MAXIMAL_LENGTH_RATIO

    # end of function printsettings


    def identifier2locusdir(self,identifier,genomedir=False):
        """
        Find identifier in AbgpDbWarehouse

        @rtype  string or False
        @return full path to AbgpGeneLocusDirectory or False
        """
        identifier = identifier.replace("'","").replace('"','').strip()
        if not identifier: return False
 
        if not genomedir: genomedir = self.identifier2genomedir(identifier)
        if not genomedir:
            print "NO GENOMEDIR!!", identifier
            return False
        ### find the locus of this identifier
        ###gclass = self.readgenomedirsettings(genomedir)['gclass']
        ###fname  = os.path.join(genomedir,"etc","*.geneloci.gff") 
        ###command = """grep '%s "%s"' %s | sed 's/.*abgp_locus [.]*/abgp_locus \1/' | awk '{ print $2 }' | tr -d '"' """ % (
        ###                gclass, identifier, fname )
        ###_locusname = os.popen(command).read().strip().replace('\x01','')
        ###if not _locusname:
        ###    print "NO LOCUSNAME!!", identifier, genomedir
        ###    return False
        # Just assume the identifier name is the locusname...
        locusname = os.path.join(genomedir,"loci",identifier)
        if os.path.isdir(locusname):   
            return locusname
        else:
            # lookup identifier in gene2transcript.tbd.gz file
            g2tfname = os.path.join(genomedir,"etc","gene2transcript.tbd.gz")
            command = """zgrep '%s' %s | awk '{ if ($2=="%s") { print $1 } }' | head -n 1""" % (
                identifier, g2tfname, identifier )
            locusname = os.popen(command).read().strip().replace('\x01','')
            if locusname:
                locusname = os.path.join(genomedir,"loci",locusname)
                if os.path.isdir(locusname):   
                    return locusname
                else:
                    print "NO LOCUSDIR!!", identifier, genomedir, "locusname::", locusname
                    return False
            else:
                print "NO LOCUSDIR!!", identifier, genomedir
                return False


    # end of function identifier2locusdir


    def identifier2genomedir(self,identifier):
        """
        Find identifier in AbgpDbWarehouse

        @rtype  string or False
        @return full path to genomedir in dbwarehouse or False
        """
        identifier = identifier.replace("'","").replace('"','').strip()
        if not identifier: return False
        # quick shortcut; do first letters of identifier correspond to a
        # known organism tag?
        searchdirs = []
        for k in ABGP_ORGANISM2FULLSPECIESNAME_MAPPING.keys():
            if identifier.find(k.upper()) == 0:
                gdpath = os.path.join(self.dbwarehouse_path,k)
                searchdirs.append(gdpath)
                break

        # extend all the genomedirs to searchdirs
        #searchdirs.extend(self.genomedirs)

        for directory in searchdirs:  #self.genomedirs:
            mapdict = {
                'id':    identifier,
                'fname': os.path.join(directory,"etc","gene2transcript.tbd.gz")
            }
            if not os.path.isfile(mapdict['fname']): continue
            command = """zgrep "%(id)s" %(fname)s | awk '{ if ($1=="%(id)s" || $2=="%(id)s") { print $0 } }' | wc -l""" % mapdict
            try:    linecount = int(os.popen(command).read())
            except: linecount = 0
            if linecount: return directory
        else:
            return False

    # end of function identifier2genomedir


    def readgenomedirsettings(self,dirname):
        """
        """
        settings = {}
        for line in open(os.path.join(dirname,'settings.txt')).readlines():
            if line[0] != "#" and line.strip():
                kw,arg = line.strip().split('==',1)
                settings[kw.strip()]=arg.strip()
        return settings

    # end of function readgenomedirsettings


    def mine(self,identifier,verbose=None):
        """ """
        # (re)set mined results to empty
        self._data = []
        self._loci = []

        # start timer
        stw = StopWatch("dbwarehouseMiner.mine('%s')" % identifier )
        if verbose: print stw.start()

        # find the current identifier in the warehouse
        identifier = identifier.replace("'","").replace('"','').strip()
        if not identifier: return False
        genomedir = self.identifier2genomedir(identifier)
        if not genomedir: return False

        # append the main/central locusdir to the loci
        locusdir = self.identifier2locusdir(identifier,genomedir=genomedir)
        if not locusdir: return False
        self._loci.append( locusdir )

        if verbose: print stw.lap(), "main locus identified"

        # now mine in the warehouse
        if self.SEARCH_METHOD != 'SIMILARITY':
            # set some column restraints as VERY strict (&&) i.s.o loose (||)
            column_restrain = "&&"
        else:
            column_restrain = "||"

        ####genomedirtag   = os.path.basename(os.path.split(genomedir)[0])
        genomedirtag   = os.path.basename(genomedir)
        blastarchpatAB = os.path.join(self.dbwarehouse_path,"_crossblastp","blast.%s_x_*.symmetrized" % (genomedirtag))
        blastarchpatBA = os.path.join(self.dbwarehouse_path,"_crossblastp","blast.*_x_%s.symmetrized" % (genomedirtag))

        basecommand = """ awk -F':' '{ print $1"\\t"$2 }' | awk """ +\
                """ '{ if (($5>=%1.3f %s $6>=%1.3f) && ($7>=%1.3f %s $8>=%1.3f) && """ % (
                    self.MINIMAL_OVERLAP_RATIO,
                    column_restrain,
                    self.MINIMAL_OVERLAP_RATIO,
                    self.MINIMAL_BITSCORE_RATIO,
                    column_restrain,
                    self.MINIMAL_BITSCORE_RATIO,
                    ) +\
                """ (($5/$6)<=%1.2f %s ($6/$5)<=%1.2f)) { print $0"\t"(($5+$6)*$4)/2 } }' """ % (
                    self.MAXIMAL_LENGTH_RATIO,
                    column_restrain,
                    self.MAXIMAL_LENGTH_RATIO,
                    ) +\
                """ | sort -gr -k 9 """

        # commands with grep and zgrep for *.symmetrized and *.symmetrized.gz files
        command_grep = """grep "%s" %s %s | sort -u | %s""" % (
                identifier,
                blastarchpatAB,blastarchpatBA,
                basecommand)
        command_zgrep = """zgrep "%s" %s %s | sort -u | %s""" % (
                identifier,
                blastarchpatAB+".gz",blastarchpatBA+".gz",
                basecommand)

        # run the grep command
        ci,co,ce = os.popen3(command_grep)
        ci.close()
        lines = co.readlines()
        co.close()
        ce.close()

        # run the zgrep command
        ci,co,ce = os.popen3(command_zgrep)
        ci.close()
        lines.extend( co.readlines() )
        co.close()
        ce.close()

        seentags = []
        ignoretags = []
        for line in lines:
            fname, idA, idB, bitscore, overlapA, overlapB, ratioA, ratioB, order = line.strip().split("\t")
            if fname.find(".symmetrized.gz") >= 0:
                # process the lines obtained with the zgrep command
                tagA,tagB = fname[0:fname.find(".symmetrized.gz")][fname.find("/blast.")+7:].split("_x_")
            else:
                # process the lines obtained with the (normal) grep command
                tagA,tagB = fname[0:fname.find(".symmetrized")][fname.find("/blast.")+7:].split("_x_")

            # ignore the line completely when a limitation on genomedirs is applied and valid
            if self.genometags_to_use:
                if not (tagA in self.genometags_to_use and tagB in self.genometags_to_use):
                    continue
            if self.genometags_to_ignore:
                if tagA in self.genometags_to_ignore or tagB in self.genometags_to_ignore:
                    continue

            # ignore this line when one of the tags are (in) ignoretags
            if tagA in ignoretags: continue
            if tagB in ignoretags: continue

            # swap tagA & tagB when the tag's are in reversed order
            # this is due to the dbwarehouse crossblastp files
            # blast.B_x_A.symmetrized.gz isa symbolic link to
            # blast.A_x_B.symmetrized.gz if B > A (in string order)
            ordered_tags = [ tagA, tagB ]
            ordered_tags.sort()
            if [ tagA, tagB ] != ordered_tags:
                # swap tagA & tagB 
                tagA,tagB = tagB,tagA


            if self.SEARCH_METHOD == 'HOMOLOGS':
                if self.ALLOW_PARALOGS:
                    pass
                else:
                    if tagA == tagB:
                        continue
                    if tagA in seentags and tagB in seentags:
                        continue
            elif self.SEARCH_METHOD == 'BDBH':
                if tagA == tagB and self.ALLOW_PARALOGS:
                    if tagA in [ tup[0] for tup in self._data ]:
                        continue # there is already a fine hit gathered
                    else:
                        pass
                else:
                    if tagA == tagB and not self.ALLOW_PARALOGS:
                        continue
                    if tagA in seentags and tagB in seentags:
                        continue


            elif self.SEARCH_METHOD == 'SAFEORTHOLOGS':
                if tagA == tagB:
                    # check if there is not a paralog in the identifier's species it self
                    # that is to close nearby this identifier (a hypothetical paralogue)
                    ratioA, ratioB = float(ratioA), float(ratioB)
                    if max([ratioA,ratioB]) > self.SAFEORTHOLOGS_RATIO:
                        # there is in its own genome a hypothetical paralogue!
                        # empty data and break out!
                        self._data = []
                        break
                    else:
                        continue

                elif tagA in seentags and tagB in seentags:
                    if idA == identifier:
                        ratio = float(ratioA)
                        thetag = tagB
                    else:
                        ratio = float(ratioB)
                        thetag = tagA
                    maxratio = self._getfromdata(self._data,thetag)[5]
                    if min([ratio/maxratio, maxratio/ratio]) > self.SAFEORTHOLOGS_RATIO:
                        # remove this tag from data -> ortholog assignment is not 100% shure!
                        self._removefromdata(self._data,thetag)
                        ignoretags.append(thetag)
                        continue
                    else:
                        continue
                else:
                    pass

            else:
                # mode similarity -> all hits are okay
                pass


            # append tags to seentags
            if tagA not in seentags: seentags.append(tagA)
            if tagB not in seentags: seentags.append(tagB)

            # if here, a similar protein is mined!
            # gather locusdir and similarity data
            bitscore = int(float(bitscore))
            overlapA = float(overlapA)
            overlapB = float(overlapB)
            ratioA   = float(ratioA)
            ratioB   = float(ratioB)
            #if idA == identifier:
            #    self._data.append(( tagB, idB, bitscore, overlapA, overlapB, ratioA, ratioB ))
            #    ###print line.strip()
            #    ###print "A", self._data[-1],"\n"
            #else:
            #    self._data.append(( tagA, idA, bitscore, overlapB, overlapA, ratioB, ratioA ))
            #    ###print line.strip()
            #    ###print "B", self._data[-1],"\n"

            if idA == identifier:
                self._data.append(( tagB, idB, bitscore, overlapA, overlapB, ratioA, ratioB ))
            elif idB == identifier:
                self._data.append(( tagA, idA, bitscore, overlapB, overlapA, ratioB, ratioA ))
            elif idA.find(identifier) == 0:
                self._data.append(( tagB, idB, bitscore, overlapA, overlapB, ratioA, ratioB ))
            elif idB.find(identifier) == 0:
                self._data.append(( tagA, idA, bitscore, overlapB, overlapA, ratioB, ratioA ))
            else:
                print "WHAT ELSE!?::", tagA,tagB,idA,idB,bitscore, overlapA, overlapB, ratioA, ratioB



        # remove the TEMPORARILY element in mode SAFEORTHOLOGS
        if self.SEARCH_METHOD == 'SAFEORTHOLOGS':
            self._removefromdata(self._data,genomedirtag)

        # order _data on bitscore
        tmpdata = []
        for item in self._data:
           tmpdata.append( ( item[2], item ) )
        tmpdata.sort()
        tmpdata.reverse()
        self._data = [ item for (s,item) in tmpdata ]

        print len(self._data), self.maximal_num_loci

        # remove _data elements when self.maximal_num_loci is exceeded
        if len(self._data) > self.maximal_num_loci -1:
            if (self.verbose and verbose==None) or verbose:
                # print the removed loci to screen
                print "# removed loci (%s): --maximal_num_loci (%s) exceeded" % (
                    len(self._data)-self.maximal_num_loci+1, self.maximal_num_loci )
                for tup in self._data[self.maximal_num_loci-1:]:
                    row = list( tup )
                    row.insert(0,genomedirtag)
                    row.insert(2,identifier)
                    print "\t".join([ str(elem) for elem in row ])
            # now actually remove the rows from _data
            # minus 1 is for the --identifier locus itself
            self._data = self._data[0:self.maximal_num_loci-1]


        # get the loci belonging to the mined similar proteins
        for ( tagB, idB, bitscore, overlapA, overlapB, ratioA, ratioB ) in self._data:
            tagBgenomedir = os.path.join(self.dbwarehouse_path,tagB)
            locusdir = self.identifier2locusdir(idB,genomedir=tagBgenomedir)
            if not locusdir: print "HEROOO...."
            self._loci.append( locusdir )

        # add genomedirtag and identifier to _data rows
        for i in range(0,len(self._data)):
            row = list( self._data[i] )
            row.insert(0,genomedirtag)
            row.insert(2,identifier)
            self._data[i] = tuple(row)

        if (self.verbose and verbose==None) or verbose:
            # print the results!
            print "# main (1th) and mined loci"
            for locus in self._loci:
                print locus
            print "# similarity data"
            #for ( tagB, idB, bitscore, overlapA, overlapB, ratioA, ratioB ) in self._data:
            #    print "\t".join([ str(elem) for elem in [genomedirtag, tagB, identifier, idB, bitscore, overlapA, overlapB, ratioA, ratioB ]])
            for row in self._data:
                print "\t".join([ str(elem) for elem in row ])
            print "# settings/options"
            print "seentags:  ", seentags
            print "ignoretags:", ignoretags
            print "use:       ", self.genometags_to_use
            print "ignore:    ", self.genometags_to_ignore
            print "# timing/performace"
            print stw.lap()

        return self._loci, self._data

    # end of function mine 


    def _removefromdata(self,data,tag):
        """
        """
        for pos in range(0,len(data)):
            if tag == data[pos][0]:
                data.pop(pos)
                break
        
    # end of function _removefromdata


    def _getfromdata(self,data,tag):
        """
        """
        for pos in range(0,len(data)):
            if tag == data[pos][0]:
                return data[pos]
        else:
            return False 
    # end of function _getfromdata 

# end of class AbgpDbwarehouseMiner



def _read_dbwarehouse_symmetrized_files(dbwarehouse_path,organism=None):
    """
    """
    symmetrized_files_primary   = []
    symmetrized_files_secondary = []
    for filename in os.listdir(os.path.join(dbwarehouse_path,"_crossblastp")):
        filename = os.path.join(dbwarehouse_path,"_crossblastp",filename)
        if filename.find("symmetrized") == -1:
            continue
        if not (filename.find("symmetrized") == len(filename) - 11 or\
        filename.find("symmetrized") == len(filename) - 14):
            continue
        parts = filename.split(".")
        orgQ,orgS = parts[1].split("_x_")
        if orgQ == orgS: continue
        orderedorgs = [ orgQ, orgS ]
        orderedorgs.sort()
        if organism and organism not in orderedorgs: continue
        if [ orgQ,orgS ] == orderedorgs:
            symmetrized_files_primary.append( filename )
        else:
            symmetrized_files_secondary.append( filename )

    return symmetrized_files_primary, symmetrized_files_secondary

# end of function _read_dbwarehouse_symmetrized_files


def remove_dbwarehouse_symlinks(dbwarehouse_path,organism=None):
    """
    """
    symmetrized_files_primary, symmetrized_files_secondary =\
        _read_dbwarehouse_symmetrized_files(dbwarehouse_path,organism=organism)
    for filenameS in symmetrized_files_secondary:
        os.remove(filenameS)

# end of function remove_dbwarehouse_symlinks


def recreate_dbwarehouse_symlinks(dbwarehouse_path,organism=None):
    """
    """
    symmetrized_files_primary, symmetrized_files_secondary =\
        _read_dbwarehouse_symmetrized_files(dbwarehouse_path,organism=organism)

    # loop over all the primary files and check if its
    # secondary file exists (asa symlink)
    for filenameP in symmetrized_files_primary:
        filenameS = mirror_symmetrized_filename(filenameP)
        if filenameS in symmetrized_files_secondary:
            # check size of the file; it should be enourmously smaller!
            if os.path.getsize(filenameP) > os.path.getsize(filenameS):
                pass
            else:
                # symlink converted to real file somewhere. Correct this
                os.remove(filenameS)
                # create a symlink
                os.system("ln -s %s %s" % (filenameP,filenameS))
            # remove filenameS from file list
            symmetrized_files_secondary.remove(filenameS)
        else:
            # create a symlink
            os.system("ln -s %s %s" % (filenameP,filenameS))
    # remove all the remaining secondary files;
    # *.symmetrized and *.symmetrized.gz can become obsolete
    # due to gzipping/gunzipping parts of the warehouse
    for filenameS in symmetrized_files_secondary:
        print "REMOVING %s" % filenameS

# end of function recreate_dbwarehouse_symlinks

symfile_pattern = re.compile("^(.*blast\.)(\w{2,8})(_x_)(\w{2,8})(\.symmetrized.*)$")

def mirror_symmetrized_filename(filename):
    """ """
    match = re.search(symfile_pattern,filename)
    parts = list( match.groups() )
    return "".join([ parts[0], parts[3], parts[2], parts[1], parts[4] ])

# end of function mirror_symmetrized_filename


if __name__ == "__main__":
    # use abgpdbwarehouseminer as standalone application
    import sys
    from lib_optparse import abgpoptparser
    from lib_optparse import generaloptions, validate_generaloptions
    from lib_optparse import dbwarehousesearchoptions, validate_dbwarehousesearchoptions
    from lib_optparse import dbwarehousesettingsoptions, validate_dbwarehousesettingsoptions
    from lib_optparse import dbwarehousefunctionoptions, validate_dbwarehousefunctionoptions
    from lib_optparse import dbwarehouseresultoptions, validate_dbwarehouseresultoptions  
 
    # construct the command line option parser
    parser = abgpoptparser()
    generaloptions(parser)
    dbwarehousesearchoptions(parser)
    dbwarehouseresultoptions(parser)
    dbwarehousesettingsoptions(parser)
    dbwarehousefunctionoptions(parser)

    # parse the command line & validate
    (options, args) = parser.parse_args()

    validate_dbwarehousefunctionoptions(parser,options) # to this one first! 
    validate_generaloptions(parser,options)
    validate_dbwarehousesettingsoptions(parser,options)
    validate_dbwarehousesearchoptions(parser,options)

    # construct the miner and mine for similar proteins
    miner = AbgpDbwarehouseMiner(optparseoptions=options)
    loci, similaritydata = miner.mine(options.identifier,verbose=options.verbose)
    options.loci = loci 

    # validate the miner results (which are updated into options attribure)
    status = validate_dbwarehouseresultoptions(parser,options)
    if not status:
        miner.printsettings()
        sys.exit()

    # create the outdir and cp all the loci to there
    if options.createoutputdir:
        # create output directory
        if os.path.isdir(options.outdir) and options.force == False:
            print "OUTPUT DIRECTORY %s ALREADY EXIST; overwrite with --force" % options.outdir
            sys.exit()
        elif os.path.isdir(options.outdir) and options.force == True:
            pass
        else:
            os.mkdir(options.outdir)
        # copy the genelocus directories
        for genelocusdir in loci:
            os.system("cp -r %s %s" % (genelocusdir,options.outdir) )
        # make similarity.csv file
        fh = open(os.path.join(options.outdir,"similarity.csv"),'w')
        # write command line to file
        fh.write("# "+" ".join([ item for item in sys.argv])+"\n")
        # write important arguments to file
        for opt in parser.option_groups[3].option_list:
            fh.write("# %s %s\n" % (opt, getattr(options,opt.dest)))
        for row in similaritydata:
            fh.write("\t".join([ str(elem) for elem in row ])+"\n")
        fh.close()


    if options.filewithloci:
        # print data to STDOUT which can serve as input filewithloci for ABFGP
        # print parameters
        for param in [
                'identifier',
                'genomedirs_to_ignore',
                'genomedirs_to_use',
                'search_mode',
                'allow_paralogs',
                'safeorthologs_ratio',
                'maximal_length_ratio',
                'maximal_num_loci',
                'minimal_bitscore_ratio',
                'minimal_num_loci',
                'minimal_overlap_ratio' ]:
            attrname = param + " "*25
            attrname = attrname[0:25]
            print "# %s: %s" % (attrname, getattr(options,param))

        # print similarity data
        print "# similarity data"
        for row in similaritydata:
            print "\t".join([ str(elem) for elem in row ])

        # print loci
        print "# target locus"
        print options.locus
        print "# informant loci"
        for genelocusdir in options.loci[1:]:
            print genelocusdir


    # print the command line options
    if options.verbose:
        print "# command line options"
        for attr in dir(options):
            if attr[0] == "_": continue
            value = getattr(options,attr)
            if str(value)[0:6] == "<bound": continue
            strattr = attr+" "*30
            print strattr[0:25], "\t", value



