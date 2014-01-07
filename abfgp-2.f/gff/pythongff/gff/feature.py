"""Holds a GFF Feature
"""

import types,string

class FeatureInputError(Exception):
    pass

class FeatureComment(Exception):    
    pass

class Feature:
    ##constructor

    def __init__(self, theentries = "", version = "2.0"):
        
        if version == "2.0":
            # Version
            self.__version = version
            
            try:
                # If input is a string then convert to list.
                if type(theentries) == types.StringType:
                    theentries = string.strip(theentries)
                    entries = string.split(theentries)
		    #entries = theentries.strip().split()
                    if entries[0][0]!='#':
                        if theentries=="":  #create empty object
                            self.__seqname = ""
                            self.__source = ""
                            self.__feature = ""
                            self.__start = 0
                            self.__end = 0
                            self.__score = '.'
                            self.__strand = '.'
                            self.__frame = "."
                            self.__attribute = ""
                        elif 8 <= len(entries) <= 10:
                            # Seqname, Soucre, Feature
                            self.__seqname = entries[0]
                            self.__source  = entries[1]
                            self.__feature = entries[2]
                            
                            # Start and End
                            try:
                                start = string.atoi(entries[3])
                                end = string.atoi(entries[4])
                            except ValueError,s:
                                raise FeatureInputError("Start and end field should be integers.")          
                            self.__start   = int(entries[3])
                            self.__end     = int(entries[4])
                            if not 1 <= self.__start <= self.__end:
                                raise FeatureInputError("Start or end index out of range.")

                            # Score
                            if not entries[5] in ['.']:
                                try:
                                    self.__score = float(entries[5])
                                except ValueError:
                                    raise FeatureInputError("Score should be a floating point value.")
                            else:
                                self.__score = '.'                                
                                # Strand
                            if not entries[6] in ['+','-','.']:
                                raise FeatureInputError("Legal strand values are '+','-' or '.'.")
                            self.__strand = entries[6]
                            # Frame
                            if not entries[7] in ['0','1','2','.']:
                                raise FeatureInputError("Legal frame values are 0,1,2 or '.'.")
                            self.__frame = entries[7]
                            # Attribute
                            if len(entries) >= 9:
                                self.__attribute = entries[8]
                            else:
                                self.__attribute = ""
                            if len(entries)>=10:
                                self.__comment = entries[9]
                            else:
                                self.__comment = ""
                        else:
                            raise FeatureInputError("Illegal number of fields.")
                    else:
                        raise FeatureComment()
            except FeatureInputError, s:
                raise FeatureInputError(str(s))
        else:
            raise ValueError("Only supports version 2.0.")
                                
   ## not strand sensitive
    def __cmp__(self,other):
        if type(self)==types.NoneType or type(other)==types.NoneType:
            return -1
        if self.start()<=other.start():
            if self.start()<other.start():
                return -1
            else: #self.start = other.start
                if self.end()==other.end():
                    return 0
                elif self.end()<other.end():
                    return -1
                else:
                    return 1                
        else:
            return 1

    def cmp(self,other):
       return  self.__cmp__(other)
    
    ## access functions ##
    def seqname(self,name = ""):
        if name!="":
            self.__seqname = name
        return self.__seqname

    def source(self,name = ""):
        if name!="":
            self.__seqname = name
        return self.__source    
    
    def feature(self,name = ""):
        if name!="":
            self.__seqname = name
        return self.__feature

    def start(self,value = -1):
        if value>=0 and value<=self.end():
            self.__start = value
        return self.__start

    def end(self,value = -1):
        if value>=0 and value>=self.start():
            self.__end = value
        return self.__end

    def score(self,value = ""):
        if value!="":
            self.__score = value
        return self.__score
    
    def strand(self,value = ""):
        if value!="" and value in ['+','-','.']:
            self.__strand = value
        return self.__strand

    def frame(self,value = ""):
        if value!="":
            self.__frame = value
        return self.__frame

    def attribute(self,value = ""):
        if value!="":
            self.__attribute = value              
        return self.__attribute

    def comment(self,value = ""):
        if value!="":
            self.__comment = value
        return self.__comment
    

    ## string representation ##
    def __str__(self, tab = '\t'):
        s = self.__seqname    + tab + \
            self.__source     + tab + \
            self.__feature    + tab + \
            str(self.__start) + tab + \
            str(self.__end)   + tab + \
            str(self.__score) + tab + \
            self.__strand     + tab + \
            self.__frame     + tab + \
            self.__attribute + tab + \
            self.__comment
        return string.strip(s)

    ##methods##

    #lenght of a feature = end-start+1
    def __len__(self):
        return self.end()-self.start()+1

    # add a offset to the features start and end coordinate
    def remap(self,offset):
        self.__start = self.start()+offset
        self.__end = self.end()+offset

    #return a copy of the feature
    def copy(self):
        return Feature(str(self))

    

    #checks for overlap between two features - if strand is 1 then the comparison is
    #strand sensitive.
    def overlap(self,other,strand = 0):
        size = 0
        if strand==0 or self.strand()==other.strand():            
            if other.start() <= self.end() <= other.end():
                size = self.end()-max(other.start(),self.start())+1
            elif self.start() <= other.end() <= self.end():
                size = other.end()-max(self.start(),other.start())+1
        return size

    # merges two features (i.e. coordinates 10-30, and 20-40 will result in 10-40)
    # if strand is 1 then the merging will only be done if the two features have the same strand
    def merge(self,other,strand = 0):
        if self.overlap(other,strand)>0:            
            newfeature = self.copy()
            newfeature.start(min(self.start(),other.start()))
            newfeature.end(max(self.end(),other.end()))
            return newfeature
        else: #no overlap - return reference to self
            return None

    
    # returns the intersection of two features (i.e. coordinates 10-30, and 20-40 will result in 20-30)
    # if strand is 1 then the intersection will only be done if the two features have the same strand
    # if the intersection is empty then none will be returned
    def intersect(self,other,strand = 0):
        newstart = 0
        newend = 0
        if self.overlap(other,strand)>0:            
            newfeature = self.copy()
            if other.start() <= self.end() <= other.end():
                newstart = max(other.start(),self.start())
                newend = self.end()
            elif self.start() <= other.end() <= self.end():
                newstart = max(other.start(),self.start())
                newend = other.end()
            newfeature.start(newstart)
            newfeature.end(newend)
            return newfeature
        else: #no overlap - return reference to self
            return None
      

