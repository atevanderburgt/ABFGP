"""
  BiRC GFF Python Library - version 1.0
"""
# Imports
import sys, string, operator, math, os

from feature import Feature, FeatureInputError, FeatureComment

class GFF_IOError(Exception):
    pass

class GFF:

   ## constructor
   def __init__(self, filename = ""):

      # Read file and generate list of features
      
      self.filename = filename
      if filename!="":
         try:
            if not os.path.isfile(filename):
               raise GFF_IOError("Could not open file '" + filename + "': Not a file.")
            f = open(filename, 'r')
         except IOError, s:
            raise GFF_IOError("Could not open file '" + filename+ "': " + str(s))
         lines = f.readlines()
         lines  = filter(string.strip, lines)
         
         features = []
         for l in lines:
             try:
                 f = Feature(l)
             except FeatureComment:
                 pass
             except FeatureInputError,s:
                 sys.stderr.write('Error in feature: ' +str(s)+"  " + l )                 
             else:               
                features.append(f)   
         self.__features = features
      else:
         self.__features = []

      

   ##methods##


   def remap(self,offset):
      for feature in self:
         feature.remap(offset)           

   def append(self,feature):
      self.__features.append(feature)

   #per default - strand is not considered
   def sort(self,func = None,strand = 0):
       if func==None:       
           self.__features.sort()
       else:
           self.__feature.sort(func)

   # returns the minimun start coordinate of all the features
   def minStart(self):
      if len(self)>0:
         min = self.__features[0].start()
         for feature in self:
            if feature.start()<min:
               min = feature.start()
         return min         
      else:
         raise IndexError

   # returns the minimun end coordinate of all the features
   def minEnd(self):
      if len(self)>0:
         min = self.__features[0].end()
         for feature in self:
            if feature.end()<min:
               min = feature.end()
         return min         
      else:
         raise IndexError

   # returns the maximum start coordinate of all the features
   def maxStart(self):
      if len(self)>0:
         max = self.__features[0].start()
         for feature in self:
            if feature.start()>max:
               max = feature.start()
         return max         
      else:
         raise IndexError

   # returns the maximum end coordinate of all the features
   def maxEnd(self):
      if len(self)>0:
         max = self.__features[0].end()
         for feature in self:
            if feature.end()>max:
               max = feature.end()
         return max         
      else:
         raise IndexError 

   # calculate the total length of the feature objects in the gff
   def featureLength(self):
       length = 0
       for feature in self:
           length = length + len(feature)
       return length
     

   #make copy of gff object
   def copy(self):
       newgff = GFF()
       for feature in self:
           newgff.__features.append(feature.copy())
       return newgff

   #remove feature from gff
   def remove(self,feature):
       self.__features.remove(feature)
       
   #concates a gff with another one
   def add(self,other):
       for feature in other:
           self.append(feature.copy())

   #overloading of '+' operator
   def __add__(self,other):
       newgff = self.copy()
       newgff.add(other)
       return newgff

    #if not entries[6] in ['+','-','.']:
   #returns a new gff object with the features that matches strand
   def strands(self,strand):
       newgff = GFF()
       for feature in self:
           if feature.strand()==strand:
               newgff.append(feature.copy())
       return newgff

   #returns a new gff object containing copies of features from self containing
   #feature names that matches one of the names in the namelist
   def featureNames(self,namelist = []):
       newgff = GFF()
       if namelist!=[]:
           for feature in self:
               if feature.feature() in namelist:
                   newgff.append(feature.copy())
       else:
           newgff = self.copy()
       return newgff

   ## returns a gff that contains no overlapping feature - i.e. overlapping features are merged.
   def mergeself(self,strand = 0):
       res = GFF()
       second = 1
       newgff = self.copy()
       newgff.sort()
       if len(newgff)<2:
           return newgff
       curFeature = newgff[0]
       while second<len(newgff):        
           f2 = newgff[second]
           if curFeature.overlap(f2,strand):  #overlap - so try next feature
               curFeature = curFeature.merge(f2,strand)                   
           else: #no overlap
               res.append(curFeature.copy())               
               curFeature = f2
           second = second + 1           
       return res

   
   def intersect(self,other,strand = 0):
       res = GFF()
       first = 0
       second = 0
       gff1 = self.copy()
       gff1.sort()
       gff2 = other.copy()
       gff2.sort()
       while first<len(gff1) and second<len(gff2):
           f1 = gff1[first]
           f2 = gff2[second]
           feature = f1.intersect(f2,strand)
           if feature!=None:
               res.append(feature.copy())
           if f1.end()<f2.end():
               first = first+1
           else:
               second = second+1
       return res
   

   def overlappingFeatures(self,strand = 0):
       overlaps = 0
       notdone = 1
       first = 0
       second = 1
       newgff = self.copy()
       newgff.sort()   
       if len(newgff)<2:
           return 0
       while notdone:        
           f1 = newgff[first]
           f2 = newgff[second]
           if f1.overlap(f2,strand):  #overlap - so try next feature              
               overlaps = overlaps + 1
               second = second + 1
           else:
               if strand==0:
                   first = first + 1
                   second = first + 1
               else:
                   second = second + 1
           if second>=len(newgff):  #check for end            
               if first>=len(newgff)-2:  # last feature
                   notdone = 0
               else:
                   first = first + 1
                   second = first + 1                    
       return overlaps   

   #save the gff to a file - overwrites.
   def save(self,filename,overwrite = 0):
       if (not os.path.isfile(filename)) or overwrite:
           oldsys = sys.__stdout__
           sys.stdout = open(filename,"w+")
           print self
           sys.stdout.close()
           sys.__stdout__ = oldsys
           return 1
       else:
           return 0

   def filter(self, seqname=[],source=[],features=[],strand=".",frame="."):
       newgff = self.copy()
       tempgff = GFF()
       
       if strand!=".":
           newgff = newgff.strands(strand)
       if features!=[]:
           newgff = newgff.featureNames(features)
       if seqname!=[]:
           for feature in newgff:
               if feature.seqname() in seqname:
                   tempgff.append(feature.copy())
           newgff = tempgff.copy()
       tempgff = GFF()
       if source!=[]:
           for feature in newgff:
               if feature.source() in source:
                   tempgff.append(feature.copy())
           newgff = tempgff.copy()
       tempgff = GFF()

       if frame!=".":
           for feature in newgff:
               if feature.frame()==frame:
                   tempgff.append(feature.copy())
           newgff = tempgff.copy()       
       return newgff
       
       
   def measures(self,gffPred,featurenames = [],overlaps = 0):
       (TN,TP,FP,FN) = (0,0,0,0)
       overlapcount = 0
       
       #calculate TN,TP,FP and FN

       if featurenames!=[]:
           gffT = self.featureNames(featurenames)
           gffP = gffPred.featureNames(featurenames)
       else:
           gffT = self.copy()
           gffP = gffPred.copy()
       gffT.sort()
       gffP.sort()
       #calculate start and end.
       if gffT.overlappingFeatures() or gffP.overlappingFeatures():
           print "In gff.measures : one of the gff objects contain overlapping features."
           sys.exit(0)
       if len(gffT)==0 or len(gffP)==0:
           return (TN,TP,FP,FN)
           
       minSeq = min(gffT.minStart(),gffP.minStart())       
       maxSeq = max(gffT.maxEnd(),gffP.maxEnd())
       #print minSeq
       #print maxSeq
                    
       tempgff = self.intersect(gffP)
       span = maxSeq-minSeq+1
 
       TP = tempgff.featureLength()
       FP = gffP.featureLength()-TP
       FN = gffT.featureLength()-TP   
       TN = span-(TP+FN+FP)+overlapcount
       return (TN,TP,FP,FN)
                  
   def calculate(self,gffPred,featurenames = [],overlaps = 0):
       (sn,sp,ac,cc,acp) = (0,0,0,0,0)
       tempgff = self.featureNames(featurenames)
       (TN,TP,FP,FN) = tempgff.measures(gffPred.featureNames(featurenames),featurenames,overlaps)
       div = math.sqrt( (TP+FN)*(TN+FP)*(TP+FP)*(TN+FN) )
       #correlations
       if div!=0:   
           cc  = (TP*TN-FN*FP)/div  
       else:
           cc = "undefined"  
       div1 = float(TP+FN)   
       div2 = float(TP+FP)      
       div3 = float(TN+FP)
       div4 = float(TN+FN)
       sum = 0
       count = 0  
       if div1!=0:  
           sum = sum + TP/div1
           count = count + 1  
       if div2!=0:  
           sum = sum + TP/div2
           count = count + 1  
       if div3!=0:  
           sum = sum + TN/div3
           count = count + 1  
       if div4!=0:  
           sum = sum + TN/div4
           count = count + 1
       if count!=0:
           count = float(count)
           acp = (1/count)*sum
       else:
           acp = 0.0

       ac  = (acp-0.5)*2
       
       if div1!=0:
           sn = TP/div1           
       else:
           sn = 0
           
       if div2!=0:
           sp = TP/div2
       else:
           sp = 0

       return (sn,sp,ac,cc,acp)

   def features(self):
      return self.__features

   ## overrides built-in functions ##   
   def __len__(self):
      return len(self.__features)	

   # used for gff[i]   
   def __getitem__(self,i):
      if i<len(self):
         return self.__features[i]
      else:  
         raise IndexError

   # used for gff[i] = new
   def __setitem__(self,i,new):
      if i<len(self):
         self.__features[i] = new
      else:
         raise IndexError
   
      
   def __str__(self):
      s = "";
      for feature in self:
         s = s + str(feature)+'\n'
      return s                 
      
##end of class gff ## 		
