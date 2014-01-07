#! /usr/local/bin/python
#
#  

import os, sys, math
from gff import *


if __name__ == '__main__':
   
   print "GFF library test\n"

    
   if not len(sys.argv) in [2]:
      print "Usage: %s <GFF file>" % (os.path.basename(sys.argv[0]))
      sys.exit(1)

   filename = sys.argv[1]
   print "filename = %s " % filename 
   try:
        gff = GFF(filename)
   except GFF_IOError, s:
        print s
        sys.exit(0)

   print "Features in this file"
   print gff
   print "number of features in this file : %s \n" % (len(gff))
   truegenes = gff.filter(source = ["genbank"])
   predgenes = gff.filter(source = ["GenPred"])
   print "true exons from genbank:"
   print truegenes
   print "predicted genes:"
   print predgenes   
   (sensitivity,specificity,ac,cc,acp) = truegenes.calculate(predgenes,["exon"])


   (TN,TP,FP,FN) = truegenes.measures(predgenes,["exon"])
   print "Comparing predicted features against true features"
   print "TN = %s, TP = %s, FP = %s, FN = %s" % (TN,TP,FP,FN)
   print "Sn = %s, Sp = %s, ac = %s, cc = %s, acp = %s" % (sensitivity,specificity,ac,cc,acp)
   


     
  
 
