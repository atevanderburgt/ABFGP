Version 1.00
Python GFF (General Feature Format) project.


Contents.
1. Overview.
2. Contens of this package.
3. License & contact information.

1. Overview

  This package is made for analyzing and do calculations on gff
  files. The GFF format is used for specifying the annotation of
  biological DNA sequences. This package contains two classes. 
  The 'Feature' class encapsulates a single line of a gff file. 
  The 'GFF' class is used to contain and analyze a complete gff file.
  To use the package just put the line 'from gff import *' in
  your file. See example.py for more details.
  See www.sanger.ac.uk/Software/formats/GFF/GFF_Spec.shtml for a
  specification of the GFF file format.  

2. Contents of this package.

   readme.txt		: this file.
   example.py           : example code for using the package.
   set.gff		: A example gff file. 
   doc/feature.html	: documentation for the feature class.	      
   doc/gff.html         : documentation for the gff class
   doc/birc-logo.gif    : graphics file.
   gff/feature.py	: the feature class
   gff/gff.py           : the gff class
   gff/__init__.py      : package initialization

   type "python example.py set.gff" to run the example script.

3. License & contact information.

  This software package is released under the GNU General Public
  License. 
  The copyright of this software still belongs to BiRC (Bioinformatic
  Research Center) at Aarhus University, Denmark. See www.birc.dk for
  more info. Author of this package is Martin Knudsen.
  Comments and suggestions are always welcome at martink@birc.dk




