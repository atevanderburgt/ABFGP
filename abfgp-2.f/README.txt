#####################################################################################
#### README ABFGP version 2.f                                                    ####
#####################################################################################

# content of this document

 1. Installation
 2. Required third party software
 3. Required Python modules
 4. Adjusting the settings/executables.py file
 5. Testing ABFGP
 6. Explaining the AbgpGeneLocusDirectory

 7. Running ABFGP (examples)

 8. Explanation of command line options for ABFGP (examples)
 9. How to use ABFGP in your research
10. The other ABFGP executables (examples)
11. Abbreviations used in ABFGP code and documentation


#####################################################################################
#### 1. Installation                                                             ####
#####################################################################################
####                                                                             ####
#### ABFGP is a collection of Python libraries that do the job.                  ####
#### No system-specific compilation or installation is required.                 ####
#### Just unpack somewhere in your code tree                                     ####
####                                                                             ####
#### ABFGP was developed and tested on UNIX (ubuntu, different versions).        ####
#### Sorry, but it will not work on Windows computers.                           ####
####                                                                             ####
#####################################################################################

# place the abfgp-VERSION-DATE.tar.gz archive in any desired directory
# unpack the archive
tar -zxf abfgp-2.0.tar.gz;
# cd to the abfgp directory
cd abfgp-2.0;


#####################################################################################
#### 2. Required third party software                                            ####
#####################################################################################
####                                                                             ####
#### ABFGP uses various third party software. Please download all, respecting    ####
#### the specific conditions and licensing terms. ABFGP was only tested for      ####
#### the software version's as listed below. No guarantee for proper functioning ####
#### for older or newer versions is given.                                       ####
#### We recommend to place all executables in the ./software/ directory.         ####
#### Some in-house developed executables are placed in this directory.           ####
####                                                                             ####
#### Files named PLEASE_INSTALL_xxxx are placed in the ./software directory,     ####
#### which exactly describe their content on our local servers; this could       ####
#### be helpful in case of failure of installed third party components.          ####
####                                                                             ####
#### All absolute paths of executables must be listed in:                        ####
#### settings/executables.py                                                     ####
####                                                                             ####
#### After installing (unpacking) the ABFGP package, please obtain               ####
#### all software !!in the specified versions!! and install it somewhere         ####
#### on you system. Most tools are expected to be present already because        ####
#### they represent mainstream DNA and PROTEIN sequence analyses tools.          ####
####                                                                             ####
#### Hint: because some versions might be different from the versions you        ####
#### are currently using, just install all in the ./software directory.          ####
#### Alternatively, if the software is already installed, you might just         ####
#### provide symlinks to these executables in the ./software directory.          ####
#### This will keep the required software for ABFGP clear, uncomplicated         ####
#### and manageable on a single place                                            ####
####                                                                             ####
#### After installation, please adjust the settings/executables.py files         ####
#### by specifying absolute paths to these software                              ####
####                                                                             ####
#####################################################################################

**REQUIRED**

blast-2.2.8
    blastall                    EXECUTABLE_BLASTALL     ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/
    formatdb                    EXECUTABLE_FORMATDB
clustalw-1.83                   EXECUTABLE_CLUSTALW
emboss                                                  easily obtainable from many mirrors; probably installed already
    getorf-6.2                  EXECUTABLE_GETORF
    tcode-6.2                   EXECUTABLE_TCODE
    transeq-6.2                 EXECUTABLE_TRANSEQ
hmmer-2.3.2                                             http://hmmer.janelia.org/software/archive
    hmmbuild                    EXECUTABLE_HMMBUILD
    hmmsearch                   EXECUTABLE_HMMSEARCH
ScanForMatches                  EXECUTABLE_SFM          http://blog.theseed.org/servers/2010/07/scan-for-matches.html

**OPTIONAL**

signalp-3.0                     EXECUTABLE_SIGNALP
tmhmm-2.0                       EXECUTABLE_TMHMM

**OPTIONAL** in the ABFGP source package

cexpander-1.0                       
    prep_launch.py              EXECUTABLE_CEXPANDER_ALLVSALL   ./software/cexpander-1.0/prep_launch.py
    cbalignp                    EXECUTABLE_CEXPANDER_CBALIGNP   ./software/cexpander-1.0/cbalignp_test_2/src/cbalignp
    cexpander_dr                EXECUTABLE_CEXPANDER_CEXPANDER  ./software/cexpander-1.0/cexpander_dr
fastalength.sh                  EXECUTABLE_FASTALENGTH          ./software/fastalength.sh
gfflength.sh                    EXECUTABLE_GFFLENGTH            ./software/gfflength.sh
unigeneannotation.py            EXECUTABLE_UNIGENEANNOTATION    ./unigeneannotation.py
load_gff.pl                     EXECUTABLE_LOAD_GFF             ./gff/load_gff.pl
get_gff_sequence_from_fasta.py  EXECUTABLE_GFF2FASTA            ./software/get_gff_sequence_from_fasta.py
                        
                        
# Some notes

-- EMBOSS getorf,tcode,transeq --
Newer and older versions of EMBOSS software are expected to have no issues (although not all versions were large-scale tested yet)
On various development machines, versions 6.2 and 6.4 were used.

-- hmmer --
hmmer-3.0 is expected to have no issues (although it was not large-scale tested yet)

-- blast --
New versions of blast could have issues.
Most likely issues relate to the BioPython NCBIStandalone blast parser,
that should be in sync with the version of blastall.

Blast+ is not supported; exact content of HSP objects changed drastically (strand,frame,etc).

-- clustalw --
Newer versions of clustalw are expected to have no issues (although not all were tested)

-- signalp --
signalp is not required but optional; if NOT installed, set required settings parameters to False
Newer versions are expected to have no issues (although not tested), but its predictions might be different.

-- tmhmm --
tmhmm is not required but optional; if NOT installed, set required settings parameters to False

-- ScanForMatches --
ScanForMatches can be downloaded at http://blog.theseed.org/servers/2010/07/scan-for-matches.html


#####################################################################################
#### 3. Required python modules                                                  ####
#####################################################################################
####                                                                             ####
#### We kept dependencies to other Python modules as limited as possible         ####
#### Required modules that might not be installed by default are:                #### 
####                                                                             ####
#### Bio.Blast.NCBIStandalone version 1.60 (BioPython)                           ####
#### Numeric verison 24.0 (Numpy)                                                ####
#### graph-0.50; (c) 2007 Pedro Matiello <pmatiello@gmail.com>)                  ####
####                                                                             ####
#### ABFGP was developed and tested using python versions 2.4, 2,6 and 2.7.4     ####
####                                                                             ####
#####################################################################################

Please solve by obtaining them and adjusting your PYTHONPATH to the respective
paths where they are installed. In case these modules are installed on a
position that python doesn't know of, update your PYTHONPATH system variable

export PYTHONPATH=$PYTHONPATH:/your/path/to/basedir/of/Biopython
export PYTHONPATH=$PYTHONPATH:/your/path/to/basedir/of/Numeric

The lightweight graph class available in graph-0.50
(Copyright (c) 2007 Pedro Matiello <pmatiello@gmail.com>) is subclassed for all
graph classes used in the ABFGP method. Is is packaged within the ABFGP distribution.


#####################################################################################
#### 4. Adjusting the settings/executables.py file                               ####
#####################################################################################
####                                                                             ####
#### After all the required software is installed, the absolute paths to         ####
#### these executables should be declared.                                       ####
#### Please open the settings/executables.py with your favourite text            ####
#### editor and fill in the respective absolute paths.                           ####
#### All concerned variable names are formatted like EXECUTABLE_XXXX,            ####
#### and are listed in section 2. of this README file                            ####
####                                                                             ####
#### Some global settings are required in the ./parsers and ./gene package.      ####
#### These packages are used by the ABFGP' developers in other projects too.     ####
#### Therefore, presence of these settings is most easily solved by providing    ####
#### absolute symlinks. You could consider a hard-copy too.                      ####
####                                                                             ####
#### Apart from specifying paths to executables, we strongly recommend           ####
#### **NOT TO CHANGE ANY** other settings. ABFGP was debugged and tested         ####
#### with a limited scope in terms of parameter range.                           ####
#### Most meaningful variables are adjustable on the command line                ####
#### in the main executable abfgp.py                                             ####
####                                                                             ####
#####################################################################################

# !IMPORTANT! set the absolute path MAIN_ABGP_PATH in settings/abgp.py
# !IMPORTANT! set the absolute path ABGP_OUTDIR_PATH in settings/abgp.py
# !IMPORTANT! provide absolute path PYTHON_PATH (python itself) in settings/executables.py
#             usually 'python' will suffice.
# !IMPORTANT! provide full paths to third party executables in settings/executables.py
#             required executables are described in the previous section.
#             hint: we just provided absolute symlinks in the software/ directory
#             of the required third party software.

# cd to the MAIN_ABGP_PATH directory

# cleanup all putatively present files / symlinks  to the executable.py file(s)
# in destination directories
rm -f parsers/executables.py;
rm -f parsers/executables.pyc;
rm -rf gene/settings;

# set current machine's executables file to
# settings/executables.py and parsers/executables.py
# by providing symlinks in directories that require the settings
ln -s $(pwd)/settings/executables.py $(pwd)/parsers/executables.py;
ln -s $(pwd)/settings $(pwd)/gene/settings;

#####################################################################################
#### 5. Testing ABFGP                                                            ####
#####################################################################################
####                                                                             ####
#### Once steps 1-4 are performed, ABFGP should be ready for execution.          ####
#### However, we recommend to test first if indeed all requirements are met.     ####
#### Two testing scripts are provided for this purpose:                          ####
####                                                                             ####
####    # a listing of ALL imports done in any script or module                  ####
####    abfgp-test-imports.py                                                    ####
####                                                                             ####
####    # a checkup on required third party software                             ####
####    abfgp-test-executables.py                                                ####
####                                                                             ####
#####################################################################################

# cd to the abfgp directory
# adjust your PYTHONPATH to include ABFGP's main directory
# and run the execuatble file abfgp-test-imports.py
export PYTHONPATH=$PYTHONPATH:/your/path/to/abfgp-2.0; python abfgp-test-imports.py;

# This should print "All imports succesfully performed!" on stdout.
# If any exception is raised on stdout/stderr, this means required python modules
# are not where they are expected. Please read and reperform REAME section 3.

# cd to the abfgp directory
export PYTHONPATH=$PYTHONPATH:/your/path/to/abfgp-2.0; python abfgp-test-executables.py;

# This should print 18 lines which all start with ' True '
# In case an executable can't be located, isn't executable or might have 
# another flaw, it is reported here. Please take action prior to ABFGP usage.


#####################################################################################
#### 6. Explaining the AbgpGeneLocusDirectory                                    ####
#####################################################################################
####                                                                             ####
#### The AbgpGeneLocusDirectory is the type of input data that was used for      ####
#### whole-genome scale re-annotation of gene models in the manuscript.          ####
#### Reason not to use the most simple plain fasta DNA of target and informant   ####
#### loci is that additional data had to be provided for some basic tasks:       ####
#### - the annotated gene model(s) must be known to check for revisions          ####
#### - the annotated gene model(s) could/does provide useful information         ####
#### - (aligned) unigene data must be linked to (informant) gene loci            ####
####                                                                             ####
#### This (and other things) is all taken care of in the AbgpGeneLocusDirectory. ####
#### It is a simple folder somewhere on your file system that contains           ####
#### fasta and GFF files that contain this information. The class                ####
#### AbgpGeneLocusDirectory ( in abgpgenelocusdirectory.py ) can read the        ####
#### information in these folders and makes them accessible to the ABFGP method  ####
####                                                                             ####
#####################################################################################

A folder is recognized as a AbgpGeneLocusDirectory as soon as it contains:
- a "*.dna*.fa"     (regex!) fasta file that contains the genomic DNA locus
- a "*.locus.gff"   (regex!) gff file describing its absolute position
Optional, useful files are
- a "*.gene*gff"    (regex!) gff file describing current annotation
- a "*.unigene*gff" (regex!) gff file(s) describing aligned unigenes


All gff files must refer to *absolute* coordinates; ABFGP does the math to recalculate
the relative locus on the provided slice of genomic DNA of the informant locus.
It is best shown by example. In the testdatafolder serveral full examples
of orthologous groups of AbgpGeneLocusDirectories are provided. This is the
content of testdata/exampleAbgpLocusDirs/CFU_827306/DOTSE_072982:

    -rw-r--r--  1 avdb users 5225 Oct 17 09:46 DOTSE_072982.dna.1500.1500.fa
    -rw-r--r--  1 avdb users  103 Oct 17 09:46 DOTSE_072982.locus.gff
    -rw-r--r--  1 avdb users  971 Oct 17 09:46 DOTSE_072982.gene.gff
    -rw-r--r--  1 avdb users  369 Oct 17 09:46 DOTSE_072982.unigene.isotig06737.full.gff
    # these files are for decorative purposes only
    -rw-r--r--  1 avdb users  640 Oct 17 09:46 DOTSE_072982.locus.signalp.out
    -rw-r--r--  1 avdb users  102 Oct 17 09:46 DOTSE_072982.locus.tmhmm.out
    -rw-r--r--  1 avdb users  118 Oct 17 09:46 DOTSE_072982.signalp.out
    -rw-r--r--  1 avdb users    0 Oct 17 09:46 DOTSE_072982.tmhmm.out
    # for completeness, the protein that corresponds to the predicted gene model is provided
    -rw-r--r--  1 avdb users  683 Oct 17 09:46 DOTSE_072982.protein.fa

This is the content of the GFF files:
# testdata/exampleAbgpLocusDirs/CFU_827306/DOTSE_072982/DOTSE_072982.locus.gff
    scaffold_7      abgp_locus      abgp_locus      263918  269089  .       -       .       gene_id "DOTSE_072982"; abgp_locus "DOTSE_072982"
# testdata/exampleAbgpLocusDirs/CFU_827306/DOTSE_072982/DOTSE_072982.gene.gff
    scaffold_7      JGI     exon    265418  265659  .       -       .       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"; transcriptId 73024
    scaffold_7      JGI     CDS     265542  265659  .       -       0       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"; proteinId 72982; exonNumber 3
    scaffold_7      JGI     stop_codon      265542  265544  .       -       0       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"
    scaffold_7      JGI     exon    265726  267286  .       -       .       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"; transcriptId 73024
    scaffold_7      JGI     CDS     265726  267286  .       -       1       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"; proteinId 72982; exonNumber 2
    scaffold_7      JGI     exon    267339  267589  .       -       .       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"; transcriptId 73024
    scaffold_7      JGI     CDS     267339  267525  .       -       2       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"; proteinId 72982; exonNumber 1
    scaffold_7      JGI     start_codon     267523  267525  .       -       0       gene_id "DOTSE_072982"; name "estExt_fgenesh1_kg.C_7_t10080"
# testdata/exampleAbgpLocusDirs/CFU_827306/DOTSE_072982/DOTSE_072982.unigene.isotig06737.full.gff
    scaffold_7      GenomeThreader  UGexon  265418  265659  .       -       .       UniGene isotig06737
    scaffold_7      GenomeThreader  UGexon  265726  267286  .       -       .       UniGene isotig06737
    scaffold_7      GenomeThreader  UGexon  267339  267589  .       -       .       UniGene isotig06737
    scaffold_7      GenomeThreader  UGintron        265660  265725  .       -       .       UniGene isotig06737
    scaffold_7      GenomeThreader  UGintron        267287  267338  .       -       .       UniGene isotig06737

As you can see, all coordinates refer to the absolute locus of the in informant gene locus.
Provide only like this; ABFGP does the math to recalculate towards the provided slice of genomic DNA.
For unigenes and genes, it is not required to provide GFF for its introns.
NOTE: Use the literal GFF gclass 'gene_id' for gene models (or adjust GFF_GENE_GCLASS settings/gff/currentannotation.py).
NOTE: Use the literal GFF fmethod 'CDS' and 'exon' for exons of genes (or adjust in settings/gff/currentannotation.py).
NOTE: Use the literal GFF gclass 'UniGene' for unigene exons (or adjust UniGene in settings/gff/unigene.py).
NOTE: Use the literal GFF fmethod 'UGexon' for unigene exons (or adjust GFF_UGEXON_FMETHOD in settings/gff/unigene.py).
Internally, ABFGP will inspect the aligned unigene and annotate it structurally (coding exons, UTRs, fragment etc.).

In case available, we recommend you to provide unigene data.
This can be easily done by the following approach:
- make AbgpGeneLocusDirectories for all your annotated gene models
- map unigene data to you genome with your favourite tool (GeneSeeker, PASA, etc).
- use bedtools (clusterbed) to cluster loci with annotated genes and loci that have expression.
- abstract the corresponding GFF lines from the aligned unigene GFF file and place them in the
  respective AbgpGeneLocusDirectory.

As you can see, the example of the provided gene model is as it could have been obtained
directly from the GFF/GTF gene catalogue file. In the dataset of fungal genomes and annotations
that was used for the manuscript, it proved harder as expected to make single, uniform script that worked for
all annotations. Be ready to encounter these and many other conflicts:

- in some repositories, gene and transcript ID are the same, but mostly not
- in some repositories, FREF (first column) did not correspond to fasta
  entries in the corresponding genome files. Some accessions seemed to be
  named: 'supercont3.1 of Gibberella moniliformis',
  'supercontig_2.1 of Fusarium oxysporum f. sp. lycopersici', where even the
  usage of underscores conflicted beteen fasta and GFF files.
- in some repositories, 'gene' and/or 'transcript' containers are not provided
- in the JGI's format, alternatively spliced loci can have completely
  unrelated names.
  

Here, we provide a small piece of code that should be able to help you automating this step. 
Below example works like charm on the gff-3 annotation of Saccharomyces cerevisae,
a species that we did not use in our study. This to show that it is trivially easy
to obtain a genome-scale catalogue of AbgpGeneLocusDirectories

-------------------------------------------------------------------------------------------------------
--- Download & prepare required data
-------------------------------------------------------------------------------------------------------

# download bedtools; provide here your path to the bedtools suite
# we are a great fan of bedtools!
exe_bedtools=/home/avdb/software/bedtools-2.17.0/bin/bedtools

# specify file names
fname_gff=Schizosaccharomyces_pombe.ASM294v2.20.gff3
fname_genome=genome.spombe.mfa
fname_bedidx=genome.spombe.bed

# download genome and gff3 annotation of S.pombe
wget ftp://ftp.ensemblgenomes.org/pub/fungi/current/gff3/schizosaccharomyces_pombe/Schizosaccharomyces_pombe.ASM294v2.20.gff3.gz
unzip Schizosaccharomyces_pombe.ASM294v2.20.gff3.gz
for chr in I II III MT MTR AB325691
do
    wget ftp://ftp.ensemblgenomes.org/pub/fungi/current/fasta/schizosaccharomyces_pombe/dna/Schizosaccharomyces_pombe.ASM294v2.20.dna.chromosome.$chr.fa.gz
    gunzip Schizosaccharomyces_pombe.ASM294v2.20.dna.chromosome.$chr.fa.gz
done
cat Schizosaccharomyces_pombe.ASM294v2.20.dna.chromosome.*.fa > $fname_genome
# make genome.spombe.bed index for bedtools
cat $fname_genome | awk '/^exit$/ {exit}{print}' | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' | sed 1d | awk '{ if (substr($1,0,1)==">") { printf $1"\t" } else { printf length($1)"\n" } }' | tr -d ">" > $fname_bedidx;
# cleanup single chromosome files
rm Schizosaccharomyces_pombe.ASM294v2.20.dna.chromosome.*.fa;

-------------------------------------------------------------------------------------------------------
--- Specify litteral GFF/GTF strings for search/replace commands
-------------------------------------------------------------------------------------------------------

# specify fmethod names & settings
organism_tag=spombe                 # usefull to recognize the source of your informant loci 
fmethod_gene="protein_coding_gene"  # literal feature name in GFF's 3th column
fmethod_transcript="transcript"     # literal feature name in GFF's 3th column
fmethod_abfgpgene="gene"            # not required in ABFGP, but provide it for sake of completeness
fmethod_abfgplocus="abgp_locus"     # feature name to translate to
gclass_gene="gene_id"               # FF_GENE_GCLASS in settings/gff/currentannotation.py
gene_locus_nt_flank=1000            # extent gene loci with extra flanking sequence.
                                    # this is helpfull/required in case this
                                    # gene model is actually incorrect!

# ABFGP works with GFF2. We advise to convert GFF3 to GFF2 first.
# Here, for simplicity, we strip all the extra metadata and rewrite
# """ID=SPAC1556.08c;biotype=prot....""" to """ID SPAC1556.08c"""
mimic_gff3_to_gff2="sed 's/^\([^=]\+\)=\([^;]\+\).*/\1 \2;/g'"

-------------------------------------------------------------------------------------------------------
--- Here the actual functionality happens; loop over each gene locus and make an AbgpGeneLocusDirectory
-------------------------------------------------------------------------------------------------------

# okay, now loop over all the gene encoding loci in this annotation file
for genelocusid in `cat $fname_gff | awk -F'\t' '{ if ($3=="'$fmethod_gene'") { print $9 } }' | sed 's/\([^;]\+\);.*/\1/' | sed 's/[^=]\+=//'`
do
    # Pombase uses an uniform way of labeling transcripts.
    # The first/single transcript of a gene has the suffix '.1',
    # and so on for alternative splicing. ABFGP has for sake of
    # simplicity the same ID for gene and transcript. This is because
    # little to no alternative splicing was annotated in the respective
    # gene catalogues. Therefor, we have to get rid of the '.1' suffix
    # of all the exon and CDS tracks in the corresponding GFF.
    # This can be easily done by just overwriting the gclass/gname
    # pair with $gclass_gene and $genelocusid.
    enforce_gclass_gname="sed 's/^\(.*\t\)\([^ ]\+ [^ ]\+\)$/\1'$gclass_gene' '$genelocusid'/'"

    gffgene=$(grep "$genelocusid" $fname_gff | awk -F'\t' '{ if ($3=="'$fmethod_gene'") { print $0 } }' | eval $mimic_gff3_to_gff2 | eval $enforce_gclass_gname);
    linecnt=$(printf "$gffgene\n" | wc -l);
    if [ "$linecnt" == "1" ];
    then
        # unique genelocusid. This will become an AbgpGeneLocusDirectory
        mkdir $genelocusid;
        # extend genomic slice with bedtools; translate FMETHOD for recognition and convert to simple GFF2 format
        printf "$gffgene\n" | $exe_bedtools slop -b $gene_locus_nt_flank -i stdin -g $fname_bedidx | sed 's/'$fmethod_gene'/'$fmethod_abfgplocus'/' > $genelocusid/$genelocusid.locus.gff;
        # get genomic slice of DNA that corresponds to this locus
        cat $genelocusid/$genelocusid.locus.gff | $exe_bedtools getfasta -s -fi $fname_genome -bed stdin -fo stdout | sed 's/>/>'$organism_tag'_'$genelocusid' /' > $genelocusid/$genelocusid.dna.fa;
        gffblock=$(grep "$genelocusid" $fname_gff | eval $mimic_gff3_to_gff2 | eval $enforce_gclass_gname);
        trnscnt=$( printf "$gffblock" | awk -F'\t' '{ if ($3=="'$fmethod_transcript'") { print $0 } }' | wc -l )
        if [ "$trnscnt" -eq "1" ];
        then
            # A 1:1 relation of gene to transcript.
            printf "$gffblock\n" | sed 's/'$fmethod_gene'/'$fmethod_abfgpgene'/' > $genelocusid/$genelocusid.gene.gff;
            pi=3
        elif [ "$trnscnt" -ge "2" ];
        then
            # In case of AS, we choose here for simplicity to just store a single gene model
            # in the AbgpGeneLocusDirectory. Here, we give the example of how
            # to make the distinction on the distinct transcript level.
            # This will vary according to the type and origin of the annotation
            # you are working with. For fungi, alternative splicing is not occuring abundantly.
            for transcriptid in `grep $genelocusid $fname_gff | awk -F'\t' '{ if ($3=="'$fmethod_transcript'") { print $9 } }' | eval $mimic_gff3_to_gff2 | awk '{ print $2 }' | tr -d ";"`
            do
                echo "Warning:altsplicing :" $genelocusid "->" $transcriptid;
                enforce_gclass_transcriptname="sed 's/^\(.*\t\)\([^ ]\+ [^ ]\+\)$/\1'$gclass_gene' '$transcriptid'/'"
                gfftrblock=$(grep "$transcriptid" $fname_gff | eval $mimic_gff3_to_gff2 | eval $enforce_gclass_transcriptname);
                printf "$gfftrblock\n" | sed 's/'$fmethod_gene'/'$fmethod_abfgpgene'/' > $genelocusid/$genelocusid.gene.$transcriptid.gff;
            done
        else
            echo "Warning:notranscript: $genelocusid"
            grep "$genelocusid" $fname_gff | awk -F'\t' '{ if ($3=="'$fmethod_gene'") { print "Warning:nonunique: "$0 } }' 
        fi
        # done! gene locus directory is created for this gene!
        echo $genelocusid
        #ls -al $genelocusid;
        #cat $genelocusid/$genelocusid.*.gff;
    else
        # problem! $genelocusid is referring to >1 genes. Extend sed regular expression
        echo "Warning:nonunique   : $genelocusid"
        grep "$genelocusid" $fname_gff | awk -F'\t' '{ if ($3=="'$fmethod_gene'") { print "Warning:nonunique: "$0 } }' 
    fi
done > log.txt;

At the end, inspect the log.txt file. There is actually only 1 alternatively spliced gene annotated.
A very small numer of 3 gene loci failed because their ID's weren't uniquely occuring in this GFF3 file.
We stop this over-detailed example here; you will be able to tune this script to your own needs.


#####################################################################################
#### 7. Running ABFGP (examples)                                                 ####
#####################################################################################
####                                                                             ####
#### The following section contains some typical command lines to ABFGP          ####
#### Do not forget to: export PYTHONPATH=$PYTHONPATH:/your/path/to/abfgp-2.0;    ####
####                                                                             ####
#### For a full explanation on the command line options, including the           ####
#### various types of input data, please read the following section 8.           ####
####                                                                             ####
#####################################################################################

export PYTHONPATH=$PYTHONPATH:/usr/share/pyshared
--------------------------------------------------------------------------------------
python abfgp.py --multifasta testdata/exampleMultifasta/CFU_840409.andinformants.mfa --target cfu
python abfgp.py --dirwithloci testdata/exampleAbgpLocusDirs/CFU_840409 --target CFU_840409
python abfgp.py --filewithloci testdata/exampleFilewithloci/CFU_840409.relativepaths.bdbh.csv 
--------------------------------------------------------------------------------------

    The above three command lines are roughly similar.
    Biggest and very important and impacting difference, is that in the
    first (--multifasta) style, only the target and informant gene loci are offert.
    No prior knowledge on annotated gene models and/or aligned unigenes (as GFF)
    are offered. In the second two command lines (--dirwithloci and --filewithloci),
    AbgpGeneLocusDirectories are offered as input (see section 6. of this README).
    For proper use of ABFGP, we strongly recommend to use AbgpGeneLocusDirectories.
    In section 9. of this README, a basic example is worked out how to avoid
    AbgpGeneLocusDirectories but still use (non-aligned!) unigenes as evidence.
    
    A note on --target cfu --target CFU_840409 or no --target applied:
    In our whole-genome scale re-annotation implementation of the ABFGP method
    as described in the manuscript, we decoupled the proces of making
    AbgpGeneLocusDirectories and mining informants from an all-vs-all
    similarity matrix from the actual running of the ABFGP method.
    ABFGP runs were prepared as a `filewithloci` that point towards
    AbgpGeneLocusDirectories somewhere on our filesystem. Please read
    section 9. of this README for why and how. The target locus was
    always the first locus in the --filewithloci. So, no --target has
    to be applied. Please feel free to assign any different informant
    as target sequence! In all other input data types, --target has to 
    be applied. The difference between 'CFU_840409'  and 'cfu' is because
    of some (non-generic) implementation (protein)-header-2-organism mapping
    we applied. In settings/dbwarehouse.py resides the variable:
    ABGP_LOCUS_IDENTIFIER2ORGANISM_MAPPING
    Because we worked with a confined set of species, and we only
    used putative orthologous informants (no paralogs of the same species!),
    we had a 1:1 relationship between protein and gene / organism identifiers.
    For simplicity, we made those available as alias.
    
    Output data and all intermediary files are created in --outdir (or
    a unique folder name in the directory specified in
    settings/abgp.py: ABGP_OUTDIR_PATH
    
    In the output GFF file, the re-annotated gene model is described fairly at the bottom.
    The GFF file is highly decorated with information that could assist in the manual
    curation of the gene model (see Figure 3 in the manuscript).
    
    The re-annotated gene model's gclass is 'AbfgpGeneStructure'
    Exons and introns fmethods are named 'AbfgpExon', 'AbfgpIntron'
    The fscore as provided in the AbfgpExon and AbfgpIntron tracks
    is the `ok` (1) or `doubtful` (0) label that is assigned by the introspection module
    as described in the manuscript.
    
    To visualize the output, configure the Generic Genome Browser (http://gmod.org/wiki/GBrowse)
    of any other framework that can visualize GFF. In case you choose for the 
    Generic Genome Browser, our config file is provided in settings/GenericGenomeBrowser.abfgp.conf
    Track names in this file correspond to those named in setting/gff
    
    Various other orthologous groups are provided in testdata/.
    They comprise a variable number of informant gene loci and are of
    different input type (--filewithloci, --dirwithloci, --multifasta).
    The 7 samples provided in testdata/exampleMultifasta, testdata/exampleFilewithloci
    and testdata/exampleAbgpLocusDirs are the same samples provided in different format.

#####################################################################################
#### 8. Explanation of command line options for ABFGP (examples)                 ####
#####################################################################################
####                                                                             ####
#### The following section explains the command line options in detail.          ####
####                                                                             ####
#####################################################################################

##
## General options
##

--------------------------------------------------------------------------------------
-h, --help          show this help message and exit
--------------------------------------------------------------------------------------

--------------------------------------------------------------------------------------
-v, --verbose       print as much status messages to stdout as possible
-q, --quiet         print only mayor milestone messages to stdout
-s, --silent        print not a single status messages to stdout,
                    just the (final) outcome
--------------------------------------------------------------------------------------

    ABFGP produces a hell of information on STDOUT which mainly helped during
    development. Most of it will be highly cryptic to non-developers.
    It logs all intermediate results; it is possible to exactly trace
    in which step(s) of the execution which part(s) of gene models
    are detected or discarded.
    
    If you are disturbed by this, please use --quiet, or even --silent

--------------------------------------------------------------------------------------
-o OUTDIR, --outdir=OUTDIR
                    specify outdir where to create result file(s); auto-
                    generated when not specified
--------------------------------------------------------------------------------------

    The ABFGP method creates by default an output folder
    that stored results (and many temporarily files).
    It is automatically and uniquely time-stamped like:
    /tmp/abgpD20131018T111233R518
    where /tmp can be specified in ABGP_OUTDIR_PATH in settings/abgp.py.
    You can specify any desired different output folder here.

--------------------------------------------------------------------------------------
-f, --force         force file-overwrite of existing files and creation of
                    new directories
--------------------------------------------------------------------------------------

    In case you specified with --outdir your own (already existing) output folder,
    and you perform distinct runs of ABFGP (e.g. by changing informants),
    output files might get overwritten. This argument enforces this.
    Otherwise, unique suffixes will be added in order to prevent existing ABFGP
    output files to be overwritten.
                 
##
## Major options defining behaviour of the ABFGP method
##

--------------------------------------------------------------------------------------
--target=TARGET     use this Organism/Gene identifier asa target (and all
                    others as informants)
--------------------------------------------------------------------------------------

    For all input options (except --filewithloci), the option --target must be applied
    to assign which gene / species sequence must be used as target. If you don't provide
    it, ABFGP will prompt you a list from which to choose.

--------------------------------------------------------------------------------------
--abinitio          do not use the TARGET annotated gene structure as
                    prior knowledge
--------------------------------------------------------------------------------------

    By default, and in case AbgpGeneLocusDirectories are offered as input that include
    annotated gene models, the ABFGP method does NOT work on a strict ab initio bases.
    In the non-abinitio case, the annotated gene models are used as prior knowledge.
    - Similarity searches are prioritized for those ORFs that encode exons.
      This causes speedup in the ABFGP method
    - In case of presence of small exons in gene models, the ABFGP method is
      warned that (parts of) the all-vs-all gene model alignment (the GeneStructureGraph)
      likely has flaws at the start. Again, ORFs that contain annotated exons are prioritized in
      the search for missing fragments of gene structures.
    - This is especially true for the SmallAnnotatedFirstExonWarning (see abgp_warnings.py);
      many studies have indepentantly shown that eukaryotic genes, including yeast and fungi,
      are (slightly) enriched for introns in the first part of the CDS of genes. Small first
      exons are difficult to detect based on a similarity search (simply because they are small).
      In case the SmallAnnotatedFirstExonWarning is encountered, the search for these is
      prioritized to the already annotated start sites / exons. This causes a huge speedup.
    We recommend NOT to use the --abinitio flag.
    
##
## Options that specify input data
##

#############################################################################################
ABFGP can basically work on two different input data types:
    - fasta DNA sequences
    - AbgpGeneLocusDirectory (AbgpGeneLocus or AbgpLocus in short);  -> see Section 6. in the README

There are 5 different options that facilitate providing this input data.
An interesting possibility is to mix various input data, which is perfectly allowed.
By example, provide a basic set of informants in --multifasta / --filewithloci / --dirwithloci,
but additionally provide a few separate fasta files(s) that could be putative
useful informants too.

#############################################################################################

--------------------------------------------------------------------------------------
--loci              perform gene prediction on these AbgpGeneLocus
                    directorie(s): --loci <dirA> .. <dirX>
--------------------------------------------------------------------------------------

    Note: if --filewithloci in not applied, --target should be specified!
    
    Provide a (as long as desired) list of informant AbgpGeneLocusDirectories.
    In case AbgpGeneLocusDirectories are already grouped based on homology in a
    separate folder, it is more convenient to use the --dirwithloci options

    The command line:
    python abfgp.py --dirwitloci testdata/exampleAbgpLocusDirs/CFU_827306/ANID_07423 .... etc 25 more AbgpGeneLocusDirs etc .... testdata/exampleAbgpLocusDirs/CFU_827306/VDAG_08650

    Is equivalent to 
    python abfgp.py --dirwitloci testdata/exampleAbgpLocusDirs/CFU_827306
    
    In the provided testdata, hundreds of AbgpGeneLocusDirectories are available
    that are grouped in homologous groups in directories to be used as a --dirwithloci

--------------------------------------------------------------------------------------                    
--dirwithloci=DIRWITHLOCI
                    perform gene prediction on this (preselected)
                    LocusDirectory with AbgpGeneLocusDirectories
--------------------------------------------------------------------------------------

    Note: if --filewithloci in not applied, --target should be specified!

    Provide a folder name that includes *ONLY* homologous/orthologous AbgpGeneLocusDirectories.
    It is perfectly okay if non-AbgpGeneLocusDirectories are in this folder, 
    the *ONLY* refers to that there are no non-homologous AbgpGeneLocusDirectories listed.
    A typical ABFGP command line is:

    python abfgp.py --dirwitloci testdata/exampleAbgpLocusDirs/CFU_827306 --target CFU_827306
    
    As target, any of the AbgpGeneLocusDirectories in testdata/exampleAbgpLocusDirs/CFU_827306
    could be assigned.

--------------------------------------------------------------------------------------
--multifasta=MULTIFASTA
                    perform gene prediction on multifasta DNA file:
                    --multifasta <fileA>[length range 2kb..20kb]
--------------------------------------------------------------------------------------

    Note: if --filewithloci in not applied, --target should be specified!

    Provide a fasta DNA file with multiple entries that correspond to homologous/orthologous
    gene-encoding loci, and specify a single entry as --target.
    The command line:
    
    python abfgp.py --multifasta testdata/exampleMultifasta/CFU_827306.andinformants.mfa --target CFU_827306
    
    Will load the same informants as:

   python abfgp.py --dirwitloci testdata/exampleAbgpLocusDirs/CFU_827306 --target CFU_827306
    
    
--------------------------------------------------------------------------------------                    
--dna               perform gene prediction on dna file(s): --dna <fileA>
                    ... <fileX> [range 2kb..20kb]
--------------------------------------------------------------------------------------

    Note: if --filewithloci in not applied, --target should be specified!
    
    Provide a (as long as desired) list of informant AbgpGeneLocusDirectories.
    gene-encoding loci, and specify a single entry as --target.
    Interesting option in case you want to provide a single / few informant DNA sequences
    from a different source as where you obtained the bulk of your input data from.

--------------------------------------------------------------------------------------
--filewithloci=FILEWITHLOCI
                    perform gene prediction on this file with absolute
                    paths to AbgpGeneLocus directories
--------------------------------------------------------------------------------------

    Note: if --filewithloci in applied the FIRST occuring AbgpGeneLocus is taken as --target.
    Note: this can be simply adjusted by specifying a different --target
    
    A `filewithloci` is a (txt) file that can have any content, but should contain paths
    to AbgpGeneLocusDirectories. The file is scanned for (existing) folders, and
    those are checked for being an AbgpGeneLocusDirectory. The first occurring line
    in the `filewithloci` that passes this test, is assumed to be the --target locus.
    
    Example `filewithloci` are provided in testdata/exampleFilewithloci/
    Seven distinct files are provided in two fashions:
    
    - containing *absolute paths* to the local file system of the developer
    - containing *relative paths* to testdata/exampleAbgpLocusDirs/
    
    We strongly recommend to use absolute paths; we provided the relative paths 
    for demonstration purposes. You can use the following command line (make
    shure you are in ABFGP's main folder because of the relative paths!):
    
    python abfgp.py --dirwitloci testdata/exampleFilewithloci/CFU_827306.relativepaths.bdbh.csv
    
    This command line is equivalent to:
    
    python abfgp.py --dirwitloci testdata/exampleAbgpLocusDirs/CFU_827306 --target CFU_827306

    Please have a quick look at this file: (testdata/exampleFilewithloci/CFU_827306.relativepaths.bdbh.csv)

    #############################################################################################
    # identifier               : CFU_827306
    # genomedirs_to_ignore     : []
    # genomedirs_to_use        : []
    # search_mode              : BDBH
    # allow_paralogs           : False
    # safeorthologs_ratio      : 0.9
    # maximal_length_ratio     : 1.5
    # maximal_num_loci         : 30
    # minimal_bitscore_ratio   : 0.1
    # minimal_num_loci         : 5
    # minimal_overlap_ratio    : 0.7
    # similarity data
    cfu     dotse   CFU_827306      DOTSE_072982    3084    1.0     1.0     0.9418  0.9175
    cfu     crypa   CFU_827306      CRYPA_283070    1798    0.9871  0.9983  0.549   0.5462
    cfu     necha   CFU_827306      NECHA_069505    1745    0.979   0.9885  0.533   0.5256
    cfu     triat   CFU_827306      TRIAT_139326    1721    0.9839  0.9983  0.5255  0.5284
    cfu     mycfi   CFU_827306      MYCFI_29007     998     0.9984  1.0     0.3047  0.3039
    ...     .....
    ...     19 lines omitted here
    ...     .....
    cfu     pgt     CFU_827306      PGTT_08418      479     0.9742  0.9919  0.1463  0.1501
    cfu     cawg    CFU_827306      CAWG_03818      368     0.9694  0.992   0.1124  0.2937
    # target locus
    testdata/exampleAbgpLocusDirs/CFU_827306/CFU_827306
    # informant loci
    testdata/exampleAbgpLocusDirs/CFU_827306/DOTSE_072982
    testdata/exampleAbgpLocusDirs/CFU_827306/CRYPA_283070
    testdata/exampleAbgpLocusDirs/CFU_827306/NECHA_069505
    testdata/exampleAbgpLocusDirs/CFU_827306/TRIAT_139326
    testdata/exampleAbgpLocusDirs/CFU_827306/MYCFI_29007
    ...
    ... 18 lines omitted here
    ...
    testdata/exampleAbgpLocusDirs/CFU_827306/PGTG_08418
    testdata/exampleAbgpLocusDirs/CFU_827306/CAWG_03818
    #############################################################################################
    
    
    This file (more precise, CFU_827306.bdbh.csv in absolute path version) is exactly
    the file that was used to obtain the results as presented in the manuscript.
    It shows how the the developers implemented their approach to genome-wide re-annotation:

    - In the middle: Entries obtained from the 29-species all-vs-all protein similarity matrix (9-column width)
                     target species id (protein_id prefixes were used as species tags)
                     informant species id
                     target protein id
                     informant protein id
                     total bits
                     %length coverage of target protein on informant protein
                     %length coverage of informant protein on target protein
                     % of maximum attainable bits (target protein aligned to itself)
                     % of maximum attainable bits (informant protein aligned to itself)
    - On top:        Settings used to obtain homologs/orthologs from this matrix using CFU_827306 as --target                        
    - At the bottom: the (absolute/relative) paths to homologous/orthologous AbgpGeneLocusDirectories
                     obtained from different species.
                     
    In our own implementation, we organized informant gene loci were per genome (not per orthologous group),
    and their AbgpGeneLocusDirectories were named by their protein id. As simple query in the above
    shown snapshot of the similarity matrix (``get the best hit for each species with at least x% and y% coverage,
    z% bits etc etc...``) yielded a set of protein identifiers that could via their organism tags
    directly be mapped to AbgpGeneLocusDirectories; /out/path/to/<SPECIES_ID>/loci/<PROTEIN_ID>.
    Just an idea how to manage your project; it worked fine for us!
    New informant species can be added incrementally; the only required step is to
    describe a new informant's gene catalogue as AbgpGeneLocusDirectories (see section 6.)
    and update your own implementation of a similarity matrix with this new informant species.
    
                     
                     
##
## Options that help managing the (number of) provided informants
##

--------------------------------------------------------------------------------------
--minimal_num_loci=MINIMAL_NUM_LOCI, --min_num_loci=MINIMAL_NUM_LOCI
                    minimal number of gene loci to use as input [default:
                    5, range 2..50]
--------------------------------------------------------------------------------------
                        
    Data presented in the manuscript was obtained with default settings
    of in total five loci: at least four informants for one target locus.
    Lowering this threshold will reduce prediction accuracy, and is therefore
    not recommended. Single exception to this rule is the case of
    available unigene data for some of the sparse informant(s).

--------------------------------------------------------------------------------------
--maximal_num_loci=MAXIMAL_NUM_LOCI, --max_num_loci=MAXIMAL_NUM_LOCI
                    maximal number of gene loci to use as input; if more
                    applied, overflow of informants is discarded
                    [default: 15, range 4..50]
--------------------------------------------------------------------------------------

    Data presented in the manuscript was obtained with default settings
    of in total 15 loci, and at most --maximal_num_loci are kept *.
    The actual number of informants that is taken along after this informant
    filtering can actually be higher. Provided unigene evidence attached to
    informant AbgpGeneLocusDirectories are taken along as extra informants.

    * In the application as published in the manuscript, --filewithloci were used as input,
    which were ordered based on (estimated) similarity to the target locus.
    As such, --maximal_num_loci=15 discarded only the more distantly related
    AbgpGeneLocusDirectory informants.
    
    Lowering this threshold will increase speed at the cost of a putative
    reducting in prediction accuracy, depending on the quality of the informants.
    When ample close informants are available, preferably some of them with unigene
    evidence, lowering this threshold is an interesting speedup. 
     
    Increasing this threshold could yield performance increase.
    In case you offer plain fasta, is it recommended yo use a higher threshold,
    because some of the best from your randomly ordered informants could be removed.
    An alternative, in case you have prior knowledge on the phylogeny of offered informants,
    you could protect some of the most closely related informants from removal by
    listing them in --required_informants
    
    None: running ABFGP on (much) more than 15 informants occasionally yields a drastic
    increase in calculation time. In case a problematic area occurs in the
    GenestructureOfCodingBlockGraphs, ABFGP will try *very* hard to resolve 
    this area (by making several HMMs by applying a leave-one-out strategy in
    search for putative outliers).

--------------------------------------------------------------------------------------
--required_informants
                    force these informants to be selected (from a longer
                    list of of informant genes): --required_informants
                    <infA> ... <infX>
--------------------------------------------------------------------------------------
                        
    Any informant listed in --required_informants won't be discarded, even
    if this might result in a higher number than --maximal_num_loci.
    Useful if you want to enforce that certain informants are maintained.
                        
--------------------------------------------------------------------------------------
--omit_informants   remove these informants prior to informant
                    selection(from a longer list of of informant genes):
                    --omit_informants <infA> ... <infX>
--------------------------------------------------------------------------------------

    Any informant listed in --omit_informants will be discarded, even
    if this might result in a lower number than --minimal_num_loci.
    Useful if you want to study the influence of a certain informant
    on the outcome of ABFGP, or in case you know/discover that
    the informant is actually an outlier that disturbs the outcome.
    
--------------------------------------------------------------------------------------
--disallow_informant_deletion
                    disallow automatic deletion of (very) poor informants
--------------------------------------------------------------------------------------

    Informants can be deleted during ABFGP execution
    a) in an early stage, in the informant selection step (as indicated by --maximal_num_loci)
    b) in a later stage, when the overall contribution of all informants is compared.
    
    In case informant(s) contribute very poorly, they are succesively discarded,
    This selection filters out gene loci that are not homologous enough,
    and highly unlikely represent orthologous loci.
    It is recommended not to disable the deletion of informants.
    Using --disallow_informant_deletion doesn't alter the behaviour of the
    informant selection as imposed by --maximal_num_loci

--------------------------------------------------------------------------------------
--omit_unigene_informants
                    do NOT add unigenes as ABGP informants
--------------------------------------------------------------------------------------

    When used, unigene data provided in the AbgpGeneLocusDirectory folder(s)
    will not be used as informants. Hardly considered useful
    except for testing & development purpose.
 
##
## Options that specify output data
##

--------------------------------------------------------------------------------------
-o OUTDIR, --outdir=OUTDIR
                    specify outdir where to create result file(s); auto-
                    generated when not specified
--------------------------------------------------------------------------------------

   Please read the explanation in teh General options section
   
--------------------------------------------------------------------------------------
--storetoggbdb      store results in GGB (MYSQL) database [default:
                    False]; specify DB_CONNECTION DB_NAME DB_HOST DB_USER
                    DB_PASS in settings.ggb
--------------------------------------------------------------------------------------

    In case you choose the configure a remote MySQL database of the
    Generic Genome Browser format (http://gmod.org/wiki/GBrowse), this
    option will that the GFF output of ABFGP is immediately stored into this database.
    ABFGP uses the perl script load_gff.pl from BioPerl for this.
    It is distributed along in the software/ directory.
    http://www.math.lsa.umich.edu/~dburns/548/bioperl-1.2/scripts/Bio-DB-GFF/load_gff.pl
    
    For this functionality, please install MySQL and GGB according to specifications
    provided by these software. Additional, configure access to your database
    by describing the following variables in settings/ggb.py:
    
    DB_CONNECTION = 'dbi:mysql'
    DB_NAME       = 'XXXXXXXXXX' 
    DB_HOST       = 'XXXXXXXXXX' 
    DB_USER       = 'XXXXXXXXXX' 
    DB_PASS       = 'XXXXXXXXXX' 

    Note: if you want to use a different output database that supports the GFF format,
    you can very easily modify the code. In the final few lines of abfgp.py,
    the functions dbcleanup() and gff2database() are called:

    dbcleanup(fref=used_fref)
    ...
    status = gff2database(gff_fname,fasta_fname)
    
    The function gff2database() are is nothing more than a proxy to the gff2db() function in gff/__init__.py.
    The function ggf2db() just layouts the command line to load_gff.pl.
    The related function dbcleanup() first removes this accession by a MySQL query.

    To make your own desired link, just adjust these 2 lines of code in the main abfgp.py script,
    and make shure your own database link is properly functioning.

    from SOMEWHERE import MY_OW_FUNCTION_THAT_DOES_dbcleanup, MY_OW_FUNCTION_THAT_DOES_gff2database
    MY_OW_FUNCTION_THAT_DOES_dbcleanup(fref=used_fref)
    ...
    status = MY_OW_FUNCTION_THAT_DOES_gff2database(gff_fname,fasta_fname)
    
--------------------------------------------------------------------------------------
--storelocitoggbdb  store raw AbgpGeneLoci data in GGB (MYSQL) database
                    [default: False]; specify DB_CONNECTION DB_NAME
                    DB_HOST DB_USER DB_PASS in settings.ggb
--------------------------------------------------------------------------------------

    In case you are studying a particular target gene, and you provided
    AbgpGeneLocusDirectories as input, all loci are stored in your
    MySQL-coupled Generic Genome Browser database in case you
    did configure it. For notes on the Generic Genome Database, please
    read the explanation in --storetoggbdb above.


#####################################################################################
#### 9. How to use ABFGP in your research                                        ####
#####################################################################################
####                                                                             ####
#### Here, we ventilate our thought and advise on how to use ABFGP.              ####
#### It mainly depends on the SCALE of which you want to employ it               ####
####                                                                             ####
#### Please make shure that have read section 8. (command line arguments) first. ####
####                                                                             ####
#####################################################################################

-------------------------------------------------------------------------------------             
1)  Manual curation assistance for a single / few genes
    use --multifasta, --loci, --dirwithloci
-------------------------------------------------------------------------------------             

    In this case, creating AbgpGeneLocusDirectories might be to much work for job to be done.
    
    - Select an extensive bunch of informant species
    - Download their genome, proteome, optionally GFF documentations
    - Download unigene data of the informant species, is available.
      Unigene data, even of distant informants, can be of great contribution to ABFGP performance!
    - make BLAST databases per genome / proteome / transcriptome
    - link target proteins from your gene catalogue to predicted proteins in the informant proteomes (BLASTP)
      In the next section ( 2) Genome-wide re-annotation ) some extra notes are provided
    - obtain the genomic locus, including ~750-1500nt of flanking sequence, at which the protein is
      encoded in the informant genome. You can use the GFF gene annotation as a proxy. If you are lazy,
      you could use TBLASTN to find the locus in its genome.
      There are a dozen way to obtain the flanking sequence.
      You could use slopBed: we are a great fan of the BEDTools suite!
      http://bedtools.readthedocs.org/en/latest/
      http://bedtools.readthedocs.org/en/latest/content/tools/slop.html      
      Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.
    - blast your target protein of interest to downloaded databases with unigenes (TBLASTN)
    - the NCBI's UniGene database could be an interesting and rich resource for plenty of unigenes:
      http://www.ncbi.nlm.nih.gov/unigene
      However, the amount of fungal data in this resource is not staggering.
    - concatenate all sequences to a multi fasta DNA file
    
    The ABFGP command line will eventually look something like this:
    
    abfgp.py --max_num_loci 30 --multifasta INFORMANTS.mfa --dna TARGET.fa --target MY_TARGET_GENE
    abfgp.py --max_num_loci 30 --multifasta INFORMANTS_AND_TARGET.mfa --dna TARGET.fa --target MY_TARGET_GENE
    abfgp.py --max_num_loci 30 --multifasta GENOMIC_INFORMANTS.mfa --dna TARGET.fa UNIGENE_1.fa UNIGENE_2.fa ... UNIGENE_X.fa --target MY_TARGET_GENE

    A note on --max_num_loci: in case you have a single gene to inspect, you just want it to
    be as properly predicted as possible. By providing genomic informants and (partial) unigene
    informants as fasta, all above --max_num_loci will be discarded. Take the risk
    of a bit longer calculation time, and run ABFGP with (much) more informants.
    
    There is not much difference in effort between obtaining a genomic
    informant locus in fasta format and making a AbgpGeneLocusDirectory of it;
    please read the next subsection 2) Genome-wide re-annotation or read
    section 6. (Explaining the AbgpGeneLocusDirectory) of this README file.
    An option is to store all genomic informant loci aiming at your target gene
    as AbgpGeneLocusDirectories in a single folder. As such, you could
    place the unigene informants (form the downloaded unigene datasets and/or
    obtained by querying the NCBI UniGene database) in a fasta file.
    In this case, your command line will look like this:
    
    abfgp.py --max_num_loci 30 --multifasta path/to/LOCI_FOR_TARGET_GENE/UNIGENE_INFORMANTS.mfa --dirwithloci path/to/LOCI_FOR_TARGET_GENE --target MY_TARGET_GENE
      
-------------------------------------------------------------------------------------             
2)  Genome-wide re-annotation of a (single) target gene catalogue by the ABFGP method
    use --dirwithloci, or better --filewithloci, but we discourage --multifasta or --dna
-------------------------------------------------------------------------------------

    This is approaching the scale as presented in the manuscript.
    Be prepared for a considerable amount of preparatory work.
    The most time consuming step will be the various BLAST steps to infer informants.
    You would have to acquire informant gene loci on a gene-by-gene basis
    for each gene in your gene catalogue.
    
    - Select an extensive bunch of informant species
    - Download their genome, proteome, optionally GFF documentations
    - Download unigene data of the informant species, is available.
      Unigene data, even of distant informants, can be of great contribution to ABFGP performance!
    - A note on RNA-Seq data: ABFGP can currently only deal with (assembled) transcript data.
      If you have RNA-Seq data for target and/or informant species, convert it
      to inferred transcripts first (e.g. by using Cufflinks; http://cufflinks.cbcb.umd.edu/)
    - make BLAST databases per genome / proteome / transcriptome
    - obtain informant loci:
    
        - link target proteins from your gene catalogue to predicted proteins in the informant proteomes (BLASTP)
        
          In the method section of the manuscript, we described our approach to enrich for likely
          orthologous genes. Lowering these thresholds a bit is fine, but make shure you do enrich 
          for homologous proteins, not just similarity based on a single domain hit.
          In fact, there is no strict reason NOT to include paralogs, as soon as the overall
          gene sequences are not to divergent (e.g. the same series of domains, approximately the same length, etc).
          N.b. there are some extra remarks on retrieving ABFGP informants in the next
          subsection ( 3) Genome-wide re-annotation of several gene catalogues ).
          
        - obtain the genomic locus, including ~750-1500nt of flanking sequence, at which the protein is
          encoded in the informant genome. You can use the GFF gene annotation as a proxy. If you are lazy,
          you could use TBLASTN to find the locus in its genome.
          There are a dozen way to obtain the flanking sequence.
          As stated above, you could use slopBed: we are a great fan of the BEDTools suite!
          
        * alternatively, you could obtain genomic loci directly by using your target protein
          as a query on the informant genome as a database (with TBLASTN). This will work fine for
          (quite) closely related informant species / quite conserved genes or gene families, 
          However, we do not strongly recommend to use this option. Gene prediction software
          might yield inaccurate gene models, they are perfectly capable of assigning gene loci.
          A whole series of issues could arise when using the direct TBLASTN approach:
          longer introns (which break up hits), genes with many introns, etc.

          * in addition, search with TBLASTN in the genome for potentially missed genes;
          do not forget to add some flanking sequence too! Be cautious not to include
          loci that comprise already annotated genes! Doublet informant data is
          in the best case not useful at all, and in the worse case might negatively
          influence the outcome because its sequence context will have more weight in
          the all-vs-all comparison steps undertaken in the ABFGP method.
          
    - now it is a matter of decision: how do you organise your data?
        Option a) and b) are not recommended. See Concluding remarks ~50 lines below why.
        a) make a multi fasta file for each target gene locus, and concatenate your
           informant gene loci to this file.
           Eventually use abfgp.py --multifasta MY_FASTA_FILE
        b) make a directory per gene and fill it with single fasta files
           of target and informant gene loci.
           Eventually use abfgp.py --dna FILE1 FILE2 ... FILEX
        c) make a directory per gene and, and place each single fasta file
           of target and informant gene loci in a separate folder in this directory.
           In other words, make a directory of AbgpGeneLocusDirectories.
           For full explanation and requirements, see section 6. of this README.
           `*.dna*.fa`      for the DNA of genomic loci in fasta format
           `*.locus.gff`    describing the locus on the genome

           MY_TARGET_GENE/
                MY_TARGET_GENE/
                    genomiclocus.dna.fa
                    genomiclocus.locus.gff
                INFORMANT1/
                    genomiclocus.dna.plus100nt.fa
                    genomiclocus.locus.gff
                INFORMANT2/
                    genomiclocus.dna.including1500nt_flank.fa
                    genomiclocus.locus.gff
                INFORMANT3/
                    locus_from_direct_tblastn.dna.fa
                    genomiclocus.locus.gff
           
           Please look at testdata/exampleAbgpLocusDirs/CFU_827306 for an example.
           Just use abfgp.py --dirwithloci path/to/MY_TARGET_GENE 

    In case unigene data for target or informant species are available, handle accordingly:
    
    - map the unigenes on the genome with your preferred software (e.g. GeneSeker, PASA)
    - link the loci at which unigenes are aligned to loci at which genes are annotated.
      There are a dozen way to do this. We are a great fan of the BEDTools suite:
      http://bedtools.readthedocs.org/en/latest/
      Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics. 26, 6, pp. 841–842.
      Use the executable clusterBed; cat <FILE_WITH_GENE_ANNOTATION.GFF> <FILE_WITH_ALIGNED_UNIGENES.GFF | clusterBed -i stdin
    - depending on the choice you organized your data, handle the unigene that is linked
      that linked to your target or informant gene locus:
      a) retrieve its fasta accession (e.g. EMBOSS.seqret), and concatenate to your multi fasta.
         This unigene will be taken along as any other informant.
         Make shure you use a (much) higher --max_num_loci **
      b) retrieve its fasta accession (e.g. EMBOSS.seqret), and place as a single fasta in your
         input directory (1 unigene accession per file!). This unigene will be taken along as any other informant.
         Make shure you use a (much) higher --max_num_loci **
      c) place the GFF describing the alignment to the GENOMIC coordinates (not relative to the gene locus!)
         in the proper subdirectory. Please check section 6. of this README for a full example.

    Concluding remarks:
    
    The advantage of using --dirwithloci is that your input data is the extra functionality
    of the AbgpGeneLocusDirectories (e.g. exploiting prior knowledge on the gene models)
    and that it is very easy to manage in case of updates. An extended or new unigene dataset
    can just be added to the existing AbgpGeneLocusDirectories, and you could re-run ABFGP.
    If you weren't satisfied by the overall predictions of the ABFGP method on current
    sampled set of informant genomes, add a few more informants species and re-run ABFGP.
    It is only a little bit more effort to switch from (multi) fasta files
    to AbgpGeneLocusDirectories, and it will yield a better performance.
    In the 3th described implementation of preparing input data, it is explained
    why you could as easily have used --filewithloci, and organize al your
    informant genes per species, and use the --filewithloci as a proxy to make the
    homologous/orthogous relationship. The --filewithloci option is explained in detail
    in the following, and its option 

    ** A big disadvantage of using unigenes in the --dna / --multifasta mode is that 
       ABFGP isn't aware of the biological link between the unigene an its genomic locus.
       By (default) filtering on --max_num_loci=15, you might easily discard the precious
       unigenes that you added with al this extra effort!
    
-------------------------------------------------------------------------------------             
3)  Genome-wide re-annotation of several gene catalogues
    use --filewithloci, definately not --multifasta
-------------------------------------------------------------------------------------
    This is the scale and mode as presented in the manuscript.
    All points and preparation flow as mentioned in case 3) are valid.
    We have a few extra remarks:

    - we strongly advise to start organizing your project with a per-species
      organized collection of AbgpGeneLocusDirectories of complete gene catalogues.
      Ample reasons are given in the `Concluding remarks` ~30 lines above.
      See section 6. of this README ( AbgpGeneLocusDirectory ) for how to easily achieve this.
    - ABFGP was not tested with *several* informants that are closely related
      to the target genome. The closest related informants used were
      Dothistroma septosporum for Cladosporium fulvum and vise versa,
      and Sclerotinia sclerotiorum for Botrytis cinerea. Please read in their
      genome publications on estimated identity on the DNA level. Including
      a whole bunch of informants that are >~95% on the whole-genome scale
      (at the DNA level) will most likely not yield better / good results.
    - make an all-vs-all similarity matrix of target and informant proteins first.
      An example is worked out in the description of the command line parameter
      --filewithloci in section 8. of this README. You do not have a single,
      but multiple query proteomes, so decoupling this step will make it better
      scalable and incrementable (with extra informants).
    - make shure you have mappable / traceable protein identifiers to gene loci.
      This enables for a direct translation from a query in your similarity matrix
      to paths on your filesystem; voila, you have input that --filewithloci understands.
    - When having a proper similarity matrix, it is actually interesting to play around
      with more sophisticated methods that infer orthology or classify orthologous groups.
      Plenty examples here, but one of the most commonly used is OrthoMCL:
      http://orthomcl.org/orthomcl/
      Li Li, Christian J. Stoeckert, Jr., and David S. Roos
      OrthoMCL: Identification of Ortholog Groups for Eukaryotic Genomes
      Genome Res. 2003 13: 2178-2189.

-------------------------------------------------------------------------------------             
4)  Integrated genome-wide re-annotation by ABFGP in existing gene annotation pipelines
-------------------------------------------------------------------------------------
    We advise to work according to our recommendation in the above sections.
    The by the ABFGP method re-annotated gene model could be added as an additional track
    of evidence to offer to your favourite gene annotation combiner tool (e.g. EVM, Evigen,
    MAKER, Jigsaw). In case your favourite combiner tools allows a manually adjustable
    evidence level per track, ABFGP acquired models should be rewarded higher as
    simple protein-alignments. In case ab-initio gene predictors are used, that did not
    got any similarity- of expression-hints, ABFGP is likely to have an overall higher
    reliability as these tools.
    An interesting stategy could be an interative approach. First make `perfect` gene
    models with your combiner tool (possibly usinf ABFGP as an evidence track), and feed
    these models including their informant gene loci (again) to ABFGP. Do not
    use the --abinitio mode, but use prior knowlegde of the gene models.
    All models that are exactly re-predicted could now we with great confident be
    places aside not requiring manual annotation. Those gene models for which ABFGP
    did propose revisions, especially in case no or little transcript evidence is available,
    are a category to inspect manually.

    
#####################################################################################
#### 10. The other ABFGP executables (examples)                                  ####
#####################################################################################
####                                                                             ####
#### During the development of the ABFGP method, various helper executables      ####
#### were made. They were used to test (parts of) the code and generate          ####
#### information on (groups of) AbgpGeneLocusDirectories. This is why they could ####
#### be of interest to you.                                                      ####
####                                                                             ####
#### abgp_geneconfirmation.py                                                    ####
#### abgp-genestructure.py                                                       ####
#### abgp-intronanalyses.py                                                      ####
#### abgp-orfonanalyses.py                                                       ####
####                                                                             ####
#### Naming is either abgp_xxxx.py or abgp-xxxx.py. Underscores mean that        ####
#### these executables are as well libraries from which abfgp.py imports         ####
#### required functions and classes. Dashes mean that these executables are      ####
#### small stand-alone units that perform some tasks (that are often employed    ####
#### within the ABFGP method too).                                               ####
####                                                                             ####
#### Below follows a brief explanation-by-example.                               ####
####                                                                             ####
#####################################################################################

-------------------------------------------------------------------------------------             
abgp_geneconfirmation.py    find inconsistencies & abnormalities in (annotated)
                            gene models provided in a AbgpGeneLocusDirectory
-------------------------------------------------------------------------------------

python abgp_geneconfirmation.py -d testdata/exampleAbgpLocusDirs/CFU_827701/DOTSE_070160/ -v
# EOF geneconfirmation; if nothing printed, no abnormalities observed

python abgp_geneconfirmation.py -d testdata/exampleAbgpLocusDirs/CFU_827701/CRYPA_341138
# WARNING: Non-canonical splice site(s) for organism/fref: CRYPA_341138 CRYPA_341138
WARNING::NonCanonicalSpliceSiteWarning(Donor:GC)

-------------------------------------------------------------------------------------             
abgp-genestructure.py       gives per-exon and per-intron statistics of the gene model;
                            Warnings are given in case of `features` in gene models
-------------------------------------------------------------------------------------             

python abgp-genestructure.py --filewithloci testdata/exampleFilewithloci/CFU_829726.relativepaths.bdbh.csv 

# this prints detailed description for each of the AbgpGeneLocusDirectories;
# - length and AT% of each exon and intron
# - PSSM quality of TSS, donor and acceptor sites; donor and acceptor sequence context
# - ORFs on which the exon(s) are located and their total length
# Warnings are printed in case gene structure(s) have features that could impact ABFGP performance on these loci;
# small and tiny exons are an example for this, but as well annotated gene models that are incomplete
# - WARNING::IncompleteGeneStructureWarning 
# - WARNING::SmallAnnotatedExonWarning
# - WARNING::SmallAnnotatedFirstExonWarning
# - WARNING::SmallAnnotatedFinalExonWarning
# - WARNING::NonCanonicalSpliceSiteWarning
# - etc... see abgp_warnings.py for a full list of warnigns that could be thrown here to stdout.

-------------------------------------------------------------------------------------
abgp-intronanalyses.py      gives detailed statistics per intron in the gene model
                            provided in the AbgpGeneLocusDirectory
                            -e prints explanation on the 19 column output
                            -v prints a 'EOF intronstatistics' message; this can help
                            you in case no introns / no gene model is given in the 
                            AbgpGeneLocusDirectory to see that there is (no) output.
-------------------------------------------------------------------------------------

python abgp-intronanalyses.py -d testdata/exampleAbgpLocusDirs/CFU_827701/AO090012000150/ -v
# No unigene GFF data present in GeneLocusDirectory
AO090012000150  AO090012000150  3       1       None    1       True    2143    2167    24      0.50    -0.09   4.34    0.82    ?       2       2       1       14
AO090012000150  AO090012000150  3       2       None    2       False   2393    2442    49      0.51    5.68    4.64    0.89    ?       2       2       1       15
AO090012000150  AO090012000150  3       3       None    2       False   2811    2865    54      0.35    3.73    6.46    0.95    C       5       3       1       14
# EOF orfstatistics

python abgp-intronanalyses.py -d testdata/exampleAbgpLocusDirs/CFU_827701/CRYPA_341138/
Warning::NonCanonicalDonor GCAAGATTTACCAGACCGCACCCACCACCCGTCCCATCCTAGACATGAGAGGGCATCTGGCTCTGCCGTGATCTGTCTCTAACCAACATCTTTTCCGGGAGCAG
CRYPA_341138    CRYPA_341138    3       1       None    1       False   1615    1719    104     0.56    8.45    3.22    0.80    ?       4       3       1       26
CRYPA_341138    CRYPA_341138    3       2       None    0       False   2006    2085    79      0.44    5.99    5.95    1.08    C       3       2       1       27
CRYPA_341138    CRYPA_341138    3       3       None    1       False   2827    2892    65      0.46    3.69    5.13    0.67    N       3       1       1       17


-------------------------------------------------------------------------------------
abgp-orfonanalyses.py       gives detailed statistics per - foward strand - ORF
                            in the provided locus. ORFs that are part of the
                            annotated CDS of the gene in case provided in the
                            AbgpGeneLocusDirectory, are flagged as True in the
                            11th column of the output.
                            -e prints explanation on the 11 column output
-------------------------------------------------------------------------------------

python abgp-orfanalyses.py -d testdata/exampleAbgpLocusDirs/CFU_827701/CRYPA_341138/ | grep True
CRYPA_341138    CRYPA_341138    1       2085    2898    868.4   1       813     ?       1.068   True    0.00
CRYPA_341138    CRYPA_341138    3       2879    3386    575.0   1       507     C       1.134   True    0.01
CRYPA_341138    CRYPA_341138    4       1697    2213    560.8   2       516     C       1.087   True    0.00
CRYPA_341138    CRYPA_341138    18      1479    1683    193.1   2       204     ?       0.947   True    0.00


#####################################################################################
#### 11. Abbreviations used in ABFGP code and documentation                      ####
#####################################################################################
####                                                                             ####
#### The ABFGP source code decorated with a variable degree of documentation.    ####
#### Some sections are excellently documented, including docstrings.             ####
#### Other sections are only briefly documented.                                 ####
####                                                                             ####
#### In the code, variables are often named by an abbreviation of                ####
#### the Python class the object is an instance off. Below is a (non exlcusive)  ####
#### listing of most often used abbreviations and some minor explanation         ####
####                                                                             ####
#####################################################################################


ABFGP           Alignment Based Fungal Gene Prediction

ABGP            Alignment Based Gene Prediction;
                Subtle modified versions of AB(F)GP were successfully tested by the
                developer on non-fungal systems too. Modules in the code that generic
                for alignment-based gene prediction, and are not tailored to fungi,
                are named accordingly.
            
AbgpGeneLocus   AbgpGeneLocusDirectory (just a shortening)

# Various types of graph objects that are used throughout ABFGP
# that describe multiple (sequence) aligned entities. 
            
CB	            Coding Block
CBG	            CodingBlockGraph
CBE	            CodingBlockEnd (gff object)
CBS	            CodingBlockStart (gff object)
PCG             PacbpCollectionGraph
PACBP	        Pairwise-Aligned Coding Block PROTEIN
PACBPORF        Pairwise-Aligned Coding Block PROTEIN with Orf objects
lsrCBG	        Low Similarity Region Coding Block Graph
GSG	            GenestructureOfCodingBlockGraphs
GTG             GeneTreeGraph
partGSG	        partial GenestructureOfCodingBlockGraphs; part of the GSG that is needed for current operation(s)
OMSR	        Overall Minimal Spanning Range
MSR	            Minimal Spanning Range
SPRDIF          (overall Minimal) Spanning Range Difference
MAXSPR	        Maximum Spanning Range
MINSPR	        Minimum Spanning Range (== OMSR??)
cbgIF	        CodingBlockGraph Interface
ECG	            ExonCollectionGraph
FEG	            FirstExonGraph; an ExonCollectionGraph with TSS exons (TSS--acceptor i.s.o donor--acceptor)
TSS	            TranslationalStartSite (gff object)

firstCBG	    the CBG with TSS for each of the organism identifiers; start of the GSG
lastCBG		    the CBG which contains the AlignedStopCodonGraph; end of the GSG
finalCBG	    the CBG which contains the AlignedStopCodonGraph; end of the GSG
donorCBG	    CBG that delivers donor sites in a cbgIF (5' side / left of a acceptorCBG)
cbgD		    CBG that delivers donor sites in a cbgIF (5' side / left of a acceptorCBG)
acceptorCBG	    CBG that delivers acceptor sites in a cbgIF (3' side / right of a donorCBG)
accepCBG	    CBG that delivers acceptor sites in a cbgIF (3' side / right of a donorCBG)
cbgA		    CBG that delivers acceptor sites in a cbgIF (3' side / right of a donorCBG)
orfD	        Open Reading Frame object with / must deliver a donor site
orfA	        Open Reading Frame object with / must deliver an acceptor site

ADG	            AlignedDonorSiteGraph
ADSG	        AlignedDonorSiteGraph
AAG	            AlignedAcceptorSiteGraph
AASG	        AlignedAcceptorSiteGraph
DSCG	        DonorSiteCollectionGraph
DCG	            DonorSiteCollectionGraph
ASCG	        AcceptorSiteCollectionGraph
ACG	            AcceptorSiteCollectionGraph

this		    current CBG instance in a adjacency comparison on both sides (prev,this,next)
thisCBG		    current CBG instance in a adjacency comparison on both sides (prev,this,next)
next		    next CBG instance in a adjacency comparison (prev,(this),next)
nextCBG		    next CBG instance in a adjacency comparison (prev,(this),next)
prev		    previous CBG instance in a adjacency comparison (prev,(this),next)
prevCBG		    previous CBG instance in a adjacency comparison (prev,(this),next)

# typically, often used variables.

orf	            Open Reading Frame (gff object), defined from stop to stop
org	            Organism/gene identifier (string)
orgQ	        Organism/gene identifier (string) that serves as Query in a comparison
orgS	        Organism/gene identifier (string) that serves as Subject in a comparison
stw	            StopWatch timer object for performance logging

