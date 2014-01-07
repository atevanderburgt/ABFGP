#!/usr/bin/python
""" Quick test if all the third party software required for AB(F)GP are properly present,executable and (hopefully) okay """
__authors__ = "Ate van der Burgt"
__license__ = "MIT"

import os
from settings.executables import *

# All required executable names that should be declared in settings.executables
required_executable_names = [
    'EXECUTABLE_BLASTALL',
    'EXECUTABLE_CEXPANDER_ALLVSALL',
    'EXECUTABLE_CEXPANDER_CBALIGNP',
    'EXECUTABLE_CEXPANDER_CEXPANDER',
    'EXECUTABLE_CLUSTALW',
    'EXECUTABLE_FASTALENGTH',
    'EXECUTABLE_FORMATDB',
    'EXECUTABLE_GETORF',
    'EXECUTABLE_GFF2FASTA',
    'EXECUTABLE_GFFLENGTH',
    'EXECUTABLE_HMMBUILD',
    'EXECUTABLE_HMMSEARCH',
    'EXECUTABLE_LOAD_GFF',
    'EXECUTABLE_SFM',
    'EXECUTABLE_SIGNALP',
    'EXECUTABLE_TCODE',
    'EXECUTABLE_TMHMM',
    'EXECUTABLE_TRANSEQ',
    'EXECUTABLE_UNIGENEANNOTATION',
    ]

# list all the required executables with extentions
extention_2_basename = {
    '.py':  ('get_gff_sequence_from_fasta.py','prep_launch.py','unigeneannotation.py'),
    '.sh':  ('fastalength.sh','gfflength.sh'),
    '.pl':  ('load_gff.pl',),
    }
    
# After calling a command line with --help/-h, nearly all executable
# list their name somewhere in stdout or stderr. For the executables
# that do not, a uniquely recognizable string from stdout/stderr is
# provided here as a check for the executable's identity
executable_lookup_strings = {
    'EXECUTABLE_TRANSEQ':   'Translate nucleic acid sequences',
    'EXECUTABLE_GETORF':    'Finds and extracts open reading frames (ORFs)',
    }

for varname in dir():
    if not varname.startswith("EXECUTABLE"): continue
    if varname.endswith("_VERSION"): continue
    if varname.endswith("PATH"): continue
    # emboss tcode & getorf settings are listed in executables too;
    # ignore these here
    if varname.find("_GETORF_") > 0: continue
    if varname.find("_TCODE_") > 0: continue

    # remove varname from required_executable_names;
    # check at the end if all required_executable_names are dealth with
    if varname in required_executable_names:
        required_executable_names.remove(varname)
    
    # this is an absolute path to an required executable in ABFGP
    # test if this path exists and, check if it is executable and,
    # if possible, check if it indeed represents the corresponding executable
    variable = eval("%s" % varname)
    
    if not os.path.exists(variable):
        print "False ... %30s -> %s [PATH DOES NOT EXIST]" % (varname,variable)
        continue
    elif not os.access(variable, os.X_OK):
        print "False ... %30s -> %s [PATH IS NOT EXECUTABLE]" % (varname,variable)
        continue


    # get extention & exptected vs observed exe name for consistency check
    extention = os.path.splitext(variable)[-1]
    expected_exe_name = varname.split("_")[-1].lower()
    exe_name = os.path.basename(variable).lower()
    
    if varname in ['EXECUTABLE_CEXPANDER_CEXPANDER']:
        # no easy test (-h,--help) available.
        # But distributed as part of the ABFGP code tree so it assumed to be okay
        stdout,stderr = '',''
        is_program_name_observed = True
        pass
    elif extention in ['.py','.pl','.sh']:
        # no easy test (-h,--help) available.
        # But distributed as part of the ABFGP code tree so it assumed to be okay
        if os.path.basename(variable) in extention_2_basename[extention]:
            stdout,stderr = '',''
            is_program_name_observed = True
            pass
        else:
            print "False ... %30s -> %s [UNEXPECTED %s EXECUTABLE]" % (varname,variable,extention)
            continue
    elif expected_exe_name != exe_name and (expected_exe_name,exe_name) != ('sfm','scan_for_matches'):
        print "False ... %30s -> %s [UNEXPECTED EXECUTABLE NAME: %s / %s ]" % (varname,variable,expected_exe_name,exe_name)
    else:
        is_program_name_observed = False
        for extra_command_line in ('--help','-h',''):
            command = "%s %s" % (variable,extra_command_line)
            ci,co,ce = os.popen3(command)
            ci.close()
            stdout = co.read()
            co.close()
            stderr = ce.read()    
            ce.close()
            is_program_name_observed = (stdout.lower()+stderr.lower()).replace(" ","").find(exe_name) >= 0
            if is_program_name_observed:
                break
            elif executable_lookup_strings.has_key(varname) and (stdout+stderr).find(executable_lookup_strings[varname]) >= 0:
                is_program_name_observed = True
                break
    if is_program_name_observed:
        print "True  ... %30s -> %s" % (varname,variable)
    else:
        print "False ... %30s -> %s [CAN'T VERIFY CORRECT IDENTITY OF THIS EXECUTABLE]" % (varname,variable)
else:
    # check at the end if all required_executable_names
    # are declared in settings.executables
    if required_executable_names:
        for varname in required_executable_names:
            print "False ... %30s -> [REQUIRED EXECUTABLE is not declared in settings.executables]" % (varname)
    
    
    