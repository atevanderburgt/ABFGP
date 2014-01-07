# generic functions for python Optparse module
from optparse import OptionValueError
import os.path

def isdir_callback(option, opt_str, value, parser):
    """ is the specified path a directory? """
    if not value:
        # no path specified -> set default!
        setattr( parser.values, option.dest, None )
    elif os.path.isdir(value):
        # set os.path.abspath of directory name
        setattr( parser.values, option.dest, os.path.abspath(value) )
    else:
        if os.path.exists(value):
            raise OptionValueError("'%s' is not a directory (%s)" % (value,opt_str))
        else:
            raise OptionValueError("'%s' does not exist (%s)" % (value,opt_str))

# end of function isdir_callback


def isfile_callback(option, opt_str, value, parser):
    """ is the specified path a file? """
    if not value:
        # no file specified -> set default!
        setattr( parser.values, option.dest, option.default )
    elif os.path.isfile(value):
        # set os.path.abspath of file name
        setattr( parser.values, option.dest, os.path.abspath(value) )
    else:
        raise OptionValueError("'%s' in not a (regular) file (%s)" % (value,opt_str))

# end of function isfile_callback


def vararg_callback(option, opt_str, value, parser):
    """ capture variable number of non-digit arguments """
    assert value is None
    value = []
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a like options
        if arg[:1] == "-" and len(arg) > 1:
            break
        value.append(arg)
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

# end of function vararg_callback


def vararg_callback_choices(option, opt_str, value, parser):
    """ capture variable number of non-digit arguments from list of choices """
    #assert value is None
    value = [ value ] 
    for arg in parser.rargs:
        # stop on --foo like options
        if arg[:2] == "--" and len(arg) > 2:
            break
        # stop on -a like options
        if arg[:1] == "-" and len(arg) > 1:
            break
        if arg in value:
            choices = ", ".join([ "'%s'" % choice for choice in option.choices ])
            raise OptionValueError("option %s: duplicate choice '%s' in arguments" % (opt_str,arg))
        elif arg in option.choices:
            value.append(arg)
        else:
            choices = ", ".join([ "'%s'" % choice for choice in option.choices ])
            raise OptionValueError("option %s: invalid choice '%s' (choose from %s)" % (opt_str,arg,choices))
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

# end of function vararg_callback


def valueinrange_callback(option, opt_str, value, parser, minval=None, maxval=None):
    """ check if an int or float value lies in a specifix range [minval,maxval]"""
    if float(value) >= minval and float(value) <= maxval:
        setattr( parser.values, option.dest, value )
    else:
        raise OptionValueError(
            "keyword argument '%s' not in range %s..%s" % (
            opt_str,minval,maxval ))

# end of function valueinrange_callback

