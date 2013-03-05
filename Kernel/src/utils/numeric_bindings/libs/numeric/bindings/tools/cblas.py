#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#  Copyright (c) 2009 Rutger ter Borg
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
#

import netlib

import re, os.path, copy
from types import StringType

# for debugging purposes
import pprint

def parse_file( filename, info_map, template_map ):

    parser_mode = template_map[ 'PARSERMODE' ]
    prefix_map = { 'BLAS': 'cblas_',
                   'LAPACK': 'clapack_' }
    prefix = prefix_map[ parser_mode ]

    print prefix

    pp = pprint.PrettyPrinter( indent = 2 )
    source = open( filename ).read() 

    for match in re.compile( '(void|float|int|double|CBLAS_INDEX) +' + prefix + '([^\(]+)\(([^\)]+)\)', re.M | re.S ).findall( source ):
        print "----"
        result_type  = match[0]
        fortran_routine = match[1].split("_sub")[0].upper().strip()
        c_routine = prefix + match[1]

        print "C" + parser_mode + " routine:", c_routine , "   " + parser_mode + " equivalent:", fortran_routine
        arguments = {}
        for arg in match[2].replace('\n','').split( ',' ):
            arg = arg.strip()
            arg_name = arg.split( " " )[-1].replace( "*", "" ).strip().upper()
            arguments[ arg_name ] = {}
            arguments[ arg_name ][ "original" ] = arg
            arguments[ arg_name ][ "pointer" ] = "*" in arg

        pp.pprint( arguments )

        if fortran_routine in info_map:
            print "Found ", fortran_routine, " in Fortran info_map."
            info_map[ fortran_routine ][ prefix + "routine" ] = c_routine

            # read aliases, if they are there
            my_key = info_map[ fortran_routine ][ 'group_name' ].lower() + '.all.cblas_alias'
            alias_map = {}
            #print my_key
            if netlib.my_has_key( my_key, template_map ) != None:
                #print "Has key.."
                for line in template_map[ netlib.my_has_key( my_key, template_map ) ].splitlines():
                    #print "Line:", line
                    alias_map[ line.split( "," )[0] ] = line.split(",")[1]

            #print alias_map

            # Try to match and insert arguments
            # argument_map is the data gathered through the Fortran interface
            for arg in info_map[ fortran_routine ][ 'argument_map' ]:
                cblas_arg = ''
                if arg in arguments:
                    cblas_arg = arg
                elif arg in alias_map:
                    if alias_map[ arg ] in arguments:
                        cblas_arg = alias_map[ arg ]

                print "Looking for " + parser_mode + " argument ", arg, " CBLAS equivalent: ", cblas_arg
                if cblas_arg in arguments:
                    print "Found matching argument, inserting call_cblas_header stuff"
                    call_cblas_header = info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_blas_header" ]
                    print "Original: ", call_cblas_header
                    if not arguments[ cblas_arg ][ "pointer" ]:
                        call_cblas_header = call_cblas_header.replace( "&", "" )

                    call_cblas_header = call_cblas_header.replace( "complex_ptr", "void_ptr" );

                    print "Result:   ", call_cblas_header
                    if arg == 'UPLO':
                        info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_" + prefix + "header" ] = \
                            prefix + "option< UpLo >::value"
                    elif arg == 'DIAG':
                        info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_" + prefix + "header" ] = \
                            prefix + "option< Diag >::value"
                    elif arg == 'SIDE':
                        info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_" + prefix + "header" ] = \
                            prefix + "option< " + netlib.template_parameter[ arg ] + " >::value"
                    elif  arg == 'TRANS' or arg == 'TRANSA' or arg == 'TRANSB':
                        info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_" + prefix + "header" ] = \
                          prefix + "option< " + netlib.template_parameter[ arg ] + " >::value"
                    else:
                        info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_" + prefix + "header" ] = call_cblas_header
                else:
                    if arg == 'INFO' and result_type == 'int':
                        info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_" + prefix + "header" ] = None
                        print "INFO is the return type, adding it with code None"
                    elif arg == 'WORK' or arg == 'LWORK':
                        info_map[ fortran_routine ][ "argument_map" ][ arg ][ "code" ][ "call_" + prefix + "header" ] = None
                        info_map[ fortran_routine ][ "clapack_disable_workspace" ] = True
                    else:
                        exit(0)

            if "ORDER" in arguments:
                print "Adding order argument."
                info_map[ fortran_routine ][ "has_" + prefix + "order_arg" ] = True
            else:
                print "Not adding order argument."
                info_map[ fortran_routine ][ "has_" + prefix + "order_arg" ] = False


