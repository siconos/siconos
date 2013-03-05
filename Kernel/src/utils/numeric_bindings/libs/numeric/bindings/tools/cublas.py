#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#  Copyright (c) 2009 Rutger ter Borg
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
#

import re, os.path, copy, netlib
from types import StringType

# for debugging purposes
import pprint

def parse_file( filename, info_map, template_map ):
    pp = pprint.PrettyPrinter( indent = 2 )
    source = open( filename ).read() 

    for match in re.compile( '^(cuComplex|cuDoubleComplex|float|double|void|int) ?CUBLASAPI ?cublas([SDCZI][a-z0-9]+) ?\(([^\)]+)\)', re.M | re.S ).findall( source ):
        print "----"

        result_type  = match[0]
        blas_routine = match[1].upper().strip()
        print "CUBLAS routine:", match[1], "   BLAS equivalent:", blas_routine

        arguments = {}
        for arg in match[2].replace('\n','').split( ',' ):
            arg = arg.strip()
            arg_name = arg.split( " " )[-1].replace( "*", "" ).strip().upper()
            arguments[ arg_name ] = {}
            arguments[ arg_name ][ "original" ] = arg
            arguments[ arg_name ][ "pointer" ] = "*" in arg

        pp.pprint( arguments )

        if blas_routine in info_map:
            print "Found ", blas_routine, " in Fortran info_map."
            info_map[ blas_routine ][ "cublas_routine" ] = 'cublas' + match[1]
            #pp.pprint( info_map[ blas_routine ] )

            # read aliases, if they are there
            my_key = info_map[ blas_routine ][ 'group_name' ].lower() + '.all.cblas_alias'
            alias_map = {}
            print my_key
            if netlib.my_has_key( my_key, template_map ) != None:
                #print "Has key.."
                for line in template_map[ netlib.my_has_key( my_key, template_map ) ].splitlines():
                    #print "Line:", line
                    alias_map[ line.split( "," )[0] ] = line.split(",")[1]

            for arg in info_map[ blas_routine ][ 'argument_map' ]:
                cublas_arg = ''
                if arg in arguments:
                    cublas_arg = arg
                elif 'S' + arg in arguments:
                        cublas_arg = 'S' + arg
                # E.g., BLAS DPARAM equals CUBLAS SPARAM
                elif 'S' + arg[1:] in arguments and arg == 'DPARAM':
                        cublas_arg = 'S' + arg[1:]
                elif 'C' + arg in arguments:
                        cublas_arg = 'C' + arg
                elif arg in alias_map:
                    if alias_map[ arg ] in arguments:
                        cublas_arg = alias_map[ arg ]

                print "Looking for BLAS argument ", arg, " CUBLAS equivalent: ", cublas_arg

                if cublas_arg in arguments:
                    print "Found matching argument, inserting call_cublas_header stuff"
                    call_cublas_header = info_map[ blas_routine ][ "argument_map" ][ arg ][ "code" ][ "call_blas_header" ]

                    print "Original: ", call_cublas_header
                    if not arguments[ cublas_arg ][ "pointer" ]:
                        call_cublas_header = call_cublas_header.replace( "&", "" )

                    call_cublas_header = call_cublas_header.replace( "complex_ptr", "void_ptr" );

                    info_map[ blas_routine ][ "argument_map" ][ arg ][ "code" ][ "call_cublas_header" ] = call_cublas_header

                else:
                    print "Could not match cublas argument '" + cublas_arg + "' to known arguments."
                    print arg, arguments
                    
                    exit(0)









