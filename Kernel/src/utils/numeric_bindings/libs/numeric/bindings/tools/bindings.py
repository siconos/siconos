#!/usr/bin/python
# -*- coding: utf-8 -*-
#
#  Copyright (c) 2008 Thomas Klimpel and Rutger ter Borg
#
# Distributed under the Boost Software License, Version 1.0.
# (See accompanying file LICENSE_1_0.txt or copy at
# http://www.boost.org/LICENSE_1_0.txt)
#

import re
import netlib

#
# Desired order of routines in C++ code: float, double, single complex, double complex
#
def routine_cmp( a, b ):
  letter_a = a[0]
  letter_b = b[0]
  value_map = { 'S': 0, 'D': 1, 'C': 2, 'Z': 3 }
  result = 0
  if value_map[ a[0] ] < value_map[ b[ 0] ]:
    result = -1
  if value_map[ a[0] ] > value_map[ b[ 0] ]:
    result = 1
  return result


#
# a < b?
# should replace the function "routine_cmp"
#
def subroutine_less( a, b, global_info_map ):

   if global_info_map[ a ][ 'value_type' ] < \
      global_info_map[ b ][ 'value_type' ]:
        return True
   if global_info_map[ a ][ 'value_type' ] > \
      global_info_map[ b ][ 'value_type' ]:
        return False

   if global_info_map[ a ][ 'precision' ] < \
      global_info_map[ b ][ 'precision' ]:
        return True
   if global_info_map[ a ][ 'precision' ] > \
      global_info_map[ b ][ 'precision' ]:
        return False

   return False


#
# This routine actually does what it's named after :-).
# Changes: 
# * used regex matching on delimiter instead of .find
#
def proper_indent( input_string ):
  max_chars = 80
  all_results = []
  find_delim = re.compile( "([\,\+\-/]|\|\||>=|< |std::log\(|work\()[ ]*" )
  for input_line in input_string.splitlines():
    result = ''
    # extra indentation size is 8
    indentation = len(input_line) - len(input_line.lstrip() ) + 8

    #print "indentation: ", indentation

    indentation_string = ''
    for i in range(indentation ):
      indentation_string += ' '
    #cur_index = input_line.find( ",", 0, -1 ) + 1

    match_delim = find_delim.search( input_line, 0, len(input_line) )
    if match_delim != None:
      #print match_delim.end( 0 )- match_delim.start( 0 )
      cur_index = match_delim.end( 0 )
    else:
      cur_index = 0

    prev_index = 0
    prev_slice = 0
    prev_indentation = 0
    while cur_index != 0:
      if ( (cur_index - prev_slice) > (79-prev_indentation) ):
        result = result + input_line[ prev_slice : prev_index ].rstrip() + '\n' + indentation_string
        prev_indentation = indentation
        prev_slice = prev_index
      prev_index = cur_index
      if cur_index != len(input_line):
        match_delim = find_delim.search( input_line, cur_index+1, len(input_line) )
        if match_delim != None:
          #print match_delim.end( 0 )- match_delim.start( 0 )
          cur_index = match_delim.end( 0 )
        else:
          cur_index = 0
        #cur_index = input_line.find( ",", cur_index+1, -1 ) + 1
        if cur_index == 0:
          cur_index = len(input_line)
      else:
        cur_index = 0
    cur_index = len( input_line )
    result = result + input_line[ prev_slice : cur_index ]
    all_results += [ result ]

  final_result =  "\n".join( all_results ) + "\n"
  final_result = final_result.replace( '(  )', '()' )

  return final_result



#
# For a group of routines, write their part in the 
# blas_names.h/lapack_names.h file.
#
def write_by_value_type( properties, template_map ):
  content = ''
  group_keys = properties.keys()
  group_keys.sort()
  for g in group_keys:
    content += '// Value-type variants of ' + g.lower() + '\n'
    for k in properties[ g ]:
      # avoid duplicate definitions for real matrices occuring in more than one group
      if (k[0] == 'S' or k[0] == 'D') and '_' not in g and g not in k:
          continue
      template = template_map[ template_map[ 'PARSERMODE' ].lower() + '_names.h_function' ]
      template = template.replace( '$SUBROUTINE', k )
      template = template.replace( '$subroutine', k.lower() )
      content += template
    content += '\n'
  return content

#
# Write the blas_names.h/lapack_names.h file.
#
def write_names_header( global_info_map, routines, template_map, dest_file ):
  parsermode = template_map[ 'PARSERMODE' ].lower()
  content = ''

  for level, level_properties in routines.iteritems():
    content += '//\n'
    content += '// ' + template_map[ 'PARSERMODE' ] + ' ' + level + ' routines\n'
    content += '//\n\n'

    if template_map[ 'PARSERMODE' ] == 'BLAS':
      if level_properties.has_key( 'routines_by_value_type' ):
        content += write_by_value_type( level_properties[ 'routines_by_value_type'], template_map )
    else:
      for problem_type, problem_properties in level_properties.iteritems():
        #
        # Here, we could write an informative header about the problem types
        #
        if problem_properties.has_key( 'routines_by_value_type' ):
          content += write_by_value_type( problem_properties[ 'routines_by_value_type'], template_map )

  result = template_map[ parsermode + '_names.h' ]
  result = result.replace( "$CONTENT", content )
  result = result.replace( "$LIBRARY_INT_TYPE", "fortran_int_t" )
  #//result = result.replace( "$PARSERMODE", template_map[ "PARSERMODE" ] )

  open( dest_file, 'w' ).write( result )






def write_header_part( global_info_map, properties, template_map ):
  parsermode = template_map[ 'PARSERMODE' ].lower()
  group_keys = properties.keys()
  group_keys.sort()
  content = ''
  for g in group_keys:
    content += '// Value-type variants of ' + g.lower() + '\n'
    for k in properties[ g ]:
      # avoid duplicate definitions for real matrices occuring in more than one group
      if (k[0] == 'S' or k[0] == 'D') and '_' not in g and g not in k:
          continue

      template = template_map[ parsermode + '.h_function' ]
      arg_list = []
      for arg in global_info_map[ k ][ 'arguments' ]:
        arg_list += [ global_info_map[ k ][ 'argument_map' ][ arg ][ 'code' ][ 'lapack_h' ] ]

      template = template.replace( "$SUBROUTINE", k )
      template = template.replace( "$ARGUMENTS", ", ".join( arg_list ) )
      template = template.replace( '$RESULT_TYPE', global_info_map[ k ][ 'return_value_type' ] )
      template = template.replace( '$RETURN_STATEMENT', global_info_map[ k ][ 'return_statement' ] )
      template = template.replace( "$LIBRARY_INT_TYPE", "fortran_int_t" )
      content += proper_indent( template )

    content += '\n'
  return content



#
# Write the blas.h/lapack.h file
#
def write_header( global_info_map, routines, template_map, dest_file ):
  parsermode = template_map[ 'PARSERMODE' ].lower()
  content = ''
  
  for level, level_properties in routines.iteritems():
    content += '//\n'
    content += '// ' + template_map[ 'PARSERMODE' ] + ' ' + level + ' routines\n'
    content += '//\n\n'

    if template_map[ 'PARSERMODE' ] == 'BLAS':
      if level_properties.has_key( 'routines_by_value_type' ):
        content += write_header_part( global_info_map, \
          level_properties[ 'routines_by_value_type'], template_map )
    else:
      for problem_type, problem_properties in level_properties.iteritems():
        #
        # Here, we could write an informative header about the problem types
        #
        if problem_properties.has_key( 'routines_by_value_type' ):
          content += write_header_part( global_info_map, \
                problem_properties[ 'routines_by_value_type'], template_map )

  result = template_map[ parsermode + '.h' ]
  result = result.replace( "$CONTENT", content )
  result = result.replace( "$LIBRARY_INT_TYPE", "fortran_int_t" )
  #result = result.replace( "$PARSERMODE", template_map[ "PARSERMODE" ] )

  open( dest_file, 'w' ).write( result )



#
# Write the include hierarchy header(s)
#
def write_include_hierarchy( global_info_map, routines, template_map, dest_path ):

  parsermode = template_map[ 'PARSERMODE' ].lower()

  for level, level_properties in routines.iteritems():
    content = ''
    postfix = '.hpp'
    if template_map[ 'PARSERMODE' ] == 'LAPACK_DOC':
        postfix = '.qbk'
    dest_file = dest_path + '/' + level + postfix

    if template_map[ 'PARSERMODE' ] == 'BLAS':
      print "something"

    else:
      # problem type = general_eigen, etc.
      # problem properties = the mapped stuff
      group_keys = []
      for problem_type, problem_properties in level_properties.iteritems():
        # the key routines_by_value_type usually has 4 routines for a group of routines
        if problem_properties.has_key( 'routines_by_value_type' ):
            group_keys.extend( problem_properties[ 'routines_by_value_type' ].keys() )
      group_keys.sort()

      prev_key = ""
      for r in group_keys:
          if r == prev_key:
              continue
          prev_key = r
          if template_map[ 'PARSERMODE' ] == 'LAPACK':
              content += '#include <boost/numeric/bindings/' + parsermode + '/' + level + \
                  '/' + r.lower() + '.hpp>\n'
          if template_map[ 'PARSERMODE' ] == 'LAPACK_DOC':
              content += '[include ' + level + \
                  '/' + r.lower() + '.qbk]\n'

    result = template_map[ parsermode + '_include_hierarchy' ]
    result = result.replace( "$CONTENT", content )
    result = result.replace( "$LEVEL", level.upper() )
    result = result.replace( "$CAPITALIZED_LEVEL", level.capitalize() )

    open( dest_file, 'w' ).write( result )


#
# Generate const-overloads
#

def generate_const_variants( prefix, argument_list, template_map ):
    print "Const variant generation prefix: ", prefix
    print "Generating const variants for ", argument_list
    permute_indices = []
    result = []
    comments = []

    # use this arg list may contain a modified version of the
    # argument list, where in some overruled cases the 'const '
    # will be removed from the argument (force permutation)
    use_this_arg_list = []

    for i in range( 0, len(argument_list) ):
        argument = argument_list[i]
        if 'const' not in argument and \
           'char' not in argument and \
           'ptrdiff_t' not in argument and \
           'typename' not in argument and \
           '$LIBRARY_INT_TYPE' not in argument and \
           'fortran_bool_t' not in argument and \
           'Vector' not in argument and \
           'Matrix' not in argument and \
           '&' in argument:
            permute_indices.append( i )
        else:
            arg_name = argument.split( ' ' )[ -1 ].upper()
            my_key = prefix + '.' + arg_name + '.level2_permute'
            if netlib.my_has_key( my_key, template_map ):
                permute_indices.append( i )
                argument = argument.replace( 'const ', '' )
        use_this_arg_list.append( argument )

    print " To be permuted: ", permute_indices

    for i in range( 0, pow( 2, len( permute_indices ) ) ):
        #print "i: ", i
        new_arg_list = []
        new_arg_list += use_this_arg_list
        new_comments = []
        for j in range( 0, len( permute_indices ) ):
            if ( i & (1<<j) ):
                #print permute_indices[j], ": const " + use_this_arg_list[ permute_indices[ j ] ]
                new_arg_list[ permute_indices[ j ] ] = "const " + use_this_arg_list[ permute_indices[ j ] ]
                arg = new_arg_list[ permute_indices[ j ] ]
                new_comments.append( "// * " + 
                    arg[ :arg.find("&" ) ] + "&" )
            else:
                arg = new_arg_list[ permute_indices[ j ] ]
                new_comments.append( "// * " + 
                    arg[ :arg.find("&" ) ] + "&" )

           # else:
                #print permute_indices[j], "don't add const"
        if new_arg_list not in result:
            result.append( new_arg_list )
            comments.append( new_comments )
        #new_arg_list = []

    #print "result: ", result
    return result, comments

#
# Generate the call to the value_tye meta-func
#
def value_type( arg ):
    return 'typename $NAMESPACEvalue_type< ' + arg + ' >::type'

#
# Search replace stuff for handling exceptional cases through the
# templating system
#

def search_replace( source, key_name, template_map ):
    print "Trying key ", key_name
    if key_name not in template_map:
        return source
    result = source
    for search_replace in template_map[ key_name ].split( "--" ):
        if "->" in search_replace:
            splitted = search_replace.split("->")
            print "Replacing '" + splitted[0] + "' with '" + splitted[1] + "'"
            result = result.replace( splitted[0], splitted[1] )
    return result

