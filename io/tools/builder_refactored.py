#!/usr/bin/env python2

# ./builder.py -I/usr/local/include/siconos

# mechanics
# ./builder.py -I/usr/local/include/siconos --target=mechanics

# control
# ./builder.py -I/usr/local/include/siconos --target=control

# we use pygccxml from Roman Yakovenko.
# http://sourceforge.net/projects/pygccxml/

import sys
import os
import getopt
import re
import itertools
from pygccxml import parser
from pygccxml import declarations
from pygccxml.parser import COMPILATION_MODE

from builder_common import *


(include_paths,
 siconos_namespace,
 targets,
 generated_file,
 source_dir,
 generated_header,
 build_path) = parse_args()

all_headers = get_headers(targets)


def is_serializable(something):
    return 'serializable' in [_c_.name for _c_ in
                              something.typedefs(allow_empty=True)]


def name(t):
    if isinstance(t, declarations.class_t):
        return t.name
    elif isinstance(t, declarations.typedef_t):
        return t.decl_string[2:]  # remove ::


def replace_by_typedef(some_type):
    if str(some_type) in typedef:
        rep_typedef = typedef[str(some_type)]
        if not '<' in rep_typedef:  # replace only if not a template
            return rep_typedef
    return str(some_type)


# main loop
if 'xml_generator_configuration_t' in dir(parser):
    config = parser.xml_generator_configuration_t(include_paths=include_paths,
                                                  ignore_gccxml_output=False,
                                                  keep_xml=True)
else:
    config = parser.config_t(include_paths=include_paths, ignore_gccxml_output=True)

decls = parser.parse(all_headers, config,  compilation_mode=COMPILATION_MODE.ALL_AT_ONCE)
global_ns = declarations.get_global_namespace(decls)

# classes in siconos_namespace
class_names = dict()

# class name of classes with a least a base (for the boost archive
# registration)
with_base = []

# a typedef table to replace templated class by their typedefs in
# macros call
typedef = dict()
for t in global_ns.typedefs():
    typedef[str(t._type)] = name(t)

with open(generated_file, 'a') as dest_file:
    write_header(dest_file, ' '.join(sys.argv), generated_header)
    write_includes(dest_file, all_headers)

    for type_ in filter(lambda c: c.parent.name == siconos_namespace,
                        itertools.chain(
                            global_ns.classes(), global_ns.typedefs())):

        is_typedef = False

        if isinstance(type_, declarations.class_t):
            class_names[declarations.full_name(type_)] = type_
            class_ = type_
        elif isinstance(type_, declarations.typedef_t):
            try:
                is_typedef = True
                class_ = class_names['::' + str(type_.type)]
            except:
                class_ = None
        # with the serializable tag
        # (could not find friend functions with pygccxml)
        if class_ is not None and \
            is_serializable(class_) and \
            (is_typedef or not
             declarations.templates.is_instantiation(class_.name)):
            if not unwanted(class_.name):
                if not class_.is_abstract:
                    with_base.append(
                        (class_.name,
                         get_priority(class_.name, source_dir,
                                      type_.location.file_name,
                                      type_.location.line)))

                # print registration macros depending on inheritance
                if class_.bases == []:
                    dest_file.write(
                        'SICONOS_IO_REGISTER({0},\n'.format(name(type_)))
                else:
                    serializable_bases = \
                        reduce(lambda r, b:
                               r + [b.related_class]
                               if is_serializable(b.related_class)
                               and
                               b.related_class.parent.name == siconos_namespace
                               else r, class_.bases, [])
                    if len(serializable_bases) > 0:
                        dest_file.write(
                            'SICONOS_IO_REGISTER_WITH_BASES({0},{1},\n'
                            .format(name(type_), ''.join(['({0})'
                                .format(replace_by_typedef(c.name))
                                for c in serializable_bases])))
                    else:
                        dest_file.write('SICONOS_IO_REGISTER({0},\n'
                                        .format(name(type_)))

                variables = [v.name
                             for v in filter(lambda x: not 'void'
                                             in str(x._get_type()),
                                             class_.variables(
                                                 allow_empty=True))]
                dest_file.write('{0})\n'
                                .format('\n'
                                        .join('  ({0})'
                                              .format(vn)
                                              for vn in
                                              filter(lambda x: not unwanted(x),
                                                     variables))))

    # filtering is not correct at this point
    # some unwanted classes are necessary
    # (the ones in SiconosFull.hpp) others not (xml)
    # some leads to compilation errors.
    write_register_with_bases(dest_file, with_base)

    write_footer(dest_file)
