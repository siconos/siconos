#!/usr/bin/env python3

import os, os.path, sys
from glob import glob
import lxml.etree
import re

from builder_common import *

def get_classes_conditional(doxy_xml_files, cond):
    """Get classes and members from a list of Doxygen XML files that
       meet the given condition."""
    found = {}
    for xmlfile in doxy_xml_files:
        xml = lxml.etree.parse(xmlfile)
        classes = xml.xpath('.//compounddef[@kind="class" or @kind="struct"]')
        for cl in classes:
            if cond(cl):
                classname = cl.find('./compoundname')
                baseclasses = cl.xpath('./basecompoundref')
                membervars = cl.xpath('.//memberdef[@kind="variable"]/name')

                # An exception: Members get attached to Graph classes
                # through this macro, and is not understood by
                # Doxygen, so we have to parse it outselves.
                graphvars = cl.xpath('.//memberdef[@kind="function"]/name'
                                     +'[text()="INSTALL_GRAPH_PROPERTIES"]')
                graphmems = []

                if len(graphvars)>0:
                    r = re.compile('\(\(\w+,\s*[\w: ]+,\s*(\w+)\)\)')
                    for g in graphvars:
                        for a in g.xpath('../argsstring'):
                            graphmems += r.findall(a.text)
                    # The INSTALL_GRAPH_PROPERTIES macro also adds a
                    # bool called "dummy"
                    graphmems.append('dummy')

                location = cl.find('./location')
                found[classname.text] = (
                    {'name': classname.text,
                     'bases': [base.text for base in baseclasses],
                     'members': [mem.text for mem in membervars] + graphmems,
                     'filepath': location.attrib['file'],
                     'line': int(location.attrib['line']),
                     'abstract': cl.xpath('@abstract="yes"'),
                    })
    return found

def classes_from_build_path(build_path, targets):
    """Get classes and members from all Doxygen XML files found on the
       provided build path."""
    doxy_xml_path = os.path.join(build_path,'docs/build/doxygen/doxy2swig-xml')
    if not os.path.exists(doxy_xml_path):
        print('%s: Error, path "%s" does not exist.'%(sys.argv[0], doxy_xml_path))
        sys.exit(1)
    doxy_xml_files = [p for t in targets
                      for p in (glob(os.path.join(doxy_xml_path, '%s/class*.xml'%t))
                                + glob(os.path.join(doxy_xml_path, '%s/struct*.xml'%t)))]

    # We want only classes that contain calls to the
    # ACCEPT_SERIALIZATION macro.
    serializable = './/memberdef/name[text()="ACCEPT_SERIALIZATION"]'

    def pred(x):
        return len(x.xpath(serializable)) > 0

    return get_classes_conditional(doxy_xml_files, pred)

def assign_targets(classes, source_dir):
    """For each class, figure out which target it is associated with based
       on the file path."""
    for cl in classes.values():
        cl['target'] = get_target(source_dir, cl['filepath'])

def assign_priorities(classes, source_dir):
    """For each class, get its priority to help in ordering the
       declarations in the generated file."""
    for cl in classes.values():
        cl['priority'] = get_priority(cl['name'], source_dir,
                                      cl['filepath'], cl['line'])

def resolve_base_classes(classes):
    """For each class, find those base classes which are also in the list.
       This is to exclude STL and BOOST base classes,
       e.g. enable_shared_from_this<>.  In practice we only want to
       list base classes which are serializable."""
    for cl in classes.values():
        resolved = []
        for base in cl['bases']:
            if base in classes:
                resolved.append(base)
        cl['resolved_bases'] = resolved

def remove_unwanted_resolved(classes):
    """This is a bit of an ugly hack: For some classes, they are not in
       "unwanted" because we want them to resolve as base classes, but
       they are covered in SiconosFull.hpp, so we don't want them to
       appear in generated headers."""
    unwanted_resolved = ['_DynamicalSystemsGraph', '_InteractionsGraph']
    for u in unwanted_resolved:
        if u in classes:
            del classes[u]

def classes_from_headers(all_headers, include_paths):
    """Use compiler preprocessor to find an approximate list of classes
       referenced by a set of headers.  May return some words which
       are not classes."""
    import os, os.path, tempfile, shutil
    classes = []
    try:
        d = tempfile.mkdtemp()
        hpp = os.path.join(d, 'headers.hpp')
        compiled = os.path.join(d, 'out.cpp')
        with open(hpp,'w') as h:
            [print('#include "%s"'%i, file=h) for i in all_headers]
        cxx = 'g++'
        if 'CXX' in os.environ:
            cxx = os.environ['CXX']
        cmd = [cxx, '-o', compiled, '-E']
        for i in include_paths:
            cmd += ['-I', i]
        cmd.append(hpp)
        # print(' '.join(cmd))
        os.system(' '.join(cmd))
        with open(compiled, 'r') as out:
            for line in out:
                words = line.split()
                if len(words)>=2 and (words[0]=='class' or words[0]=='struct'):
                    classes.append(words[1])
        return classes
    finally:
        shutil.rmtree(d)

if __name__=='__main__':
    (include_paths,
     siconos_namespace,
     targets,
     generated_file,
     source_dir,
     generated_header,
     build_path) = parse_args(need_build_path=True)

    all_headers = get_headers(targets)
    doxygen_classes = classes_from_build_path(build_path, targets)

    header_classes = classes_from_headers(all_headers, include_paths)

    # Allow for inner classes
    def in_maybe_inner(k,h):
        if '::' in k:
            # (a bit loose but we can't forward-declare them so match
            # individual elements instead)
            return all([c in h for c in k.split('::')])
        return k in h

    # Find the join of the two lists
    classes = {k: v for k,v in doxygen_classes.items()
               if in_maybe_inner(k, header_classes) and not unwanted(k)}

    print('{:} classes found for targets {:}.'.format(len(classes), ', '.join(targets)))

    if len(classes) < 10:
        print('%s: Error, not enough classes found.'%sys.argv[0])
        sys.exit(1)

    assign_targets(classes, source_dir)

    assign_priorities(classes, source_dir)

    resolve_base_classes(classes)

    remove_unwanted_resolved(classes)

    with open(generated_file, 'w') as dest_file:
        write_header(dest_file, ' '.join(sys.argv), generated_header)
        write_includes(dest_file, all_headers)

        sorted_classes = sorted(classes.values(),
                                key = lambda k: (k['priority'], k['name']))

        class_list = [
            (cl['name'],
             cl['resolved_bases'],
             [m for m in cl['members'] if not unwanted(m)]
             )
            for cl in sorted_classes]

        write_classes(dest_file, class_list)

        with_base = [(cl['name'],
                      cl['priority'],
                      cl['target'])
                     for cl in sorted_classes
                     if not cl['abstract']]

        # Note: This was the condition before, but noticed that
        # builder.py did not check the number of bases before adding
        # to with_base!
        # if len(cl['bases'])>0 and not cl['abstract']]

        # filtering is not correct at this point
        # some unwanted classes are necessary
        # (the ones in SiconosFull.hpp) others not (xml)
        # some leads to compilation errors.
        write_register_with_bases(dest_file, with_base)

        write_footer(dest_file)
