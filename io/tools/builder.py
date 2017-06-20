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

myname = sys.argv[0]
include_paths = []
siconos_namespace = '::'

def usage():
    print('{0} [--namespace=<namespace>] -I<path> [-I<path> ...] \
               [--targets=<Mod1>[,Mod2[,...]]] \
               [--output=<filename>] \
               [--source=<siconos source dir>] \
               header'.format(myname))


try:
    opts, args = getopt.getopt(sys.argv[1:], 'I:', ['help', 'namespace=',
                                                    'targets=', 'output=', 'source='])
except getopt.GetoptError as err:
    print(str(err))
    usage()
    sys.exit(2)

targets = []
generated_file = None
source_dir = None

for opt, arg in opts:
    if opt == '--targets':
        targets = arg.split(',')
    if opt == '-I':
        include_paths += [arg]
    if opt == '--namespace':
        siconos_namespace = arg
    if opt == '--help':
        usage()
        sys.exit(0)
    if opt == '--output':
        generated_file = arg
    if opt == '--source':
        source_dir = arg

if generated_file is None:
    usage()
    print('{0} --output option is mandatory.'.format(myname))
    sys.exit(1)

if source_dir is None:
    usage()
    print('{0} --source  option is mandatory.'.format(myname))
    sys.exit(1)

generated_header = os.path.splitext(os.path.basename(generated_file))[0]

if len(args) != 0:
    usage()
    sys.exit(1)

input_headers = dict()

input_headers['kernel'] = ["SiconosKernel.hpp"]

input_headers['mechanics'] = ["MechanicsFwd.hpp", "SpaceFilter.hpp", "SpaceFilter_impl.hpp",
                             "ExternalBody.hpp",
                             "Disk.hpp", "Circle.hpp", "DiskDiskR.hpp",
                             "DiskMovingPlanR.hpp",
                             "DiskPlanR.hpp", "SphereLDS.hpp",
                             "SphereLDSPlanR.hpp",
                             "SphereLDSSphereLDSR.hpp", "SphereNEDS.hpp",
                             "SphereNEDSPlanR.hpp",
                             "SphereNEDSSphereNEDSR.hpp",
                             "SiconosBodies.hpp",
                             "CircleCircleR.hpp", "CircularDS.hpp",
                              "KneeJointR.hpp", "PivotJointR.hpp", "PrismaticJointR.hpp"
                             ]

# fix missing forwards for Control
input_headers['control'] = ['FirstOrderNonLinearDS.hpp', 'FirstOrderLinearDS.hpp', 'SiconosControl.hpp']

all_headers = [h for h in itertools.chain(*(input_headers[target]
                                            for target in targets))]


def is_serializable(something):
    return 'serializable' in [_c_.name for _c_ in
                              something.typedefs(allow_empty=True)]


# un processed classed or attributes : to be defined explicitely in
# SiconosFull.hpp
def unwanted(s):
    m = re.search('xml|XML|Xml|MBlockCSR|fPtr|SimpleMatrix|SiconosVector|SiconosSet|SiconosGraph|SiconosSharedLibrary|numerics|computeFIntPtr|computeJacobianFIntqPtr|computeJacobianFIntqDotPtr|PrimalFrictionContact|FrictionContact|Lsodar|_moving_plans|_err|Hem5|_bufferY|_spo|_measuredPert|_predictedPert|_blockCSR', s)
    # note _err,_bufferY, _spo, _measuredPert, _predictedPert -> boost::circular_buffer issue with serialization
    # _spo : subpluggedobject
    #_blockCSR -> double * serialization needed by hand (but uneeded anyway for a full restart)
    return m is not None


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

# try to provide an ordering for registering a class
# The main issue is with the Model and the NonSmoothDynamicalSystem, Topology,
# DSG, ...
# If this is not done, the compiler may fail (like Clang) because of too many
# recursive template instantiations
# The logic is the following one:
# - give priority based on the module
# - for the Kernel give priority based on the folder
# - some hacks are needed to get the Topology & co done properly
# the returned priority is used to sort the list of classes to register
def get_priority(type_, name):
    module_prio = (('numerics', 100),
                   ('kernel', 200),
                   ('mechanics', 300),
                   ('control', 400))

    kernel_prio = (('utils/SiconosException', 0),
                   ('utils/Memory', 1),
                   ('utils/SiconosAlgebra', 2),
                   ('utils/SiconosTools', 3),
                   ('utils', 4),
                   ('plugin', 5),
                   ('modelingTools', 6),
                   ('simulationTools', 7),
                   ('model', 8),
                   (r'.*', 9))

    big_hack_prio = {'GraphProperties': 1e-3,
                     'DynamicalSystemProperties': 2e-3,
                     'InteractionProperties': 3e-3,
                     'MatrixIntegrator': 4e-3,
                     'DynamicalSystemsGraph': 5e-3,
                     'InteractionsGraph': 6e-3,
                     'Topology': 7e-3}

    if name in big_hack_prio.keys():
        return 200 + 5 + big_hack_prio[name]
    header = type_.location.file_name
    prio = type_.location.line/10000.
    for e in sorted(module_prio, key=lambda k: k[1]):
        if e[0] in header:
            prio += e[1]
            if e[0] is 'kernel':
                header_file = os.path.basename(header)
                # print('walking {:}, looking for {:}'.format(source_dir + '/kernel/src', header_file))
                for path, dirlist, filelist in os.walk(source_dir + '/kernel/src'):
                    files = [i for i in filelist if i == header_file]
                    if len(files) == 1:
                        # print('Found {:} in {:}'.format(header_file, path))
                        for ee in sorted(kernel_prio, key=lambda k: k[1]):
                            if ee[0] in path:
                                prio += ee[1]
                                return prio
                    elif len(files) >1:
                        print(files)
                print('Error while finding {:}, found no match'.format(header_file))
            else:
                return prio
    print('Error proccessing header {:}'.format(header))



# main loop
if 'xml_generator_configuration_t' in dir(parser):
    config = parser.xml_generator_configuration_t(include_paths=include_paths, ignore_gccxml_output=False,
                                                  keep_xml=True, xml_generator='castxml')
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

    dest_file.write('// generated with the command : {0}\n'
                    .format(' '.join(sys.argv)))
    dest_file.write('#ifndef {0}_hpp\n'.format(generated_header))
    dest_file.write('#define {0}_hpp\n'.format(generated_header))
    dest_file.write('#include <SiconosConfig.h>\n'.format(generated_header))
    dest_file.write('#ifdef WITH_SERIALIZATION\n'.format(generated_header))

    for header in all_headers:
        dest_file.write('#include "{0}"\n'.format(header))

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
                    with_base.append((class_.name, get_priority(type_, class_.name)))

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
    with_base_s = zip(*sorted(with_base, key=lambda k: k[1]))[0]
    dest_file.write('\n')
    dest_file.write('template <class Archive>\n')
    dest_file.write('void siconos_io_register_generated(Archive& ar)\n')
    dest_file.write('{{\n{0}\n}}\n'
                    .format('\n'
                            .join(
                                '  ar.register_type(static_cast<{0}*>(NULL));'
                                .format(x) for x in with_base_s)))
    dest_file.write('#endif\n')
    dest_file.write('#endif\n')
