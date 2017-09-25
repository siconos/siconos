
__all__ = ['unwanted', 'get_priority', 'parse_args', 'get_headers',
           'write_header', 'write_footer', 'write_includes',
           'write_register_with_bases', 'write_classes']

import os
import os.path
import sys
import re
import itertools
import getopt

input_headers = {
    'kernel': ["SiconosKernel.hpp"],

    'mechanics': ["MechanicsFwd.hpp", "SpaceFilter.hpp",
                  "SpaceFilter_impl.hpp",
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
                  "KneeJointR.hpp", "PivotJointR.hpp",
                  "PrismaticJointR.hpp",
                  "FixedJointR.hpp",
                  "CylindricalJointR.hpp",
                  "NewtonEulerJointR.hpp",
                  "JointStopR.hpp",
                  "JointFrictionR.hpp",
                  "BodyDS.hpp", "SiconosShape.hpp", "ContactR.hpp",
                  "SiconosCollisionQueryResult.hpp",
                  "SiconosContactor.hpp",
                  "SiconosCollisionManager.hpp"],

    # fix missing forwards for Control
    'control': ['FirstOrderNonLinearDS.hpp',
                'FirstOrderLinearDS.hpp',
                'SiconosControl.hpp'],
    }


def unwanted(s):
    """ un processed classed or attributes : to be defined explicitely in SiconosFull.hpp"""
    m = re.search('xml|XML|Xml|MBlockCSR|fPtr|SimpleMatrix|SiconosVector|SiconosGraph|SiconosSharedLibrary|numerics|computeFIntPtr|computeJacobianFIntqPtr|computeJacobianFIntqDotPtr|PrimalFrictionContact|FrictionContact|Lsodar|_moving_plans|_err|Hem5|_bufferY|_spo|_measuredPert|_predictedPert|_blockCSR', s)
    # note _err,_bufferY, _spo, _measuredPert, _predictedPert -> boost::circular_buffer issue with serialization
    # _spo : subpluggedobject
    # _blockCSR -> double * serialization needed by hand (but uneeded anyway for a full restart)
    return m is not None

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
def get_priority(name, source_dir, header_path, header_line):
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
    prio = header_line/10000.
    for e in sorted(module_prio, key=lambda k: k[1]):
        if e[0] in header_path:
            prio += e[1]
            if e[0] == 'kernel':
                path, name = resolve_path(source_dir, header_path)
                for ee in sorted(kernel_prio, key=lambda k: k[1]):
                    if ee[0] in path:
                        prio += ee[1]
                        return prio
            else:
                return prio
    print('Error proccessing header {:}'.format(header_path))


def resolve_path(source_dir, header_path):
    """Find a header file by name under the given source path."""
    header_file = os.path.basename(header_path)
    # print('walking {:}, looking for {:}'.format(source_dir + '/kernel/src', header_file))
    for path, dirlist, filelist in (
            os.walk(os.path.join(source_dir, 'kernel/src'))):
        files = [i for i in filelist if i == header_file]
        if len(files) == 1:
            return path, files[0]
        elif len(files) >1:
            raise ValueError('Multiple matches found for '+header_file)
    raise ValueError('Error while finding {:}, no match found.'
                     .format(header_file))


def usage():
    myname = sys.argv[0]
    print(' '.join("""{0} [--namespace=<namespace>] -I<path> [-I<path> ...]
    [--targets=<Mod1>[,Mod2[,...]]]
    [--output=<filename>]
    [--source=<siconos source dir>]
    [--build=<siconos build dir>]
    header""".format(myname).split()))

def parse_args(need_build_path=False):
    """Parse command-line arguments for generated header builder utility."""
    include_paths = []
    siconos_namespace = '::'
    myname = sys.argv[0]

    try:
        longopts = ['help', 'namespace=', 'targets=', 'output=', 'source=']
        if need_build_path:
            longopts.append('build=')
        opts, args = getopt.getopt(sys.argv[1:], 'I:', longopts)
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        sys.exit(2)

    targets = []
    generated_file = None
    source_dir = None
    build_path = None

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
        if opt == '--build':
            build_path = arg

    if generated_file is None:
        usage()
        print('{0} --output option is mandatory.'.format(myname))
        sys.exit(1)

    if source_dir is None:
        usage()
        print('{0} --source  option is mandatory.'.format(myname))
        sys.exit(1)

    if need_build_path and build_path is None:
        usage()
        print('{0} --build  option is mandatory.'.format(myname))
        sys.exit(1)

    generated_header = os.path.splitext(os.path.basename(generated_file))[0]

    if len(args) != 0:
        usage()
        sys.exit(1)

    return (include_paths,
            siconos_namespace,
            targets,
            generated_file,
            source_dir,
            generated_header,
            build_path)


def get_headers(targets):
    # ensure a certain ordering between these three targets
    order = ['mechanics','kernel','control']
    ordered_targets = []
    for o in order:
        if o in targets:
            ordered_targets.append(o)
    # they must come after all other targets
    ordered_targets = [t for t in targets if t not in order] + ordered_targets
    all_headers = [h for h in itertools.chain(*(input_headers[target]
                                                for target in ordered_targets))]
    return all_headers


def write_header(dest_file, cmd, generated_header):
    dest_file.write('// generated with {0}\n'.format(os.path.split(sys.argv[0])[-1]))
    dest_file.write('#ifndef {0}_hpp\n'.format(generated_header))
    dest_file.write('#define {0}_hpp\n'.format(generated_header))
    dest_file.write('#include <SiconosConfig.h>\n'.format(generated_header))
    dest_file.write('#ifdef WITH_SERIALIZATION\n'.format(generated_header))


def write_footer(dest_file):
    dest_file.write('#endif\n')
    dest_file.write('#endif\n')


def write_includes(dest_file, all_headers):
    for header in all_headers:
        dest_file.write('#include "{0}"\n'.format(header))


def write_register_with_bases(dest_file, with_base):
    with_base_s = [c for c, p in sorted(with_base, key=lambda k: (k[1], k[0]))]
    dest_file.write('\n')
    dest_file.write('template <class Archive>\n')
    dest_file.write('void siconos_io_register_generated(Archive& ar)\n')
    dest_file.write('{{\n{0}\n}}\n'
                    .format('\n'
                            .join(
                                '  ar.register_type(static_cast<{0}*>(NULL));'
                                .format(x) for x in with_base_s)))


def write_classes(dest_file, classes):
    """Write registration macros for classes, according to whether or
       not they have serializable bases.

       Input is a list of tuples of (class name, list of
       serializable bases, list of wanted members)."""
    for clname, bases, members in classes:

        # Write classes according to whether they have serializable bases
        if len(bases) > 0:
            dest_file.write(
                'SICONOS_IO_REGISTER_WITH_BASES({0},{1},\n'
                .format(clname, '('+')('.join(bases)+')'))
        else:
            dest_file.write('SICONOS_IO_REGISTER({0},\n'
                            .format(clname))

        # Write (wanted) member variables
        dest_file.write('\n'.join(sorted('  ({0})'.format(m) for m in members
                                         if not unwanted(m))) + ')\n')
