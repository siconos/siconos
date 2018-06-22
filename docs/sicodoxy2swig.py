"""Tools to generate docstrings for siconos python interface.

Process:

#1. doxygen generates xml outputs from comments in C/C++ headers
#2. SiconosDoxy2Swig class (derived from doxy2swig https://github.com/m7thon/doxy2swig) is
used to generates .i files (swig docstring feature) from xml
#3. .i files are used by swig to create docstrings in .py

This file is to be copied into CMAKE_BINARY_DIR/share using configure_file
"""

""" Siconos is a program dedicated to modeling, simulation and control
 of non smooth dynamical systems.

 Copyright 2018 INRIA.

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
"""

from breathe import node_factory as nf
from itertools import count
import pickle
import os, sys
import textwrap

try:
    from doxy2swig import Doxy2SWIG
except:
    pass # 



type_map = {
    'SiconosMatrix': 'array_like (np.float64, 2D)',
    'SimpleMatrix': 'array_like (np.float64, 2D)',
    'SiconosVector': 'array_like (np.float64, 1D)',
    'SP::SiconosMatrix': 'array_like (np.float64, 2D)',
    'SP::SimpleMatrix': 'array_like (np.float64, 2D)',
    'SP::SiconosVector': 'array_like (np.float64, 1D)',
    r'double*': 'array_like (np.float64, 1D)',
    r'double *': 'array_like (np.float64, 1D)',
    'Index': 'array_like (int, 1D)',
    'const' : '',
    r'&' : '',
    'unsigned': '',
    r'std::': '',
    r'SP::': '',
    'friend': '',
    'string': 'str',
    'void': 'None',
    'Siconos::': '',
    'UBLAS_TYPE': 'int',
    'false': 'False',
    'true' : 'True',
    }

class SiconosDoxy2Swig(Doxy2SWIG):
    """Converts Doxygen generated XML files into a file containing
    docstrings that can be used by swig (feature("docstring")).
    Once the data is parsed it is stored in
    self.pieces.
    """

    """Global (static var!) counter for latex formulas in doctrings/swig comments.
    """
    __ids = count(0)

    def __init__(self, xmlfile, component, swig_workdir):

        
        with_function_signature = True # Replace 'c' like arg list in function
        # with something more pythonic
        with_type_info = True
        with_constructor_list = True
        with_attribute_list = True
        with_overloaded_functions = True
        textwidth = 80
        quiet = False

        super().__init__(xmlfile,
                         with_function_signature,
                         with_type_info, with_constructor_list,
                         with_attribute_list,
                         with_overloaded_functions,
                         textwidth, quiet)

        # dict to export latex formula, to be processed later
        self.latex_forms = {} 
        # dict to export enums, to be processed later
        self.enums = {} 
        # dict to trace list of all features, used as input for autodoc
        # to sort functions doc by file of origin
        self.features = []
        # xml source file name
        self.xml_filename = xmlfile

        self.name = os.path.basename(self.xml_filename).split('.')[0]
        # name (full path) of the .i file that will be created
        self.swig_outputname = os.path.join(swig_workdir, self.name + '.i')
        # pickle file used to save latex forms
        self.swig_outputname = os.path.join(swig_workdir, self.name + '.i')
        self.latex_filename = os.path.join(swig_workdir, 'latex_' + self.name + '.pickle')

        self.ignores.append('author')
        self.ignores.append('date')
        self.ignores.append('version')

    def add_separator(self):
        """add horizontal line in doc
        """
        self.add_text(['\n', 15 * '-', '\n'])
        
    def parse_typemap(self, source):
        """Replace C types with python types in docstrings (func. prototypes)
        """
        for item in type_map:
            source = source.replace(item, type_map[item])
        return source
    
    def parse_enum(self, node):
        """Parse memberdef node with kind=enum
        """
        for n in self.get_specific_subnodes(node, 'enumvalue', recursive=4):
            # Get name of the enum
            ename = self.extract_text(self.get_specific_subnodes(n, 'name'))
            # Value of the enum
            evalue = self.extract_text(self.get_specific_subnodes(n, 'initializer'))
            # Description
            edescr = self.extract_text(self.get_specific_subnodes(n, 'briefdescription'))
            self.enums[ename] = (evalue, edescr)

    def parse_typedefs(self, node):
        """Parse memberdef node with kind=typedef
        """
        # Get name of the typedef
        ename = self.extract_text(self.get_specific_subnodes(node, 'name'))
        # type
        old_type = self.extract_text(self.get_specific_subnodes(node, 'type'))
        # Description
        edescr = self.extract_text(self.get_specific_subnodes(node, 'briefdescription'))
        if len(edescr) > 0:
            self.enums[ename] = (old_type, edescr)
            
    def write(self, fname):
        with open(fname, 'w') as o:
            if(sys.version_info > (3, 0)):
                o.write(''.join(self.pieces))
            else:
                o.write(str(u''.join(self.pieces).encode('utf-8').strip()))
            o.write('\n')

    def do_includes(self, node):
        # Get header file path
        locnode = self.get_specific_subnodes(self.xmldoc, 'location', recursive=3)[0]
        headername = locnode.attributes['file'].value
        refname = headername.replace('/', '_')
        self.add_text(['\n', 'Generated class (swig), based on C++ header ',
                       ':ref:`pgm_' + refname + '`.', '\n\n'])
        # self.subnode_parse(node)
        
    def do_formula(self, node):
        """Read latex formula from xml file
        and post-process it for swig/sphinx-apidoc.
        
        Parameters
        ----------
        content : xml node
    
        The current formula is appended to latex_forms.
        
        Note
        ----
        * latex formula in doxygen comments leads to a real mess 
        during xml->doxy2swig-->swig ... process.
        
        """    
        # ret = doctools.xml_formula_to_swig(node.firstChild.data.strip())
        # self.add_text(' ')
        # self.add_text(ret[0])

        content = node.firstChild.data.strip()
        __id = next(self.__ids)
        node_factory = nf.create_node_factory()
        rst_node = None
        latex = content
        # Either inline
        if latex.startswith("$") and latex.endswith("$"):
            latex = latex[1:-1]
            latex = r':math:`' + latex.strip() + r'`'
            # If we're inline create a math node like the :math: role
            rst_node = node_factory.math()
        else:
            # Else we're multiline
            rst_node = node_factory.displaymath()

        # Or multiline
        if latex.startswith("\[") and latex.endswith("\]"):
            latex = latex[2:-2:]

        # Here we steal the core of the mathbase "math" directive handling code from:
        #    sphinx.ext.mathbase
        rst_node["latex"] = latex

        # Required parameters which we don't have values for
        rst_node["label"] = 'inline'
        rst_node["nowrap"] = False
        rst_node["docname"] = None
        rst_node["number"] = None
        self.latex_forms[__id] = rst_node
        # !!  '_' is required 
        # e.g. replace('FORMULA1') leads to wrong behavior because of FORMULA10
        # while 'FORMULA1_' works.
        id_formula = 'FORMULA' + str(__id) + '_'
        #self.add_text(' ')
        self.add_text(id_formula)

    def do_verbatim(self, node):
        """Read latex or other verbatim rst from xml file
        and post-process it for swig/sphinx-apidoc.

        Parameters
        ----------
        node : xml node

        The current formula is appended to latex_forms.

        Note
        ----
        * latex formula in doxygen comments leads to a real mess 
        during xml->doxy2swig-->swig ... process. 
        Use \rst / \endrst in doxygen files to make it work!

        """    
        self.start_new_paragraph()

        if not node.firstChild.data.strip().startswith("embed:rst"):
            self.subnode_parse(node, pieces=[''], indent=4)
        else:
            content = node.firstChild.data
            if content.strip().startswith("embed:rst:leading-asterisk"):

                lines = content.splitlines()
                # Replace the first * on each line with a blank space
                lines = map(lambda text: text.replace("*", " ", 1), lines)
                content = "\n".join(lines)

            # Remove the first line which is "embed:rst[:leading-asterisk]"
            text = '\n'.join(content.split(u'\n')[1:])

            # Remove starting whitespace
            text = textwrap.dedent(text)
            __id = next(self.__ids)
            node_factory = nf.create_node_factory()
            rst_node = node_factory.displaymath()
            rst_node["latex"] = text
            rst_node["label"] = None
            rst_node["nowrap"] = False
            rst_node["docname"] = None
            rst_node["number"] = None

            self.latex_forms[__id] = rst_node
            id_formula = 'FORMULA' + str(__id) + '_'
            #self.add_text(' ')
            self.add_text(id_formula)

    def get_memberdef_nodes_and_signatures(self, node, kind):
        """Collects the memberdef nodes and corresponding signatures that
        correspond to public function entries that are at most depth 2 deeper
        than the current (compounddef) node. Returns a dictionary with 
        function signatures (what swig expects after the %feature directive)
        as keys, and a list of corresponding memberdef nodes as values."""
        sig_dict = {}
        sig_prefix = ''
        if kind in ('file', 'namespace'):
            ns_node = node.getElementsByTagName('innernamespace')
            if not ns_node and kind == 'namespace':
                ns_node = node.getElementsByTagName('compoundname')
            if ns_node:
                sig_prefix = self.extract_text(ns_node[0]) + '::'
        elif kind in ('class', 'struct'):
            # Get the full function name.
            cn_node = node.getElementsByTagName('compoundname')
            sig_prefix = self.extract_text(cn_node[0]) + '::'

        md_nodes = self.get_specific_subnodes(node, 'memberdef', recursive=2)
        for n in md_nodes:
            
            if n.attributes['prot'].value != 'public':
                continue
            if n.attributes['kind'].value in ['variable']:
                continue
            if not self.get_specific_subnodes(n, 'definition'):
                continue
            name = self.extract_text(self.get_specific_subnodes(n, 'name'))
            if name[:8] == 'operator':
                continue
            sig = sig_prefix + name
            if sig in sig_dict:
                sig_dict[sig].append(n)
            else:
                sig_dict[sig] = [n]
        return sig_dict

            
    def handle_typical_memberdefs(self, signature, memberdef_nodes):
        """Overload doxy2swig method to complete features list
        """
        if len(memberdef_nodes) == 1 or not self.with_overloaded_functions:
            self.handle_typical_memberdefs_no_overload(signature, memberdef_nodes)
            return

        self.add_text(['\n', '%feature("docstring") ', signature, ' "'])
        #for n in memberdef_nodes:
        #    self.add_line_with_subsequent_indent(self.get_function_signature(n))
        self.add_text([signature, '(*args)'])
        self.add_text('\n')
        self.add_text('\n')
        text = 'multiple signatures available, check prototypes below.'
        self.add_text(['.. centered:: **Warning -** :ref:`Overloaded function <overloaded_functions>` : ' + text, '\n','\n'])
        #self.add_text('.. rubric:: Function prototypes')
        self.add_text('\n')
        self.add_text('\n')
        
        for n in memberdef_nodes:
            self.add_text('\n')
            # Function name and prototype
            self.add_text(['.. py:function:: ', self.get_function_signature(n), '\n'] )
            #self.add_text(['* **', self.get_function_signature(n), '**'])
            #            self.add_line_with_subsequent_indent('* **' + self.get_function_signature(n) + '**')

            # Parameters and descriptions
            self.subnode_parse(n, pieces=[], indent=0, ignore=['definition', 'name'])
            self.add_text('\n')
        self.add_text(['";', '\n'])


    def get_function_signature(self, node):
        """Returns the function signature string for memberdef nodes."""
        name = self.extract_text(self.get_specific_subnodes(node, 'name'))
        if self.with_type_info:
            argsstring = self.extract_text(self.get_specific_subnodes(node, 'argsstring'))
        else:
            argsstring = []
            param_id = 1
            for n_param in self.get_specific_subnodes(node, 'param'):
                declname = self.extract_text(self.get_specific_subnodes(n_param, 'declname'))
                if not declname:
                    declname = 'arg' + str(param_id)
                defval = self.extract_text(self.get_specific_subnodes(n_param, 'defval'))
                if defval:
                    defval = '=' + defval
                argsstring.append(declname + defval)
                param_id = param_id + 1
            argsstring = '(' + ', '.join(argsstring) + ')'
        rtype = self.extract_text(self.get_specific_subnodes(node, 'type'))
        argsstring =  self.parse_typemap(argsstring)
        function_definition = name + argsstring
        if rtype != '' : #and type != 'void':
            rtype =  self.parse_typemap(' -> ' + rtype)
            function_definition = function_definition + rtype
        return function_definition

    
    def make_constructor_list(self, constructor_nodes, classname):
        """Produces the "Constructors" section and the constructor signatures
        (since swig does not do so for classes) for class docstrings."""
        if constructor_nodes == []:
            return
        self.add_text(['\n', '.. rubric:: Constructors', '\n', '\n'])
        for n in constructor_nodes:
            self.add_text('\n')
            self.add_text(['.. py:function:: ', self.get_function_signature(n), '\n'] )
            self.add_text('\n')
            self.subnode_parse(n, pieces = [], indent=0, ignore=['definition', 'name'])
            self.add_text('\n')

    def make_attribute_list(self, node):
        """Produces the "Attributes" section in class docstrings for public
        member variables (attributes).
        """
        atr_nodes = []
        for n in self.get_specific_subnodes(node, 'memberdef', recursive=2):
            if n.attributes['kind'].value == 'variable' and n.attributes['prot'].value == 'public':
                atr_nodes.append(n)
        if not atr_nodes:
            return
        self.add_text(['\n', 'Attributes',
                       '\n', '----------'])
        for n in atr_nodes:
            name = self.extract_text(self.get_specific_subnodes(n, 'name'))
            self.add_text(['\n ', name, ' : '])
            argstring = self.parse_typemap(self.extract_text(self.get_specific_subnodes(n, 'type')))
            self.add_text(['`', argstring, '`'])
            self.add_text('  \n')
            restrict = ['briefdescription', 'detaileddescription']
            self.subnode_parse(n, pieces=[''], indent=4, restrict=restrict)


            
    def handle_typical_memberdefs_no_overload(self, signature, memberdef_nodes):
        """Produce standard documentation for memberdef_nodes."""
        for n in memberdef_nodes:
            self.add_text(['\n', '%feature("docstring") ', signature, ' "'])
            self.add_text([self.get_function_signature(n), '\n'])
            self.subnode_parse(n, pieces=[], ignore=['definition', 'name'])
            self.add_text(['";', '\n'])

    def do_simplesect(self, node):
        kind = node.attributes['kind'].value
        # We do not need these informations in python doc
        # (neither anywhere else indeed ...)
        if kind in ('date', 'rcs', 'version', 'author'):
            return
        self.start_new_paragraph()
        if kind == 'warning':
            self.subnode_parse(node, pieces=['**Warning**: ',''], indent=4)
        elif kind == 'see':
            self.subnode_parse(node, pieces=['See also: ',''], indent=4)
        elif kind == 'return':
            if self.indent == 0:
                pieces = ['Returns', '\n', len('Returns') * '-', '\n', '']
            else:
                pieces = ['Returns:', '\n', '']
            self.subnode_parse(node, pieces=pieces)
        else:
            self.subnode_parse(node, pieces=[kind + ': ',''], indent=4)

    def do_compounddef(self, node):
        """Overload doxy2swig method to complete features list

        "This produces %feature("docstring") entries for classes, and handles
        class, namespace and file memberdef entries specially to allow for 
        overloaded functions. For other cases, passes parsing on to standard
        handlers (which may produce unexpected results).

        Remark : the only differences with base class method are the
        call to self.features.append.
        """
        kind = node.attributes['kind'].value
        # -- Prepare doctstring for classes and structs --
        if kind in ('class', 'struct'):
            # Check public/private status
            prot = node.attributes['prot'].value
            if prot != 'public':
                return
            self.add_text('\n\n')

            # Get class name
            classdefn = self.extract_text(self.get_specific_subnodes(node, 'compoundname'))
            classname = classdefn.split('::')[-1]
            self.features.append(classdefn)
            self.add_text('%feature("docstring") {0}"'.format(classdefn))
            if self.with_constructor_list:
                # Parse memberdef nodes and look for constructors
                # match if definition == classdefn::classdefn and status (prot) == public
                constructor_nodes = []
                for n in self.get_specific_subnodes(node, 'memberdef', recursive=2):
                    if n.attributes['prot'].value == 'public':
                        if self.extract_text(self.get_specific_subnodes(n, 'definition')) == classdefn + '::' + classname:
                            constructor_nodes.append(n)
                # For each constructor, append signature to self.pieces.
                #for n in constructor_nodes:
                #    self.add_line_with_subsequent_indent(self.get_function_signature(n))
                #    self.add_text('\n')

            # Get description of the class
            names = ('briefdescription','detaileddescription')
            sub_dict = self.get_specific_nodes(node, names)
            for n in ('briefdescription','detaileddescription'):
                if n in sub_dict:
                    self.parse(sub_dict[n])

            # Name of the original (C/C++) header.
            sub_list = self.get_specific_subnodes(node, 'includes')
            if sub_list:
                self.parse(sub_list[0])


            # List and prototypes for constructors
            if self.with_constructor_list:
                # Create a list of constructors
                self.make_constructor_list(constructor_nodes, classname)

            if self.with_attribute_list:
                self.make_attribute_list(node)

            self.add_text(['";', '\n'])

            names = ['compoundname', 'briefdescription','detaileddescription', 'includes']
            self.subnode_parse(node, ignore = names)

        elif kind in ('file', 'namespace'):
            nodes = node.getElementsByTagName('sectiondef')
            for n in nodes:
                self.parse(n)

        # now explicitely handle possibly overloaded member functions.
        if kind in ['class', 'struct','file', 'namespace', 'typedef']:
            # Get a dictionnary of objects (memberdef nodes) in the node
            # - exclude private and protected
            # - exclude variables
            # - dict.keys = name of the member, dict.value = the memberdef node
            md_nodes = self.get_memberdef_nodes_and_signatures(node, kind)
            for sig in md_nodes:
                self.features.append(sig)
                self.handle_typical_memberdefs(sig, md_nodes[sig])
                # typedef
                #for n in md_nodes[sig]:
                #    if n.attributes['kind'].value == 'typedef':
                #        self.parse_typedefs(n)

        # Process enums
        self.parse_enum(node)
                
    def do_parameteritem(self, node):
        # Overload doxy2swig to remove '*'
        self.subnode_parse(node, pieces=['\n',''])
        
    def do_parametername(self, node):
        # if self.pieces != [] and self.pieces != ['* ', '']:
        #     self.add_text(', ')
        data = self.extract_text(node)
        self.add_text([data])

    def run(self):
        """Parses the file set in the initialization.  The resulting
        data is stored in `self.pieces`.
        """
        self.generate()
        
        # Save latex forms into python dict.
        if len(self.latex_forms) > 0:
            with open(self.latex_filename, 'wb') as current:
                pickle.dump(self.latex_forms, current)

        # write swig file
        self.write(self.swig_outputname)
