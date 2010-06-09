# Python GetFEM++ interface
#
# Copyright (C) 2004-2009 Yves Renard, Julien Pommier.
#                                                       
# This file is a part of GETFEM++                                         
#                                                                         
# GetFEM++  is  free software;  you  can  redistribute  it  and/or modify it
# under  the  terms  of the  GNU  Lesser General Public License as published
# by  the  Free Software Foundation;  either version 2.1 of the License,  or
# (at your option) any later version.
# This program  is  distributed  in  the  hope  that it will be useful,  but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or  FITNESS  FOR  A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# You  should  have received a copy of the GNU Lesser General Public License
# along  with  this program;  if not, write to the Free Software Foundation,
# Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
#

# this is a rough module for graphical vizualisation using
# the getfem python interface
# It may change in the future
# examples of use can be found in the tests/python directory
# see for example demo_plasticity, demo_stokes_3D_tank_draw.py

# It requires installation of the TVTK module from enthought
# https://svn.enthought.com/enthought/wiki/TVTK

try:
    from enthought.tvtk.api import tvtk
except:
    print "\n\n** Could not load tvtk. Did you install it ?\n"
    print "   ( https://svn.enthought.com/enthought/wiki/TVTK ) **\n\n"
    raise

import os
import sys
import getfem
import numpy
import scipy

def gf_colormap(name):
    if name == 'tripod':
        s=64; s1=20; s2=25; s3=48; s4=55;
        c = []
        for i in range(1,s):
            c1 = max(min((i-s1)/(s2-s1),1),0);
            c2 = max(min((i-s3)/(s4-s3),1),0);
            c += [(1-c2)*((1-c1)*0.7 + c1) + c2,
                  (1-c2)*((1-c1)*0.7) + c2*.8,
                  (1-c2)*((1-c1)*0.7) + c2*.2]
    elif name == 'chouette':
        c = [.8,  1, .8,
             .7, .9, .4,
             .3, .8, .2,
             .1, .7, .4,
             .2, 0.7, 1.0000,
             .3, 0.3, 1.0000,
             1.0, .8, .1,
             1.0, .6, .1,
             1.0, .45, .1,
             1.0, 0.3, .1]
    elif name == 'froid':
        c = [.8, 1, .8,
             .7, .9, .4,
             .3, .8, .2,
             .1, .7, .4,
             .2, 0.7, 1.0000,
             .3, 0.3, 1.0000]
    elif name == 'tank':
        c = [0, 0, 1,
             0, .5, 1,
             0, 1, .5,
             0, 1, 0,
             .5, 1, 0,
             1, .5, 0,
             1, .4, 0,
             1, 0, 0,
             1, .2, 0,
             1, .4, 0,
             1, .6, 0,
             1, .8, 0];
    elif name == 'earth':
        c = [252, 233, 79, #   Butter 1
             247, 222, 30,
             237, 212,  0, #   Butter 2
             216, 180,  0,
             196, 160,  0, #   Butter 3
             138, 226, 52, #   Chameleon 1
             115, 210, 22, #   Chameleon 2
             78, 154,   6]
        c = numpy.array(c) / 255.0;
    c = numpy.array(c);
    c.shape = (-1,3)
    return c

def _getfem_to_tvtk_points(points):
    (N,nbpt) = points.shape
    if N<3:
        points=numpy.concatenate((points,
                                  numpy.zeros([3-N, nbpt])), axis=0)
    points=numpy.array(points.transpose(), 'd')
    return points

class FigureItem:
    def __init__(self, fig):
        self.fig = fig
        self.sl = None
        self.nrefine = 3
        self.actors = None
        self.show_edges = False
        self.show_faces = True
        self.use_scalar_bar = False
        self.mapper = None
        self.lookup_table = None
        self.scalar_data = None
        self.scalar_data_name = None
        self.scalar_data_range = (0,1)
        self.scalar_bar = None
        self.vector_data = None
        self.edges_color = None
        self.edges_width = None
        self.glyph_name = None
        self.glyph_nb_pts = 1000
        self.glyph_scale_factor = 1
        self.tube_color = None
        self.set_colormap('tripod')
    def set_nrefine(self, nr):
        self.nrefine = nr
    def set_scalar_bar(self,v):
        self.use_scalar_bar = v
    def scalar_range(self, *args):
        if (len(args)==0):
            return self.scalar_data_range
        if (len(args)==1):
            self.scalar_data_range = (args[0][0], args[0][1])
        else:
            self.scalar_data_range = (args[0], args[1])
        if self.mapper is not None:
            self.mapper.scalar_range = self.scalar_data_range;
    
    def build_from_mesh(self,m, **args):        
        dim = m.dim();
        if (dim == 2):
            self.sl=getfem.Slice(('none',),m,self.nrefine)
        elif (dim == 3):
            self.sl=getfem.Slice(('boundary',),m,self.nrefine);
        else:
            raise Exception('%d-D Meshes are not supported'%(dim,))
        self.build_from_slice(self.sl, **args)
        
    def build_from_slice(self, sl, **args):
        self.sl = sl
        self.show_faces = args.get('faces', True)
        self.show_edges = args.get('edges', True)
        self.edges_color = args.get('edges_color', (0.1, 0.1, 0.1))
        self.edges_width = args.get('edges_width', 0.7)
        self.glyph_name = args.get('glyph', None)
        self.glyph_scale_factor = args.get('glyph_scale', 1.0)
        self.glyph_nb_pts = args.get('glyph_nb_pts', 1000)
        self.tube_color = args.get('tube_color',(1,1,1))
        self.actor = None

    def dfield_on_slice(self, data):
        mf = None
        if (isinstance(data, tuple)):
            if (len(data) == 2):
                mf = data[0]
                U = data[1]
            elif (len(data) == 1):
                U = data[0]
            else:
                raise Exception, "wrong data tuple.."
        else:
            U = data
        if mf is not None:
            return getfem.compute(mf, U, 'interpolate on', self.sl);
        else:
            return U

    def set_scalar_data(self, data, name='scalars'):
        self.scalar_data = self.dfield_on_slice(data)

        self.scalar_data_name = name
        #self.scalar_data_range = (min(self.scalar_data),
        #                          max(self.scalar_data))
        m = self.scalar_data.mean()
        s = self.scalar_data.std()
        self.scalar_data_range = (max(m - s, self.scalar_data.min()),
                                  min(m + s, self.scalar_data.max()))

    def set_vector_data(self, vdata, name='vectors'):
        d = self.dfield_on_slice(vdata)
        n = self.sl.nbpts()
        if d.size % n != 0:
            raise Exception, "non consistent dimension for data"
        if d.size > n:
            d = d.transpose()
            d.shape = (n,-1)
        self.vector_data = d
        if self.glyph_name is None:
            self.glyph_name = 'default'
        
    def deformation_from_mf(self, mf, U, scale):
        P=self.sl.pts()
        deform = getfem.compute(mf, U, 'interpolate on', self.sl)
        try:
            scale=float(scale)
        except ValueError:
            if scale.endswith('%'):
                a = max(abs(P.max()),abs(P.min()),1e-10);
                b = max(abs(deform.max()),abs(deform.min()));
                scale = float(scale[:-1]) * 0.01 * a/b;
        P=P + scale * deform
        self.sl.set_pts(P)
        #print "deformation!", repr(mf), U.size(), repr(self.sl), "\nDEFORM=",self.deform,"\n"
        sys.stdout.flush()


    def set_colormap(self, c):
        if isinstance(c,str):
            lut = tvtk.LookupTable()
            c=gf_colormap(c)
            lut.number_of_table_values=c.shape[0]
            for i in range(c.shape[0]):
                lut.set_table_value(i,c[i,0],c[i,1],c[i,2],1)
        elif isinstance(c, tvtk.LookupTable):
            lut = c
        else:
            raise Exception, "expected a string or a tvtk.LookupTable"
        self.lookup_table = lut
        if (self.mapper is not None):
            self.mapper.lookup_table = self.lookup_table
        if (self.scalar_bar is not None):
            self.scalar_bar.lookup_table = self.lookup_table
    def vtk_actors(self):
        if (self.actors is None):
            self.actors = []
            points=_getfem_to_tvtk_points(self.sl.pts())
            (triangles,cv2tr)=self.sl.splxs(2);
            triangles=numpy.array(triangles.transpose(), 'I');
            data = tvtk.PolyData(points=points, polys=triangles)
            if self.scalar_data is not None:
                data.point_data.scalars = numpy.array(self.scalar_data)
            if self.vector_data is not None:
                data.point_data.vectors = numpy.array(self.vector_data)
            if self.glyph_name is not None:
                mask = tvtk.MaskPoints()
                mask.maximum_number_of_points = self.glyph_nb_pts
                mask.random_mode = True
                mask.input = data

                if self.glyph_name == 'default':
                    if self.vector_data is not None:
                        self.glyph_name = 'arrow'
                    else:
                        self.glyph_name = 'ball'
                
                glyph = tvtk.Glyph3D()
                glyph.scale_mode = 'scale_by_vector'
                glyph.color_mode = 'color_by_scalar'
                #glyph.scale_mode = 'data_scaling_off'
                glyph.vector_mode = 'use_vector' # or 'use_normal'
                glyph.input = mask.output
                if self.glyph_name == 'arrow':
                    glyph.source = tvtk.ArrowSource().output
                elif self.glyph_name == 'ball':
                    glyph.source = tvtk.SphereSource().output
                elif self.glyph_name == 'cone':
                    glyph.source = tvtk.ConeSource().output
                elif self.glyph_name == 'cylinder':
                    glyph.source = tvtk.CylinderSource().output
                elif self.glyph_name == 'cube':
                    glyph.source = tvtk.CubeSource().output
                else:
                    raise Exception, "Unknown glyph name.."
                #glyph.scaling = 1
                #glyph.scale_factor = self.glyph_scale_factor
                data = glyph.output
                
            if self.show_faces:
##                if self.deform is not None:
##                    data.point_data.vectors = array(numarray.transpose(self.deform))
##                    warper = tvtk.WarpVector(input=data)
##                    data = warper.output
##                lut = tvtk.LookupTable()
##                lut.hue_range = 0.667,0
##                c=gf_colormap('tripod')
##                lut.number_of_table_values=c.shape[0]
##                for i in range(0,c.shape[0]):
##                    lut.set_table_value(i,c[i,0],c[i,1],c[i,2],1)
                

                   
                self.mapper = tvtk.PolyDataMapper(input=data);
                self.mapper.scalar_range = self.scalar_data_range;
                self.mapper.scalar_visibility = True
                # Create mesh actor for display
                self.actors += [tvtk.Actor(mapper=self.mapper)]
            if self.show_edges:
                (Pe, E1, E2)=self.sl.edges();
                if Pe.size:
                    E = numpy.array(numpy.concatenate((E1.transpose(),
                                                       E2.transpose()),
                                                      axis=0), 'I')
                    edges=tvtk.PolyData(points=_getfem_to_tvtk_points(Pe),
                                        polys=E)
                    mapper_edges = tvtk.PolyDataMapper(input=edges);
                    actor_edges = tvtk.Actor(mapper=mapper_edges)
                    actor_edges.property.representation = 'wireframe'
                    #actor_edges.property.configure_traits()
                    actor_edges.property.color = self.edges_color
                    actor_edges.property.line_width = self.edges_width
                    actor_edges.property.ambient = 0.5
                    self.actors += [actor_edges];
            if self.sl.nbsplxs(1):
                # plot tubes
                (seg,cv2seg)=self.sl.splxs(1)
                seg=numpy.array(seg.transpose(),'I')
                data=tvtk.Axes(origin=(0,0,0), scale_factor=0.5, symmetric=1)
                data=tvtk.PolyData(points=points, lines=seg)
                tube = tvtk.TubeFilter(radius=0.4, number_of_sides=10,
                                       vary_radius='vary_radius_off',
                                       input=data)
                mapper = tvtk.PolyDataMapper(input=tube.output)
                actor_tubes = tvtk.Actor(mapper=mapper)
                #actor_tubes.property.representation = 'wireframe'
                actor_tubes.property.color = self.tube_color
                #actor_tubes.property.line_width = 8
                #actor_tubes.property.ambient = 0.5
                    
                self.actors += [actor_tubes]
    
            if self.use_scalar_bar:
                self.scalar_bar = tvtk.ScalarBarActor(title=self.scalar_data_name,
                                                 orientation='horizontal',
                                                 width=0.8, height=0.07)
                self.scalar_bar.position_coordinate.coordinate_system = 'normalized_viewport'
                self.scalar_bar.position_coordinate.value = 0.1, 0.01, 0.0
                self.actors += [self.scalar_bar]
                
            if (self.lookup_table is not None):
                self.set_colormap(self.lookup_table)

        return self.actors

            
    

class Figure:
    def __init__(self, gui='tvtk'):
        self.actors = []
        self.gui = None
        self.renderer = None
        self.items = []
        if gui == 'tvtk':
            self._create_tvtk_window()
        else:
            self._create_ivtk_window()
        self.renderer.background = (1,1,1)
    def _create_tvtk_window(self, size=(500,500)):
        # create a renderer
        self.renderer = tvtk.Renderer()
        # create a render window and hand it the renderer
        self.render_window = tvtk.RenderWindow(size=size)
        self.render_window.add_renderer(self.renderer)
        # create interactor and hand it the render window
        # This handles mouse interaction with window.
        self.interactor = tvtk.RenderWindowInteractor(render_window=self.render_window)
        self.gui = None

    def _create_ivtk_window(self, size=(800,800)):
        from enthought.tvtk.tools import ivtk
        from enthought.pyface.api import GUI

        # Create a GUI instance.
        self.gui = GUI()
        window = ivtk.IVTKWithCrustAndBrowser(size=size)  # Size is optional.
        # Open the window.
        window.open()
        self.renderer = window.scene
        self.render_window = window


    def show_mesh(self, m, **args):
        it = FigureItem(self)
        it.build_from_mesh(m, **args)
        self.actors += it.vtk_actors()
        self.items.append(it)

    def show_mesh_fem(self, mf, **args):
        it = FigureItem(self)
        it.build_from_mesh(mf.linked_mesh(), **args)
        if args.has_key('deformation'):
            it.deformation_from_mf(args.get('deformation_mf',mf),
                                   args['deformation'],
                                   args.get('deformation_scale','10%'));
        if args.has_key('data'):
            it.set_scalar_data(args.get('data'),
                               args.get('scalar_label', 'data'));
        it.set_scalar_bar(args.get('scalar_bar', False))

        if args.has_key('vdata'):
            it.set_vector_data(args.get('vdata'))
            
        self.actors += it.vtk_actors()
        self.items.append(it)

        it.set_colormap(args.get('colormap','earth'));

    def show_slice(self, sl, **args):
        it = FigureItem(self)
        it.build_from_slice(sl, **args)
        
        if args.has_key('data'):
            it.set_scalar_data(args.get('data'),
                               args.get('scalar_label', 'data'));

        it.set_scalar_bar(args.get('scalar_bar', False))

        if args.has_key('vdata'):
            it.set_vector_data(args.get('vdata'))
                               
        self.actors += it.vtk_actors()
        self.items.append(it)

        it.set_colormap(args.get('colormap','chouette'));

    def scalar_range(self, *args):
        if len(self.items):
            if len(args)==0:
                return self.items[-1].scalar_range()
            else:
                for i in self.items:
                    i.scalar_range(*args)
        else:
            raise Exception, "plot something before changing its scalar range!"

##    def scalar_bar(self):
##        if len(self.items):
##            self.items[-1].set_scalar_bar(True)

    def set_colormap(self, c):
        if (len(self.items)):
            self.items[-1].set_colormap(c)
    def show(self, mf, **args):
        if isinstance(mf, getfem.MeshFem):
            self.show_mesh_fem(mf, **args)
        elif isinstance(mf, getfem.Mesh):
            self.show_mesh(mf, **args)
        elif isinstance(mf, getfem.Slice):
            self.show_slice(mf, **args)
        else:
            raise TypeError, "argument must be a drawable getfem object"
    def loop(self):
        for a in self.actors:
            self.renderer.add_actor(a)
        if self.gui:
            self.renderer.reset_zoom()
            self.gui.start_event_loop()
        else:
            self.interactor.start()

    def export_picture(self, filename):
        w2if = tvtk.WindowToImageFilter()
        w2if.magnification = 2
        w2if.input = self.render_window
        ex = tvtk.PNGWriter()
        ex.file_name = filename
        ex.input = w2if.output
        ex.write()


##mf=getfem.MeshFem('load','tripod.mf');
##m=mf.linked_mesh()
##mfvm=getfem.MeshFem('load','tripod.mfe',m);
##U = numarray.fromfile('tripod.U','d')
##VM = numarray.fromfile('tripod.VM','d')
##fig = Figure()
##fig.show(mfvm, data=VM, deformation_mf=mf, deformation=U,
##         scalar_bar=True, scalar_label='Von Mises Stress')
##fig.set_colormap('chouette')
##fig.loop()
##sys.exit(1)

###plot_mesh(m);
###def plot_mesh(m):
##if m:
##    p = tvtk.Property(representation='wireframe') 
##    p.representation = 's' 
##    p.representation 
##    # -> 'surface' 
##    #p.configure_traits()

##    sl=getfem.Slice(('boundary',),m,2);
##    (Pe, E1, E2)=sl.edges();
    
    
##    points=sl.pts(); points.transpose();
##    points=array(points);
##    (triangles,cv2tr)=sl.splxs(2);
##    triangles.transpose();
##    triangles=array(triangles);
##    mesh = tvtk.PolyData(points=points, polys=triangles)

##    Pe.transpose();
##    E1.transpose()
##    edges=tvtk.PolyData(points=array(Pe),polys=array(E1))

##    print mesh.get()


##    #data = array([[0,0,0,10], [1,0,0,20],
##    #              [0,1,0,20], [0,0,1,30]], 'f')
##    #triangles = array([[0,1,3], [0,3,2],
##    #                   [1,2,3], [0,2,1]])
##    #points, temperature = data[:,:3], data[:,-1]
##    #mesh = tvtk.PolyData(points=points, polys=triangles)
##    #mesh.point_data.scalars = temperature


##    ### TVTK PIPELINE
##    if 0:
##        # create a renderer
##        renderer = tvtk.Renderer()
##        # create a render window and hand it the renderer
##        render_window = tvtk.RenderWindow(size=(400,400))
##        render_window.add_renderer(renderer)

##        # create interactor and hand it the render window
##        # This handles mouse interaction with window.
##        interactor = tvtk.RenderWindowInteractor(render_window=render_window)
##    else:
##        # Create a GUI instance.
##        gui = GUI()

##        # Create and open an IVTK application window that has an embedded TVTK
##        # pipeline browser and an interactive Python interpreter shell via
##        # PyCrust.  If you don't want all these you can choose between the
##        # following classes in ivtk -- IVTK, IVTKWithCrust, IVTKWithBrowser
##        # and IVTKWithCrustAndBrowser.
##        window = ivtk.IVTKWithCrustAndBrowser(size=(800,600))  # Size is optional.

##        # Open the window.
##        window.open()
        
##        viewer = window
##        renderer = viewer.scene

 
##    # Set the mapper to scale temperature range
##    # across the entire range of colors
##    mapper = tvtk.PolyDataMapper(input=mesh);
##    mapper_edges = tvtk.PolyDataMapper(input=edges);
##    print mapper
 
##    #mapper = tvtk.PolyDataMapper(input=mesh)
##    #mapper.scalar_range = min(temperature), max(temperature)
    
##    # Create mesh actor for display
##    actor = tvtk.Actor(mapper=mapper)
##    actor_edges = tvtk.Actor(mapper=mapper_edges)
##    actor_edges.property.representation = 'wireframe'
    
##    # Create a scalar bar
##    scalar_bar = tvtk.ScalarBarActor(title="Temperature",
##                                     orientation='horizontal',
##                                     width=0.8, height=0.17,
##                                     lookup_table = mapper.lookup_table)
##    scalar_bar.position_coordinate.coordinate_system = 'normalized_viewport'
##    scalar_bar.position_coordinate.value = 0.1, 0.01, 0.0
    
##    # Use the ScalarBarWidget so we can drag the scalar bar around.
##    #sc_bar_widget = tvtk.ScalarBarWidget(interactor=interactor,
##    #                                     scalar_bar_actor=scalar_bar)
    
##    # Now add the actors to the renderer and start the interaction.
##    renderer.add_actor(actor)
##    renderer.add_actor(actor_edges)
##    #interactor.initialize()
##    # Enable the widget so the scalar bar can be seen.  Press 'i' to
##    # disable the widget.
##    #sc_bar_widget.enabled = True
##    #interactor.start()
##    gui.start_event_loop()


##    print "finished!"


##    #f=mlab.figure()
##    #f.add(mlab.TriMesh(points,triangles));
    


