"""Functions and tools to create, save and/or load a model to/from hdf files
"""

from guitar import StringDS, Fret, Guitar
# For matlab to python input file conversion
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os


def create_model(n_modes, max_coords=(7.8e-3, .64),
                 fe=15680, initial_time=0., final_time=0.1,
                 output_freq=1, from_matlab=None,
                 frets_file='./donnees_siconos/pb2_h.mat',
                 filt_frets=True, enable_frets_output=False,
                 visu=False, restitution_coeff=0.9):
    """Build string and model

    Parameters
    ----------

    n_modes : int
        spatial discretisation
    max_coords : tuple
        coordinates of the shifted point at t0
    fe : double
        sampling freq
    final_time : double
    from_matlab : string, optional
        radix of input files (matlab) to set frets positions, freq and damping.
        If set, frets_file is ignored.
    frets_file : string
        input file (matlab) to set frets positions on the neck
    filt_frets : bool
        true set interactions only on frets.
        If false, set contact at each dof on the neck.
    enable_frets_output : string
        if 'light' saves interactions data at each time step (only y[0])
        if 'all' saves all interactions data at each time step
        else saves nothing
    visu: bool
        if true plot frets positions on the neck.
    """
    # ======== Description of the string(s) ==========
    # -- Geometry and material --
    G_string = {
        'length': 0.863,
        # diameter = equivalent diameter (A5)
        'diameter': 1.14e-3,
        'density': 6.69e-3,
        'B': 3.5e-5,
        'tension': 191.6,
    }

    # A dictionnary with parameters required to compute quality factor
    damping_parameters = {
        'eta_air': 1.8e-5,
        'rho_air': 1.2,
        'delta_ve': 0.01,  # if fretless : 0.014
        '1/qte': 6e-6}

    # -- Spatial discretisation (modal proj) and initial conditions --
    # position (index) of the max value of initial state
    # --> middle of the string
    #imiddle = int((n_modes + 2) / 2)

    # -- The dynamical system(s) --
    # Warning: 'real dofs' numbers start from 0 to number_of_modes + 1
    # but DS size is number_of_modes, boundary points are ignored.
    string = StringDS(n_modes, geometry_and_material=G_string,
                      damping_parameters=damping_parameters,
                      max_coords=max_coords,
                      from_matlab=from_matlab)

    # -- The interaction(s) between strings and frets --
    # One Interaction is needed for each contact point.

    if from_matlab is not None:
        frets_file = from_matlab + '_h.mat'
        msg = 'Read data from files :\n'
        msg += '- neck profile:' + frets_file
        msg += '\n- eigenfrequencies: ' + from_matlab + '_frequs.mat\n'
        msg += '- damping: ' + from_matlab + '_amortissements.mat\n'
    else:
        msg = 'Read data from file ' + frets_file + ' for neck profile.'
    print(msg)
    all_frets_positions = scipy.io.loadmat(frets_file)['h'][:, 0]
    x = np.linspace(0, string.length, n_modes + 2)
    x = x[1:-1]

    if filt_frets:
        delt = all_frets_positions[1:]
        diff = np.abs(delt - all_frets_positions[:-1])
        ii = np.where(diff > 1e-4)[0][1::2]
        xp = x[ii]
        frets_positions = all_frets_positions[ii]
        nb_frets = frets_positions.size
        frets_indices = ii
    else:
        frets_positions = all_frets_positions
        nb_frets = frets_positions.size
        frets_indices = np.arange(nb_frets)
        xp = x
        
    #dx = string.space_step
    
    #array(np.round(frets_positions / dx), np.int32)
    #frets_y = np.linspace(-0.5e-4, -1e-4, 20)

    frets = []
    interactions = {}
    for ifret in range(nb_frets):
        frets.append(Fret(string,
                          contact_positions=(frets_indices[ifret],
                                             frets_positions[ifret]),
                          restitution_coeff=restitution_coeff))
        interactions[frets[-1]] = string
    nb_frets = len(frets)
    # contact at a point close to left boundary.
    # guitar_fret_left = Fret(guitar_string, contact_positions=(3, -1.e-4),
    #                         restitution_coeff=0.9)

    # -- Model and simulaton --
    # sample freq and time-discretisation
    assert fe >= 14000.
    # lower frequencies require quadruple prec for exp computation.

    print('Ready to start simulation for frequency {0}.'.format(fe))
    print('Save output every {0} time steps.'.format(output_freq))

    # Plot frets position, visual check ...
    if visu:
        plt.plot(x, all_frets_positions, '-')
        plt.plot(xp, frets_positions, 'x')
        plt.title('frets positons on the neck')

    model = Guitar(interactions, [initial_time, final_time],
                   fe, output_freq,
                   enable_interactions_output=enable_frets_output)

    # -- Save inital state --
    # Note about savings:
    #  For each iteration k, we need :
    #  - the current time value, saved in guitar_model.time[k]
    #  - data for ds (positions, velocities ...):
    #     use save_ds_state_modal(k, ds) for each required ds
    #  - data for interactions (impulsion at impact, distance ...)
    #     use save_interaction_state(k, interaction)

    # Warning : during simulation (save_ds_state_modal),
    # values saved are 'modal' values,
    # i.e. to recover real positions, one needs to multiply data with smat.
    # This is done during post-processing to reduce simulation time.
    model.time[0] = initial_time
    model.save_ds_state_modal(0, string)
    if enable_frets_output is not None:#'all' or 'light':
        for j in range(nb_frets):
            model.save_interaction_state(0, frets[j])
    return model, string, frets


def save_model_to_hdf5(model, ds, filename, matlab_data, filt_frets):
    """Save ds states and time instants in hdf5 file

    Parameters
    ----------
    model : Guitar
        the complete model (nsds + simu)
    ds : StringDS
        the string to be saved
    filename : string
        output file name
    matlab_data : string
        path to matlab files used to build the model
    filt_frets: bool
        use all points (if false) or only 'real' frets (if true) to
        create interactions
    """
    # Before saving, ensure all saved displacements are in the same 'state'
    assert ((model._convert == True).all() or (model._convert == False).all())
    mode = 'w'
    filedir = os.path.dirname(filename)
    if not os.path.exists(filedir):
        os.makedirs(filedir)
    h5file = h5py.File(filename, mode)
    # First, save all attributes required to create the current model
    h5file.attrs['number_of_modes'] = ds.n_modes
    h5file.attrs['max_coords'] = ds.max_coords
    h5file.attrs['frequency'] = model.fs
    h5file.attrs['initial_time'] = model.t0()
    h5file.attrs['final_time'] = model.finalT()
    h5file.attrs['output_freq'] = model.output_freq
    h5file.attrs['matlab_data'] = matlab_data
    h5file.attrs['filter frets'] = filt_frets
    h5file.attrs['frets output'] = model.save_interactions
    interactions = model.interactions_linked_to_ds(ds)
    nb_inter = len(interactions)
    h5file.attrs['number of frets'] = nb_inter
    if (model._convert == False).all():
        h5file.attrs['converted'] = True
    time_shape = (model.time.size, )
    time_steps = h5file.create_dataset('times', time_shape,
                                       dtype=np.float64)
    dof_shape = model.data_ds[ds].shape
    ds_states = h5file.create_dataset('dof', dof_shape,
                                      dtype=np.float64)
    time_steps[...] = model.time
    ds_states[...] = model.data_ds[ds]
    if model.save_interactions:
        inter_states = [None, ] * nb_inter
        for ic in range(nb_inter):
            interaction = interactions[ic]
            ishape = model.data_interactions[interaction].shape
            inter_states[ic] = h5file.create_dataset('contact_' + str(ic),
                                                     ishape,
                                                     dtype=np.float64)
            inter_states[ic][...] = model.data_interactions[interaction]
    h5file.close()


def load_model(filename, from_matlab=None):
    """Read hdf file to:
    * load parameters and create model
    * load simulation results

    Parameters
    ----------
    filename : string
         hdf5 file (relative or full path)
    from_matlab : string, optional
         matlab config path + radix (neck shape, frets positions ...)
         if the one in h5file is wrong (e.g. simulation done on another cluster)

    Remark: the hdf file should have been created at the
    end of a simulation using "save_hdf5" method.
    """

    # Load hdf attributes, to create the model
    print('Load model from file ' + filename)
    mode = 'r'
    h5file = h5py.File(filename, mode)
    n_modes = h5file.attrs['number_of_modes']
    max_coords = h5file.attrs['max_coords']
    fe = h5file.attrs['frequency']
    initial_time = h5file.attrs['initial_time']
    final_time = h5file.attrs['final_time']
    output_freq = h5file.attrs['output_freq']
    if from_matlab is None:
        from_matlab = h5file.attrs['matlab_data']
    filt_frets = h5file.attrs['filter frets']
    enable_frets_output = h5file.attrs['frets output']
    # Create a new model object from file content
    guitar_model, guitar_string, frets = create_model(
        n_modes=n_modes, max_coords=max_coords,
        fe=fe, initial_time=initial_time,
        final_time=final_time,
        output_freq=output_freq,
        from_matlab=from_matlab,
        filt_frets=filt_frets,
        enable_frets_output=enable_frets_output,
        visu=False)
    # load hdf data (simulation results) to fill model variables
    # for post-processing
    guitar_model.time[...] = h5file['times']
    guitar_model.data_ds[guitar_string][...] = h5file['dof']
    if guitar_model.save_interactions is not None:
        interactions = guitar_model.interactions_linked_to_ds(guitar_string)
        nb_inter = len(interactions)
        for ic in range(nb_inter):
            inter = interactions[ic]
            ref = 'contact_' + str(ic)
            guitar_model.data_interactions[inter][...] = h5file[ref]
    h5file.close()
    return guitar_model, guitar_string, frets


def load_convert_and_save(filename, from_matlab):
    # Load hdf5 file to set model and string
    ref_model, ref_string, ref_frets = load_model(filename, from_matlab)
    # Convert (modal to real) displacements
    ref_model.convert_modal_output(ref_string)
    mode = 'w'
    outputfilename = 'converted_' + os.path.basename(filename)
    dirname = os.path.dirname(filename)
    outputfilename = os.path.join(dirname, outputfilename)
    h5source = h5py.File(filename, 'r')
    filt_frets = h5source.attrs['filter frets']
    h5source.close()
    save_model_to_hdf5(ref_model, ref_string, outputfilename, from_matlab, filt_frets)

    
