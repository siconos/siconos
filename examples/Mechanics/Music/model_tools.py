"""Functions and tools to create, save and/or load a model to/from hdf files
"""

from guitar import StringDS, Fret, Guitar
# For matlab to python input file conversion
import scipy.io
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os


def build_frets_from_file(string, restitution_coeff, matlab_frets_file, filt, visu=False):

    all_frets_positions = scipy.io.loadmat(matlab_frets_file)['h'][:, 0]
    # Compute pos[i] - pos[i-1] to find 'real' frets
    delt = all_frets_positions[1:]
    diff = np.abs(delt - all_frets_positions[:-1])
    # indices of 'real' frets:
    if filt:
        i_frets = np.where(diff > 1e-4)[0][1::2]
    else:
        i_frets = np.arange(all_frets_positions.size)

    x_frets = string.x[i_frets]
    y_frets = all_frets_positions[i_frets]
    nb_frets = y_frets.size
    frets = []
    interactions = {}
    for ifret in range(nb_frets):
        frets.append(Fret(string,
                          contact_positions=(i_frets[ifret],
                                             y_frets[ifret]),
                          restitution_coeff=restitution_coeff))
        interactions[frets[-1]] = string
   
    # Plot frets position, visual check ...
    if visu:
        plt.plot(string.x, all_frets_positions, '-')
        plt.plot(x_frets, y_frets, 'x')
        plt.legend(["neck", "frets"])
        plt.title('frets positons on the neck')

    return interactions

        
    return x_frets, y_frets, i_frets


def save_simu_to_hdf5(model, ds, filename, matlab_data, filt_frets, restit):
    """Save model, ds, ... to h5 file

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
    restit : restitution coeff (assume all interactions have the same e)
    """
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
    h5file.attrs['interactions output'] = model.interactions_output
    h5file.attrs['length'] = ds.length
    interactions = model.interactions_linked_to_ds(ds)
    nb_inter = len(interactions)
    h5file.attrs['number of frets'] = nb_inter
    h5file.attrs['modal'] = model.modal_values
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
        ishape = (model.nb_time_steps_output + 1, model.interactions_output) # assume inter size == 1
        for ic in range(nb_inter):
            interaction = interactions[ic]
            nbc = interaction.dimension()
            inter_states[ic] = h5file.create_dataset('contact_' + str(ic),
                                                     ishape,
                                                     dtype=np.float64)
            for i in range(model.interactions_output):
                inter_states[ic][:, i] = model.data_interactions[interaction][i]
    h5file.attrs['restit'] = restit
    h5file.close()


def load_model(filename, visu=True):
    """Read hdf file to:
    * load parameters and create model
    * load simulation results

    Parameters
    ----------
    filename : string
         hdf5 file (relative or full path)

    Remark: the hdf file should have been created at the
    end of a simulation using "save_simu_to_hdf5" method.
    """

    # Load hdf attributes, to create the model
    print('Load model from file ' + filename)
    mode = 'r'
    h5file = h5py.File(filename, mode)
    n_modes = h5file.attrs['number_of_modes']
    max_coords = h5file.attrs['max_coords']
    fs = h5file.attrs['frequency']
    initial_time = h5file.attrs['initial_time']
    final_time = h5file.attrs['final_time']
    output_freq = h5file.attrs['output_freq']
    matlab_input = h5file.attrs['matlab_data']
    filt_frets = h5file.attrs['filter frets']
    interactions_output = h5file.attrs['interactions output']
    restit = h5file.attrs['restit']
    length = h5file.attrs['length']
    G_string = {
        'length': length,
        }
    # -- The dynamical system --
    string = StringDS(n_modes, geometry_and_material=G_string,
                      max_coords=max_coords,
                      matlab_input=matlab_input)

    # -- The interactions --
    frets_file = matlab_input + '_h.mat'
    #filt_frets=True
    interactions = build_frets_from_file(string, restit, frets_file, filt_frets, visu)

    frets = list(interactions.keys())
    nb_frets = len(frets)
    
    # -- The nsds --
    guitar_model = Guitar(interactions, [initial_time, final_time],
                          fs, output_freq,
                          interactions_output=interactions_output)
    guitar_model.time[...] = h5file['times']
    guitar_model.data_ds[string][...] = h5file['dof']
    guitar_model.modal_values = h5file.attrs['modal']
    if guitar_model.save_interactions:
        for ic in range(nb_frets):
            interaction = frets[ic]
            ref = 'contact_' + str(ic)
            for i in range(guitar_model.interactions_output):
                guitar_model.data_interactions[interaction][i][...] = h5file[ref][:, i]
    h5file.close()
    return guitar_model, string, frets, restit


def load_convert_and_save(filename):
    # Load hdf5 file to set model and string
    ref_model, ref_string, ref_frets, restit = load_model(filename)
    # Convert (modal to real) displacements
    ref_model.convert_modal_output(ref_string)
    outputfilename = 'converted_' + os.path.basename(filename)
    dirname = os.path.dirname(filename)
    outputfilename = os.path.join(dirname, outputfilename)
    source = h5py.File(filename,'r')
    #restit = source.attrs['restit']
    filt_frets = source.attrs['filter frets']
    matlab_data = source.attrs['matlab_data']
    print("Write new 'converted' file " + outputfilename)
    save_simu_to_hdf5(ref_model, ref_string, matlab_data=matlab_data,
                      filename=outputfilename,
                      filt_frets=filt_frets, restit=restit)
    source.close()
    
def load_data(filename):
    h5source = h5py.File(filename, 'r')

    data_ds = np.asarray(h5source['dof'])
    data_inter = []
    for ic in h5source.items():
        if ic[0].find('contact') > 0:
            data_inter.append(np.asarray(h5source[ic]))
    time = np.asarray(h5source['times'])
    h5source.close()
    return data_ds, data_inter, time
