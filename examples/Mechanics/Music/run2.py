"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.


Bass guitar with frets.
Restitution coeff = 0.

"""
import sys
import time
import os
from model_tools import create_model, save_model_to_hdf5

visu = False

if __name__ == "__main__":
    # sample freq and time-discretisation

    current_path = os.path.dirname(os.path.realpath(__file__))
    # if freq is set as input arg ...
    if len(sys.argv) > 1:
        fs = float(sys.argv[1])
        output_freq = int(sys.argv[2])
    else:
        fs = 51200.
        output_freq = 1

    assert fs > 14000.
    # lower frequencies require quadruple prec for exp computation.

    final_time = 1.
    number_of_modes = 862
    filt_frets = True
    # Data (from_matlab parameter), choose between:
    # - bass_guitar/pb2 : bass with frets
    # - fretless_bass_guitar/bsf
    matlab_input = os.path.join(current_path, 'bass_guitar/pb2')
    #matlab_input = os.path.join(current_path, 'fretless_bass_guitar/bsf')
    guitar_model, guitar_string, frets = create_model(
        n_modes=number_of_modes, max_coords=(7.8e-3, .64),
        fe=fs, final_time=final_time,
        output_freq=output_freq,
        from_matlab=matlab_input,
        filt_frets=filt_frets,
        enable_frets_output='all', visu=False,
        restitution_coeff=0.)

    simu = guitar_model.simulation()

    # -- Model setup with dynamics, interaction and simulation --

    k = 1
    print("Start simulation ...")
    start_time = time.clock()
    pos = 1
    nb_frets = len(frets)
    while simu.hasNextEvent():
        if k % 100000 == 0:
            print('step = ', k, '---- time = ',
                  simu.nextTime(),
                  '------------------------')
        simu.computeOneStep()
        if k % guitar_model.output_freq == 0:
            guitar_model.time[pos] = simu.nextTime()
            start = time.clock()
            guitar_model.save_ds_state_modal(pos, guitar_string)

            if guitar_model.save_interactions:
                for i in range(nb_frets):
                    guitar_model.save_interaction_state(pos, frets[i])
            pos += 1
        k += 1
        simu.nextStep()
    print('End of simulation process. Duration: ', time.clock() - start_time)
    print("nb steps", k)
    # --- Output dir for results ---
    result_dir = os.getcwd()# + '/temp'
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    filename = 'g_' + str(number_of_modes) + '_' + str(int(fs)) + '.h5'
    filename = os.path.join(result_dir, filename)
    start = time.clock()
    save_model_to_hdf5(guitar_model, guitar_string,
                       matlab_data=matlab_input,
                       filename=filename, filt_frets=filt_frets)
    print('write (hdf5)', time.clock() - start)
