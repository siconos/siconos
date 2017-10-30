"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.
"""
import sys
import time
import os
from model_tools import create_model, save_model_to_hdf5

visu = False

if __name__ == "__main__":
    # sample freq and time-discretisation

    # if freq is set as input arg ...
    if len(sys.argv) > 1:
        fs = float(sys.argv[1])
        output_freq = int(sys.argv[2])
    else:
        fs = 51200.
        output_freq = 1

    assert fs > 14000.
    # lower frequencies require quadruple prec for exp computation.

    final_time = 1.00
    number_of_modes = 864
    filt_frets = True
    matlab_input_file = './donnees_siconos/pb2_h.mat'
    guitar_model, guitar_string, frets = create_model(
        n_modes=number_of_modes, max_coords=(7.8e-3, .64),
        fe=fs, final_time=final_time,
        output_freq=output_freq,
        frets_file=matlab_input_file,
        filt_frets=filt_frets,
        enable_frets_output=True, visu=False)

    simu = guitar_model.simulation()

    # -- Model setup with dynamics, interaction and simulation --

    k = 1
    print("Start simulation ...")
    start_time = time.clock()
    pos = 1
    nb_frets = len(frets)
    while simu.hasNextEvent():
        if k % 1000 == 0:
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

    # --- Output dir for results ---
    result_dir = 'results'
    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
    filename = 'g_' + str(number_of_modes) + '_' + str(int(fs)) + '.h5'
    filename = os.path.join(result_dir, filename)
    start = time.clock()
    save_model_to_hdf5(guitar_model, guitar_string,
                       matlab_data=matlab_input_file,
                       filename=filename, filt_frets=filt_frets)
    print('write (hdf5)', time.clock() - start)
