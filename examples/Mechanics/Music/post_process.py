"""Tools and functions to post-process bass guitar simulations results
"""

from model_tools import load_model
import matplotlib.pyplot as plt
import numpy as np


def compute_errors(filelist, indices=None):
    """Compute relative error defined in 3.36 from Clara's manuscript.
    """

    # Load reference and current simulations
    reference_file = filelist[-1]
    ref_model, ref_string, ref_frets = load_model(reference_file)
    ref_model.convert_modal_output(ref_string, indices)
    sref = ref_model.data_ds[ref_string][:, indices]
    sum_ref = sref.sum(0) ** 2

    nbfiles = len(filelist) - 1
    nbpoints = len(indices)    
    error = np.zeros((nbfiles, nbpoints), dtype=np.float64)
    for i in range(nbfiles):
        current_model, current_string, current_frets = load_model(filelist[i])
        current_model.convert_modal_output(current_string, indices)
        scurrent = current_model.data_ds[current_string][:, indices]
        for j in range(nbpoints):
            buff = (sref[:, j] - scurrent[:, j] ) ** 2
            error[i, j] = buff.sum(0) / sum_ref[j]
            error[i, j] = error[i, j] ** 0.5
    return error



if __name__ == "__main__":

    # Get results files names from a campaign
    from simulation_campaigns import campaign_14112017 as campaign09, campaign_16112017 as campaign0
    freqs = campaign09['freqs']
    filelist = campaign09['files']
    indices = [i for i in range(3, 820)]
    errors09 = compute_errors(filelist, indices)

    freqs0 = campaign09['freqs']
    assert (freqs0 == freqs).all()
    filelist = campaign0['files']
    errors0 = compute_errors(filelist, indices)

    plt.figure()
    #plt.plot(np.log10(freqs[:-1]), np.log10(errors), 'o-')
    for j in range(len(indices)):
        plt.loglog(freqs[:-1], errors09[:, j], 'o-')
        plt.loglog(freqs[:-1], errors0[:, j], 'x-')
    plt.xlabel("$F_e$(Hz)")
    plt.ylabel("$L^2$ error")
    plt.savefig("convergence_study.pdf")
