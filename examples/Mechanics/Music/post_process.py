"""Tools and functions to post-process bass guitar simulations results
"""

from model_tools import load_model
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import scipy.io as sio


def save_dof(fileslist, dof, outputfile):
    """Save displacement for a given dof, in mat file
    Extract results from a list of simulations 
    (e.g. convergence study, to plot dof = f(frequency))
    
    Parameters
    ----------
    fileslist : list of string
        list of files (h5) from which displacements must be loaded
    dof : int
        degree of freedom of interest
    outputfile : string
        name of resulting file (.mat)
    """

    files = list(fileslist.values())
    freqs = list(fileslist.keys())
    nbfiles = len(files)
    current_model, current_string, current_frets, e = load_model(files[0])
    # Number of time instants taken into account
    nbtimes = current_model.data_ds[current_string].shape[1]
    print(nbtimes)
        
    #files = [f for f in files if f != reference_file]
    # Number of freqs taken into account
    nbfiles = len(files)
    # Number of dofs where errors are computed
    nbpoints = 1
    # Number of time instants taken into account
    nbtimes = current_model.data_ds[current_string].shape[1]
    # Number of contact points
    nbfrets = len(current_frets)

    # Results buffers:
    dof_val= np.zeros((nbfiles, nbtimes), dtype=np.float64)

    # Compute errors for all freqs
    for i in range(nbfiles):
        current_model, current_string, current_frets, e = load_model(files[i])
        dof_val[i, :] = current_model.data_ds[current_string][dof, :]

    dofname = 'dof_' + str(dof)
    matdict = {'freqs': freqs, dofname : dof_val}
    sio.savemat(outputfile, matdict)
    return dof_val



def compute_errors(fileslist, dofs=None, shift=1, savedofs=None, ifref=-1):
    """Compute relative error defined in 3.36 from Clara's manuscript.

    Parameters
    ----------
    fileslist : list of string
        list of h5 files, each of them holding results of a simulation for 
        a given frequency.
    dofs : list of int
        degrees of freedom (ds) of interest
    savedofs : list of int
        list of dofs to be saved in matlab file (default : none)
    ifref : id (in filelist) of the simu (i.e. freq) used as reference
    for error computation. Default = max freq.
    """
    #fref = max(list(fileslist.keys()))
    fref = list(fileslist.keys())[ifref]
    print("compute errors with reference freq = " + str(fref))
    # Load reference and current simulations,
    reference_file = fileslist[fref]
    ref_model, ref_string, ref_frets, e = load_model(reference_file)
    # if indices is none, errors are computed for all degrees of freedom.
    if dofs is None:
        dofs = [i for i in range(ref_string.dimension())]

    sref = ref_model.data_ds[ref_string][dofs, ::shift]
    sum_ref = (sref ** 2).sum(1)
    freqs = []
    tref = ref_model.time[::shift]
    files = list(fileslist.values())[:ifref]
    files = [f for f in files if f != reference_file]
    # Number of freqs taken into account
    nbfiles = len(files)
    # Number of dofs where errors are computed
    nbpoints = len(dofs)
    # Number of contact points
    nbfrets = len(ref_frets)
    # Number of time instants taken into account
    nbtimes = ref_model.data_ds[ref_string].shape[1]
    xref = ref_string.x[dofs]
    if savedofs is not None:
        dofs_val = {}
        for dof in dofs:
            dofs_val[dof]= np.zeros((nbfiles + 1, nbtimes), dtype=np.float64)
            dofs_val[dof][-1, :] = ref_model.data_ds[ref_string][dof, :]
        
    # Results buffers:
    # errors[i, j] = error for freq number i at dof j
    errors = np.zeros((nbfiles, nbpoints), dtype=np.float64)
    # ymin[i, j] = minimal value (through time instants) of distance at contact j for freq i
    ymin = np.zeros((nbfiles + 1, nbfrets), dtype=np.float64)
    # Compute ymin for reference model
    for j in range(nbfrets):
        ymin[-1, j] = (ref_model.data_interactions[ref_frets[j]][0][:]).min()

    # Compute errors for all freqs
    for i in range(nbfiles):
        current_model, current_string, current_frets,e = load_model(files[i])
        scurrent = current_model.data_ds[current_string][dofs, :]
        # Ensure time instants are the same for both models (ref an current)
        time_step = current_model.time_step
        tcurrent = current_model.time
        assert np.allclose(tref, tcurrent), 'Error: time instants are different.'
        xcurrent = current_string.x[dofs]
        assert np.allclose(xref, xcurrent)
        #errors[i, :] = np.sqrt(time_step * (((sref - scurrent) ** 2).sum(1)) / sum_ref)
        errors[i, :] = np.sqrt((((sref - scurrent) ** 2).sum(1)) / sum_ref)
        if savedofs is not None:
            for dof in dofs:
                dofs_val[dof][i, :] = current_model.data_ds[current_string][dof, :]

        for j in range(nbfrets):
            ymin[i, j] = (current_model.data_interactions[current_frets[j]][0][:]).min()
        freqs.append(current_model.fs)

    if savedofs is not None:
        frs = list(fileslist.keys())
        matdict = {'freqs': frs}
        outputfile = savedofs
        for dof in dofs:
            dofname = 'dof_' + str(dof)
            matdict[dofname] = dofs_val[dof]
            
        sio.savemat(outputfile, matdict)
    return errors, ymin, freqs, xref

def compute_errors_rel(fileslist, dofs=None, shift=1):
    """Compute relative error defined in 3.36 from Clara's manuscript.
    """
    fref = max(fileslist.keys())

    # Load reference and current simulations,
    reference_file = fileslist[fref]
    ref_model, ref_string, ref_frets,e = load_model(reference_file)
    # if dofs is none, errors are computed for all degrees of freedom.
    if dofs is None:
        dofs = [i for i in range(ref_string.dimension())]

    sref = ref_model.data_ds[ref_string][dofs, ::shift]
    #sum_ref = (sref ** 2).sum(1)
    freqs = []
    tref = ref_model.time[::shift]
    files = list(fileslist.values())
    files = [f for f in files if f != reference_file]
    # Number of freqs taken into account
    nbfiles = len(files)
    # Number of dofs where errors are computed
    nbpoints = len(dofs)
    # Number of contact points
    nbfrets = len(ref_frets)

    # Results buffers:
    # errors[i, j] = error for freq number i at dof j
    errors = np.zeros((nbfiles, nbpoints), dtype=np.float64)
    # ymin[i, j] = minimal value (through time instants) of distance at contact j for freq i
    ymin = np.zeros((nbfiles + 1, nbfrets), dtype=np.float64)

    # Compute ymin for reference model
    for j in range(nbfrets):
        ymin[-1, j] = (ref_model.data_interactions[ref_frets[j]][0][:]).min()

    # Compute errors for all freqs
    for i in range(nbfiles):
        current_model, current_string, current_frets,e = load_model(files[i])
        scurrent = current_model.data_ds[current_string][dofs, :]
        # Ensure time instants are the same for both models (ref an current)
        time_step = current_model.time_step
        tcurrent = current_model.time
        assert np.allclose(tref, tcurrent), 'Error: time instants are different.'
        
        #errors[i, :] = np.sqrt(time_step * (((sref - scurrent) ** 2).sum(1)) / sum_ref)
        errors[i, :] = ((((sref - scurrent) / sref ) ** 2).sum(1)) ** 0.5

        
        for j in range(nbfrets):
            ymin[i, j] = (current_model.data_interactions[current_frets[j]][0][:]).min()
        freqs.append(current_model.fs)
    return errors, ymin, freqs



def check_time_vectors(filelist):
    """Ensure compatibility between time vectors from a list of simulation
    """

    # Load reference and current simulations,
    reference_file = filelist[-1]
    ref_model, ref_string, ref_frets,e = load_model(reference_file)
    tref = ref_model.time
    nbfiles = len(filelist) - 1
    for i in range(nbfiles):
        if os.path.exists(filelist[i]):
            current_model, current_string, current_frets,e = load_model(filelist[i])
            tcurrent = current_model.time[:]
            print("check i ... ", i)
            assert np.allclose(tcurrent, tref)
        else:
            print('Missing file ' + filelist[i] + ' - Skip')


def plot_campaign(campaign, dofs, fig=39):
    freqs = campaign['freqs']
    filelist = campaign['files_converted']
    myticks = ['o-', 'x-', '^-', '*-', '+-'] * 200 
    errors = compute_errors(filelist, dofs)
    print(errors)
    print(freqs)
    plt.figure(fig)
    nbpoints = errors.shape[1]
    nbfreqs = errors.shape[0]
    for j in range(nbpoints):
        plt.loglog(freqs[:nbfreqs], errors[:, j], myticks[j])
    plt.xlabel("$F_e$(Hz)")
    plt.ylabel("$L^2$ error")
    plt.savefig("convergence_study.pdf")
 


def plot_errors(errors, freqs, dof, xref, iplot=0, figures_path=None):
    plt.figure(iplot, figsize=(17, 17))
    tickslabels = ['o:', '^:', '*:'] * 3
    i=0
    leg = []
    for name in errors:
        ref_error = errors[name]
        ref_freqs = freqs[name]
        plt.loglog(ref_freqs, ref_error[:, dof], tickslabels[i])
        leg.append(name)
        plt.grid(True,which="both",ls="-")
        i+=1
    plt.title("Convergence, x=" + str(xref[name][dof])) # ...
    #ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
    #ax.yaxis.set_ticks(np.arange(-2, 2, 0.025))
    #    major_ticks = np.arange(0, 1, 0.2)                                              
    #minor_ticks = np.arange(0, 1, 0.05)                                               
    #ax = plt.gca()
    #ax.set_yticks(major_ticks)                                                       
    #ax.set_yticks(minor_ticks, minor=True)                                           
    #ax.grid(which='both')   
    plt.legend(leg)
    #plt.grid('on')
    plt.xlabel('Fe(Hz)')
    plt.ylabel('error')
    if figures_path is not None:
        plt.savefig(os.path.join(figures_path, 'cvg_x='  + str(xref[name][dof]) + '.pdf'))
    return plt

def plot_y(ymin, freqs, iplot=0):

    plt.figure(iplot, figsize=(17, 8))
    #leg = []
    for name in ymin:
        print(ymin[name].shape)
        for j in range(ymin[name].shape[1]):
            #plt.plot(freqs[name], ymin[name][:,j])
            plt.semilogx(freqs[name][:], ymin[name][:-1,j],'x:')

    plt.grid('on')
    plt.xlabel('Fe(Hz)')
    plt.ylabel('min value of distances at contact')
    return plt
            
            
def plot_durations(timers):
    tick1 = ['o--', '^--', 'x--']
    tick2 = ['o', '^', 'x']
    fig0, ax1 = plt.subplots(figsize=(17, 8))
    
    ax2 = ax1.twinx()
    leg =[]
    i = 0
    for name in timers:
        timers[name].sort(0)
        ax1.semilogx(timers[name][:, 0], timers[name][:, 1]/60., tick1[i])
        ax2.semilogx(timers[name][:, 0], timers[name][:, 2], tick2[i])
        i += 1
        leg.append(name)
    ax1.legend(leg)
    ax1.grid('on')
    ax1.set_xlabel('Fe (Hz)')
    ax1.set_ylabel('duration (min)')
    ax2.set_ylabel('host number')
    return plt
    #fig.tight_layout()
    
if __name__ == "__main__":

    # Get results files names from a campaign
    from simulation_campaigns import campaign_14112017 as campaign09, campaign_16112017 as campaign0
    freqs = campaign09['freqs']
    filelist = campaign09['files']
    dofs = [i for i in range(3, 820)]
    errors09 = compute_errors(filelist, dofs)

    freqs0 = campaign09['freqs']
    assert (freqs0 == freqs).all()
    filelist = campaign0['files']
    errors0 = compute_errors(filelist, dofs)

    plt.figure()
    #plt.plot(np.log10(freqs[:-1]), np.log10(errors), 'o-')
    for j in range(len(dofs)):
        plt.loglog(freqs[:-1], errors09[:, j], 'o-')
        plt.loglog(freqs[:-1], errors0[:, j], 'x-')
    plt.xlabel("$F_e$(Hz)")
    plt.ylabel("$L^2$ error")
    plt.savefig("convergence_study.pdf")
