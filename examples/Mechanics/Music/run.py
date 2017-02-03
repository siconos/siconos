"""Implementation of vibrating string model, described
in JSV paper (Issanchou 2017) and using Siconos for contact
simulation.
"""
from string_ds import StringDS
from fret import Guitar, Fret
import matplotlib.pyplot as plt
import time

# ---- Description of the string ---
# -- Geometry and material --
G_string = {
    'length': 1.002,
    'diameter': 0.43e-3,
    'density': 1.17e-3,
    'B': 0.0, # 1.78e-5,
    'tension': 180.5,
}

# A dictionnary with parameters required to compute quality factor
# damping_parameters = {
#     'nu_air': 1.8e-5,
#     'rho_air': 1.2,
#     'delta_ve': 4.5e-3,
#     '1/qte': 2.03e-4}

damping_parameters = {
    'nu_air': 0.,
    'rho_air': 0.,
    'delta_ve': 0.,
    '1/qte': 0.}

# -- Spatial discretisation (modal proj) and initial conditions --
number_of_modes = 201
ndof = number_of_modes + 2
imax = int(ndof / 2)
# -- The dynamical system(s) --
# ie guitar strings
guitar_string = StringDS(ndof, geometry_and_material=G_string,
                         damping_parameters=damping_parameters,
                         umax=1.8e-3, imax=imax,
                         modal_form=False)
guitar_string_m = StringDS(ndof, geometry_and_material=G_string,
                           damping_parameters=damping_parameters,
                           umax=1.8e-3, imax=imax, use_sparse=False,
                           modal_form=True)
guitar_string_sparse = StringDS(ndof, geometry_and_material=G_string,
                                damping_parameters=damping_parameters,
                                umax=1.8e-3, imax=imax, use_sparse=True,
                                modal_form=True)

# -- The obstacles (Interactions) --
# ie guitar fret(s)
fret = Fret(guitar_string, position=[imax, -0.00])
fret_m = Fret(guitar_string_m, position=[imax, -0.00])
fret_sp = Fret(guitar_string_sparse, position=[imax, -0.00])

# -- The model to gather frets and strings and simulate the dynamics --
t0 = 0.
tend = 0.11
guitar_model = Guitar(guitar_string, fret, [t0, tend], fs=2.01e5)
guitar_model_m = Guitar(guitar_string_m, fret_m, [t0, tend], fs=2.01e5)
guitar_model_sparse = Guitar(guitar_string_sparse, fret_sp,
                             [t0, tend], fs=2.01e5)

# -- Run the simulation --


def run_simu(model):
    model.save_state(0)
    k = 1
    print("Start simulation ...")
    while model.simu.hasNextEvent():
        if (k%100==0):
            print ('step = ', k, '---- time = ',  model.simu.nextTime(), '------------------------')
        model.simu.computeOneStep()
        model.save_state(k)
        k += 1
        model.simu.nextStep()
    print 'End of simulation process.'


start_time = time.clock()
run_simu(guitar_model)
print 'duration (std form): ', time.clock() - start_time

# start_time = time.clock()
# run_simu(guitar_model_m)
# print 'duration (modal form, no sparse): ', time.clock() - start_time

# start_time = time.clock()
# run_simu(guitar_model_sparse)
# print 'duration (modal form, sparse): ', time.clock() - start_time


plt.ioff()
# -- plot state vs time --
fig1 = guitar_model.plot_state(1, pdffile='model_std.pdf')
#fig2 = guitar_model_m.plot_state(2, pdffile='model_mod.pdf')
#fig3 = guitar_model_sparse.plot_state(3, pdffile='model_mod2.pdf')

fig1.show()

# -- create string animation --
#guitar_model.plot_modes('string.mp4')
