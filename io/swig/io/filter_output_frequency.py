#!/usr/bin/env python
import sys
import os
import vtk
import numpy

from siconos.io.mechanics_io import Hdf5

## the best way to dump all data
#  $ h5dump toto.hdf5 > toto.txt

import h5py

import getopt

def usage():
    """
    {0} <hdf5 file>
    """.format(sys.argv[0])
    print '{0}: Usage'.format(sys.argv[0])
    print """
    {0} [--help] [--output_frequency=n] [--output_filename=]  <hdf5 file>
    """

try:
    opts, args = getopt.gnu_getopt(sys.argv[1:], '',
                                   ['output_frequency=',
                                    'output_filename=',
                                    'cf-scale='])
except getopt.GetoptError, err:
        sys.stderr.write('{0}\n'.format(str(err)))
        usage()
        exit(2)

min_time = None
max_time = None
scale_factor = 1
output_frequency=10
output_filename=None
for o, a in opts:

    if o == '--help':
        usage()
        exit(0)

    elif o == '--output_frequency':
        output_frequency = float(a)

    elif o == '--output_filename':
        out_filename=a

    elif o == '--cf-scale':
        scale_factor = float(a)


if len(args) > 0:
    io_filename = args[0]
    if output_filename == None:
        out_filename=''.join(io_filename.rsplit('.')[:-1])+'-filtered.hdf5'
else:
    usage()
    exit(1)



with Hdf5(io_filename=io_filename, mode='r') as io:
    with Hdf5(io_filename=out_filename, mode='w') as out:



        hdf1 = io._out
        hdf2 = out._out

        # copy /data/input
        hdf2.__delitem__('data/input')
        h5py.h5o.copy(hdf1.id, "data/input", hdf2.id, "data/input")
        # copy /data/nslaws
        hdf2.__delitem__('data/nslaws')
        h5py.h5o.copy(hdf1.id, "data/nslaws", hdf2.id, "data/nslaws")
        # copy /data/ref
        hdf2.__delitem__('data/ref')
        h5py.h5o.copy(hdf1.id, "data/ref", hdf2.id, "data/ref")

        print('***************************************************** ')
        print('************ Parsing simulation data ****************')
        print('***************************************************** ')

        def load():

            ispos_data = io.static_data()
            idpos_data = io.dynamic_data()

            icf_data = io.contact_forces_data()[:]

            isolv_data = io.solver_data()

            return ispos_data, idpos_data, icf_data, isolv_data

        spos_data, dpos_data, cf_data, solv_data = load()

        #print('io._data',io._data)
        #print('static position data : spos_data',spos_data)
        #print('spos_data.value',spos_data.value)
        #print('dynamic position data : dpos_data',dpos_data)
        print('dpos_data.value',dpos_data.value)
        print('cf_data',cf_data)
        #print('solv_data',solv_data)

        times = list(set(dpos_data[:, 0]))
        times.sort()
        print('len(times)',len(times))

        if (len(times) ==0 ):
            print('no results in the hdf5 file')
        else:


            print('Results for ',len(times),' steps in the hdf5 file')
            #ndyna = len(numpy.where(dpos_data[:, 0] == times[0])) does not work
            ndyna = len(dpos_data[:, 0])/len(times)
            print('ndyna =', ndyna)
            if len(spos_data) > 0:
                nstatic = len(numpy.where(spos_data[:, 0] == times[0]))
                nstatic = spos_data.shape[0]
            else:
                nstatic = 0
            print('nstatic =', nstatic)
            #    instances = set(dpos_data[:, 1])



        # filtering


        p=0
        current_line=0
        for k in range(len(times)):
            #print(times[k])
            if (k+1  < len(times) ):
                time_step=times[k+1]-times[k]
                #print('time_step',time_step)

            if (k%output_frequency==0):
                if k==0 :
                    print('filter for k',k,'at times', times[k], 'p', p)
                    out._dynamic_data.resize((p+1)*ndyna,0)
                    out._dynamic_data[p*ndyna:(p+1)*ndyna,:] = numpy.array(dpos_data[k*ndyna:(k+1)*ndyna,:])
                    #print('times',dpos_data[k*ndyna:(k+1)*ndyna,0])

                    out._static_data.resize((p+1)*nstatic,0)
                    out._static_data[p*nstatic:(p+1)*nstatic,:] = numpy.array(spos_data[0:nstatic,:])

                    out._solv_data.resize((p+1),0)
                    out._solv_data[p:(p+1),:] = numpy.array(solv_data[0:1,:])


                    id_f = numpy.where(abs(cf_data[:, 0] - times[k]) < time_step*1e-5)[0]
                    if len(id_f) == 0:
                        print('no contact data at time',times[k])
                    else:
                        #print('index of contact :', min(id_f), max(id_f))
                        out._cf_data.resize(max(id_f)+1,0)
                        out._cf_data[min(id_f):max(id_f),:] = numpy.array(cf_data[min(id_f):max(id_f),:])
                        current_line = max(id_f)
                else:
                    print('filter for k',k,'at times', times[k], 'p', p)
                    out._dynamic_data.resize((p+1)*ndyna,0)
                    #print( dpos_data[k*ndyna:(k+1)*ndyna,:])
                    #print('times',dpos_data[(k+1)*ndyna:(k+2)*ndyna,0])
                    out._dynamic_data[p*ndyna:(p+1)*ndyna,:] = dpos_data[(k+1)*ndyna:(k+2)*ndyna,:]

                    # out._static_data.resize((p+1)*nstatic,0)
                    # #print( dpos_data[k*nstatic:(k+1)*nstatic,:])
                    # out._static_data[p*nstatic:(p+1)*nstatic,:] = spos_data[k*nstatic:(k+1)*nstatic,:]

                    out._solv_data.resize((p+1),0)
                    out._solv_data[p:(p+1),:] = numpy.array(solv_data[k:k+1,:])

                    id_f = numpy.where(abs(cf_data[:, 0] - times[k]) < time_step*1e-5)[0]
                    if len(id_f) == 0:
                        print('no contact data at time',times[k])
                    else:
                        #print('index of contact :', min(id_f), max(id_f))
                        new_line = current_line+max(id_f)-min(id_f)+1
                        #print('new_line',new_line)
                        #print('current_line',current_line)
                        #print('size of contact data', max(id_f)-min(id_f)+1)
                        out._cf_data.resize(new_line,0)
                        #print('fill out._cf_data indices', current_line, new_line-1)
                        out._cf_data[current_line:new_line,:] = numpy.array(cf_data[min(id_f):max(id_f)+1,:])
                        current_line=new_line
                        #print('current_line',current_line)

                p = p+1
        #print(dpos_data)

        print(out._dynamic_data.shape)
        print(out._static_data.shape)
        print(out.static_data().value)
        print(out._cf_data.shape)
        print(out._solv_data.shape)
