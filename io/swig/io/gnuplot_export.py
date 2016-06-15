import os,sys
import h5py
import numpy

filename = '{0}.hdf5'.format(os.path.splitext(os.path.basename(sys.argv[1]))[0])

withPlot=False


print filename
out= h5py.File(filename, 'r')
def group(h, name):
    try:
        return h[name]
    except KeyError:
        return h.create_group(name)

def data(h, name, nbcolumns):
    try:
        return h[name]
    except KeyError:
        return h.create_dataset(name, (0, nbcolumns),
                                maxshape=(None, nbcolumns))

    
_data = group(out, 'data')
ref = group(_data, 'ref')
joints = group(_data, 'joints')
static_data = data(_data, 'static', 9)
velocities_data = data(_data, 'velocities', 8)
dynamic_data = data(_data, 'dynamic', 9)
cf_data = data(_data, 'cf', 15)
solv_data = data(_data, 'solv', 4)
input = group(_data, 'input')
nslaws = group(_data, 'nslaws')

dpos_data = dynamic_data
max_time = max(dpos_data[:, 0])
times = list(set(dpos_data[:, 0]))
times.sort()
ndyna = len(numpy.where(dpos_data[:, 0] == times[0]))
ntime=len(times)


print('time range :', times[0], times[-1])
print('ndyna :', ndyna)
print('ntime:', ntime)
instances = set(dpos_data[:, 1])

#output_dict = {}

#output_dict[1]= [1,2,3]


######## position output ########

nvalue = ndyna*7+1

position_output = numpy.empty((ntime,nvalue))
#print('position_output shape', numpy.shape(position_output))
position_output[:,0] = times[:]
for t in range(len(times)):
    for i in range(ndyna):
        position_output[t,1+i*7:1+(1+i)*7] = dpos_data[t*ndyna+ndyna, 2:9]
#print('position_output', position_output)
filename_output = '{0}_position.dat'.format(os.path.splitext(os.path.basename(sys.argv[1]))[0])
print('output file:', filename_output)
numpy.savetxt(filename_output, position_output)

######## position output ########
nvalue = ndyna*6+1

velocity_output = numpy.empty((ntime,nvalue))
#print('position_output shape', numpy.shape(position_output))
velocity_output[:,0] = times[:]
for t in range(len(times)):
    for i in range(ndyna):
        velocity_output[t,1+i*6:1+(1+i)*6] = velocities_data[t*ndyna+ndyna, 2:8]
#print('position_output', position_output)
filename_output = '{0}_velocity.dat'.format(os.path.splitext(os.path.basename(sys.argv[1]))[0])
print('output file:', filename_output)
numpy.savetxt(filename_output, velocity_output)



if withPlot:
    import matplotlib
    havedisplay = "DISPLAY" in os.environ
    if not havedisplay:
        matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    plt.subplot(411)
    plt.title('position x')
    plt.plot(position_output[:, 0], position_output[:, 1])
    plt.subplot(412)
    plt.title('position y')
    plt.plot(position_output[:, 0], position_output[:, 2])
    plt.subplot(413)
    plt.title('position z ')
    plt.plot(position_output[:, 0], position_output[:, 3])


    plt.figure()
    plt.subplot(411)
    plt.title('orientation q0')
    plt.plot(position_output[:, 0], position_output[:, 4])
    plt.subplot(412)
    plt.title('orientation q1')
    plt.plot(position_output[:, 0], position_output[:, 5])
    plt.subplot(413)
    plt.title('orientation q2 ')
    plt.plot(position_output[:, 0], position_output[:, 6])
    plt.subplot(414)
    plt.title('orientation q3 ')
    plt.plot(position_output[:, 0], position_output[:, 7])



    if havedisplay:
        plt.show()
    else:
        plt.savefig("bbts.png")
