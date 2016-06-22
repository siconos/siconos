import numpy

#results = numpy.loadtxt('simulation_results.dat')

results = []
results_ref = []

with open ("simulation_results.dat", "r") as myfile:
    data=myfile.readlines()
    #print data[1:]
    for i in data[1:]:
        data_l_i = i.split('\t')
        data_l_i_f = []
        #for val in i:
        #    print val
        for val in data_l_i :
            try:
                float(val)
                data_l_i_f.append(float(val))
            except:
                print (' -- ')
        results =  numpy.array(data_l_i_f)
        #print i,float(i.split('\t'))
        #print ' '

with open ("simulation_results.ref", "r") as myfile:
    data=myfile.readlines()
    #print data[1:]
    for i in data[1:]:
        data_l_i = i.split('\t')
        data_l_i_f = []
        #for val in i:
        #    print val
        for val in data_l_i :
            try:
                float(val)
                data_l_i_f.append(float(val))
            except:
                print (' -- ')
        results_ref =  numpy.array(data_l_i_f)
        #print i,float(i.split('\t'))
        #print ' '

error_= numpy.linalg.norm(results-results_ref)
print 'Error w.r.t reference file = ', error_

assert (error_< 1e-14)
