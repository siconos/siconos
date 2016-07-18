#!/usr/bin/env python

from __future__ import print_function
from pylab import *
import os, time

MAX_THREADS = 12
SAMPLES = 10

times = [list() for _ in range(MAX_THREADS)]

# Progressive updates
ion()

def time_n(n):
    t = time.time()
    cmd = 'env OMP_NUM_THREADS=%d siconos chute.py 2>&1 >/dev/null'%(n)
    print(cmd)
    res = os.system(cmd)
    print(res)
    #time.sleep(randint(3)+1)
    elapsed = time.time() - t
    return elapsed

def update_graph():
    clf()
    for x,result in enumerate(times):
        if len(result)>0:
            plot([x+1]*len(result), result, 'xk')
    xlim(0,MAX_THREADS+1)
    xticks(arange(1,MAX_THREADS+1))
    xlabel('OMP_NUM_THREADS')
    ylabel('time (s)')
    ylim(0,ylim()[1]*1.1)
    draw()

def boxplot_graph():
    clf()
    boxplot(times)
    xlim(0,MAX_THREADS+1)
    xticks(arange(1,MAX_THREADS+1))
    xlabel('OMP_NUM_THREADS')
    ylabel('time (s)')
    ylim(0,ylim()[1]*1.1)

# Do measurements
total = SAMPLES*MAX_THREADS
left = 'unknown'
totaltime = 0
for j in range(SAMPLES):
    for i in range(MAX_THREADS):
        current = j*MAX_THREADS+i+1
        print(current,'/',total,'-- left:',left,'mins')
        t = time_n(i+1)
        totaltime += t
        left = (totaltime / current * total - totaltime)/60.0
        times[i].append(t)
        update_graph()
        pause(0.1)

# Display final results
boxplot_graph()
ioff()
show()
