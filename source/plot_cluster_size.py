import matplotlib.pyplot as plt
import sys
import numpy as np
import matplotlib.mlab as mlab

if len(sys.argv) < 2:
	print "\nUsage : "+sys.argv[0]+" [input filename]"
	exit(0)

f = open(sys.argv[1],'r')
mass = []
while 1:
    s = f.readline()
    if s :
        data = s.split()
        if len(data) < 2:
            break
        mass.append(float(data[1]))
    else:
        break

fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(mass,50, facecolor='green', alpha=0.75)

ax.set_xlabel('Mass')
ax.set_ylabel('Count')
#ax.set_title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')

ax.grid(True)

plt.show()