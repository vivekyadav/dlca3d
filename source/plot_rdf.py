import math
import os
import sys
import matplotlib
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
	print "\nUsage : "+sys.argv[0]+" [directory name]"
	exit(0)


fig = plt.figure()
ax = fig.add_subplot(111)

for i in os.listdir(sys.argv[1]):
    if 'rdf_' in i:
        n = 0
        sum_x = 0
        sum_y = 0
        sum_x_squared = 0
        sum_xy = 0

        f=open(sys.argv[1]+i,'r')
        vector_x=[];vector_y = []
        while 1:
            s=f.readline()
            if(s):
                data = s.split()
                if float(data[0]) == 0 or float(data[1]) == 0:
                    continue
                x = (float(data[0]))
                y = (float(data[1]))
                vector_x.append(x);vector_y.append(y)
                n += 1
                sum_x += x
                sum_y += y
                sum_x_squared += x*x
                sum_xy += x*y
            else:
                break

        slope = (sum_xy - sum_x*sum_y/n) / (sum_x_squared - sum_x**2/n)
        intercept = (sum_y - slope*sum_x) / n
        bestline_x = [vector_x[0],vector_x[-1]]
        bestline_y = [vector_x[0]*slope+intercept,vector_x[-1]*slope+intercept]

        ax.plot(vector_x,vector_y, '-')
        #ax.plot(bestline_x,bestline_y,'-')
#ax.set_title("BOx Length Vs Number of particles.\nline of best fit: slope= %f  intercept= %f" % (slope, intercept))
plt.show()
