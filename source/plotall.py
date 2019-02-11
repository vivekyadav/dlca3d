import math
import os
import matplotlib.pyplot as plt
import sys
import re

fig = plt.figure()
ax = fig.add_subplot(111)
lines = []
location_dir = "."
if (len(sys.argv)>1):
        location_dir = sys.argv[1]
if location_dir[-1] != '/':
        location_dir = location_dir + '/'
for i in os.listdir(location_dir):
    if 'result' in i:
        n = 0
        sum_x = 0
        sum_y = 0
        sum_x_squared = 0
        sum_xy = 0

        f=open(location_dir+i+"/RvsN.txt",'r')
        vector_x=[];vector_y = []
        while 1:
            s=f.readline()
            if(s):
                data = s.split()
                if len(data) != 2:
                    continue
                if float(data[0])<=0 or float(data[1])< 10:
                    continue
                x = math.log(float(data[0]))
                y = math.log(float(data[1]))
                if y/x < 2:
                    continue
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
        bestline_x = [0,10]
        bestline_y = [intercept,10*slope+intercept]

        #ax.plot(vector_x,vector_y, 'o')
        m = re.match(r'result_(\d+)_([0-9.]*)_([0-9.]*)_(\d+)',i)
        l, =ax.plot(bestline_x,bestline_y,'-',label = str(slope)+",conc="+m.group(3))
        lines.append(l)
#ax.set_title("BOx Length Vs Number of particles.\nline of best fit: slope= %f  intercept= %f" % (slope, intercept))
ax.legend(lines,[l.get_label() for l in lines],loc='upper left')
plt.show()
