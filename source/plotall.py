import math
import os
import matplotlib
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111)
lines = []

for i in os.listdir('.'):
    if 'result_100_0.000' in i:
        n = 0
        sum_x = 0
        sum_y = 0
        sum_x_squared = 0
        sum_xy = 0

        f=open(i+"/RvsN.txt",'r')
        vector_x=[];vector_y = []
        while 1:
            s=f.readline()
            if(s):
                data = s.split()
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
        l, =ax.plot(bestline_x,bestline_y,'-',label = str(slope)+",conc="+i[-8:-3])
        lines.append(l)
#ax.set_title("BOx Length Vs Number of particles.\nline of best fit: slope= %f  intercept= %f" % (slope, intercept))
ax.legend(lines,[l.get_label() for l in lines],loc='upper left')
plt.show()
