import math
import sys
import matplotlib
import matplotlib.pyplot as plt

n = 0
sum_x = 0
sum_y = 0
sum_x_squared = 0
sum_xy = 0

if len(sys.argv) < 2:
	print "\nUsage : "+sys.argv[0]+" [input filename]"
	exit(0)
f=open(sys.argv[1],'r')
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
bestline_x = [vector_x[vector_x.index(min(vector_x))],vector_x[vector_x.index(max(vector_x))]]
bestline_y = [bestline_x[0]*slope+intercept,bestline_x[-1]*slope+intercept]

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(vector_x,vector_y, 'co')
ax.plot(bestline_x,bestline_y,'-')
#ax.set_yscale('log')
#ax.set_xscale('log')
ax.set_title("BOx Length Vs Number of particles.\nline of best fit: slope= %f  intercept= %f" % (slope, intercept))
plt.show()
