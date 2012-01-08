import matplotlib.pyplot as plt
import sys

location_dir = "."
if (len(sys.argv)>1):
        location_dir = sys.argv[1]
if location_dir[-1] != '/':
        location_dir = location_dir + '/'
s = [float(i) for i in open(location_dir+"moment_ratios.txt").read().split()]
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(range(0,len(s)),s)
ax.set_xlabel("Time")
ax.set_ylabel("Ratio of 2nd and 1st moments")
plt.show()
