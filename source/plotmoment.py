import matplotlib.pyplot as plt

s = [float(i) for i in open("moment_ratios.txt").read().split()]
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(range(0,len(s)),s)
ax.set_xlabel("Time")
ax.set_ylabel("Ratio of 2nd and 1st moments")
plt.show()
