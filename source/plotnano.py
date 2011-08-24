from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
#IF MATPLOTLIB version >= 1 then
#ax = fig.add_subplot(111, projection='3d')
#Else for version < 1
ax = Axes3D(fig)
for fname, c in [('nano.txt','b'),('plmr.txt','r')]:
	f = open(fname,'r')
	x=[];y=[];z=[]
	while 1:
		s = f.readline()
		if (s):
			data = s.split()
			x.append(float(data[0]));y.append(float(data[1]));z.append(float(data[2]));
		else:
			break

	ax.scatter(x,y,z,c='b',marker='o',color=c)
	f.close()

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()

