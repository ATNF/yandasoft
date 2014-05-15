# plot 44 GHz fluxes vs. 36 GHz ones
import numpy
import matplotlib
import matplotlib.pyplot as plt

data = numpy.loadtxt('result.dat',usecols=[0,1,2,3]).transpose()
indices = [i for i in range(len(data[0])) if data[3][i]>12]

x = numpy.arange(-90,60,1.5)
y = numpy.arange(-12.,12.,0.24)
el = matplotlib.mlab.griddata(data[0],data[1],data[3],x,y)  

data=data[:][:,indices]
z = matplotlib.mlab.griddata(data[0],data[1],data[2],x,y)  

plt.rc('legend',fontsize=10)

ax1 = plt.subplot(1,1,1)
#ax1.yaxis.set_label_coords(ylabelx,0.5)
cs = plt.contour(x,y,z, colors='k')
plt.contour(x,y,el,[12], colors='r',linewidths=3)
plt.clabel(cs, inline=1, fontsize=10,fmt='%i')
plt.ylabel('Hour angle, h')
plt.xlabel('Declination, deg')
#plt.title("Residual w-term for a single offset beam")
#plt.title("Residual w-term for 4 offset beams")
plt.title("Residual w-term for 36 beams")
#plt.xlim(-1,88)

#plt.gcf().subplots_adjust(hspace=0,wspace=0)
#plt.show()
#plt.savefig("plot.eps",transparent=True)
plt.savefig("plot.png",transparent=False)
