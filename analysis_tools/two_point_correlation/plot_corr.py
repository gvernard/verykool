import matplotlib.pyplot as plt
import numpy as np
import json


#plt.style.use('../plot_style')


f = open("options.json","r")
input_str = f.read()
pars = json.loads(input_str)


r,curv,modgauss,gauss = np.loadtxt("corr_theo.dat",unpack=True)
#curv /= np.amax(curv)

#modgauss *= 0.2
#gauss *= 0.2

r_true_merger,s_true_merger = np.loadtxt("corr_true_merger.dat",unpack=True)
s_true_merger /= np.amax(s_true_merger)
#s_true -= np.unique(s_true)[1]

r_true_m83,s_true_m83 = np.loadtxt("corr_true_m83.dat",unpack=True)
s_true_m83 /= np.amax(s_true_m83)


#r_model,s_model = np.loadtxt("corr_model.dat",unpack=True)
#s_model /= np.amax(s_model)



mycolor='lightgrey'
plt.plot(r,curv,drawstyle='default',color=mycolor,ls='-',label="curvature")
plt.plot(r,modgauss,drawstyle='default',color=mycolor,ls=':',label=r"modgauss $\sigma="+str(pars['modgauss'][1]['val'])+"$")
plt.plot(r,gauss,drawstyle='default',color=mycolor,ls='--',label=r"gauss $\sigma="+str(pars['gauss'][1]['val'])+"$")
#plt.plot(r,reconstruction,drawstyle='default',color='red',ls='-')

plt.plot(r_true_merger,s_true_merger,drawstyle='default',color='red',ls='-',label="merger true")
plt.plot(r_true_m83,s_true_m83,drawstyle='default',color='blue',ls='-',label="M83 true")

#plt.plot(r_model,s_model,drawstyle='default',color='red',ls='-',label="model")


#plt.plot(r,np.exp(-r/0.14427))
#plt.axvline(0.1)
#plt.axhline(0.5)

plt.xlim(0,1.25)
plt.ylim(0,1.1)
plt.xlabel('r')
plt.ylabel(r'$\xi(r)$')

plt.legend()


plt.savefig('corr.png')
