import matplotlib.pyplot as plt
import numpy as np
import json

f = open("options.json","r")
input_str = f.read()
pars = json.loads(input_str)


r,curv,modgauss,gauss = np.loadtxt("corr_theo.dat",unpack=True)
curv /= np.amax(curv)

r_true,s_true = np.loadtxt("corr_true.dat",unpack=True)
s_true /= np.amax(s_true)
#s_true -= np.unique(s_true)[1]


r_model,s_model = np.loadtxt("corr_model.dat",unpack=True)
s_model /= np.amax(s_model)



mycolor='lightgrey'
plt.plot(r,curv,drawstyle='steps-post',color=mycolor,ls='-',label="curvature")
plt.plot(r,modgauss,drawstyle='steps-post',color=mycolor,ls=':',label=r"modgauss $\sigma="+str(pars['modgauss'][1]['val'])+"$")
plt.plot(r,gauss,drawstyle='steps-post',color=mycolor,ls='--',label=r"gauss $\sigma="+str(pars['gauss'][1]['val'])+"$")
#plt.plot(r,reconstruction,drawstyle='steps-post',color='red',ls='-')

plt.plot(r_true,s_true,drawstyle='steps-post',color='black',ls='-',label="true")

plt.plot(r_model,s_model,drawstyle='steps-post',color='red',ls='-',label="model")


#plt.plot(r,np.exp(-r/0.14427))
#plt.axvline(0.1)
#plt.axhline(0.5)

plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('r')
plt.ylabel(r'$\xi(r)$')

plt.legend()


plt.savefig('corr.png')
