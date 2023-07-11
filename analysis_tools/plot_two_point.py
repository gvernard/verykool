import matplotlib.pyplot as plt
import numpy as np
from scipy.special import kv,gamma


def myexp(r,l):
    return np.exp(-r/l)

def mygauss(r,l):
    return np.exp(-r*r/(2*l*l))

def mymatern(d,e,l):
    arg = np.sqrt(2*e)*d/l
    K = kv(e,arg)
    G = gamma(e)
    return pow(2,1.0-e)*pow(arg,e)*K/G




clength = 0.2
eta = 0.01
r = np.arange(0.01,1.8,0.01)

exp = np.zeros((len(r)))
gauss = np.zeros((len(r)))
matern = np.zeros((len(r)))
for i in range(0,len(r)):
    exp[i] = myexp(r[i],clength)
    gauss[i] = mygauss(r[i],clength)
    matern[i] = mymatern(r[i],eta,clength)



fig,sub = plt.subplots(1,figsize=(8,5))
art_fit_exp,    = plt.plot(r,exp,drawstyle='default',color='red',ls='-')
art_fit_gauss,  = plt.plot(r,gauss,drawstyle='default',color='blue',ls='-')
art_fit_matern, = plt.plot(r,matern,drawstyle='default',color='black',ls='--')


plt.xlim(0,0.8)
plt.ylim(0,1.1)
plt.xlabel('r [arcsec]')
plt.ylabel(r'$\xi(r)$')
plt.setp(sub.get_xticklabels()[0], visible=False)

plt.savefig('corr.pdf')


