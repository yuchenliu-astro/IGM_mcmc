"""
MCMC for IGM transmission
"""


import astropy.units as u
import numpy as np
from astropy.table import Table
from astropy.io import ascii
import datetime
import time


Wlen = []
Os_st = []
Gamma = []

T = []
with open('VPFIT/atom.dat', 'r') as file:
    for line in file:
        items = line.split()
        T.append(items)

LS = T[8:40]
for ls in LS:
    Wlen.append(float(ls[2]))
    Os_st.append(float(ls[3]))
    Gamma.append(float(ls[4]))

LS_tb = Table(np.c_[Wlen,Os_st,Gamma],names=('Wlen','Os_st','Gamma'))



#optical depth for LyC photons
def tau_lyc(NHI,wlen):
    nu = (3e10*u.cm/u.s / (wlen * u.angstrom)).to(u.Hz).value
    nu0 = (3e10*u.cm/u.s / (911.8 * u.angstrom)).to(u.Hz).value
    tau_array = np.zeros_like(wlen)
    tau_array[np.where(wlen>911.8)[0]] = 0
    tau_array[np.where(wlen<=911.8)[0]] = NHI * 6.3e-18 * (wlen[np.where(wlen<=911.8)]/911.8)**3
    if wlen>911.8 : 
        tau = 0
    else:
        tau = NHI * 6.3e-18 * (wlen/911.8)**3
    return tau_array




def line_profile(wlen,b,LS_tb,i):
    c = 3e5 #km/s
    b = b #km/s
    lambda_i = LS_tb['Wlen'][i] #A
    lambda_d = (b/c)*lambda_i #A
    x = (wlen - lambda_i)/lambda_d  #[0,10]
    a = lambda_i**2*LS_tb['Gamma'][i]/4/np.pi/(c * 10**13)/lambda_d #[1e-4,1e-8]
    K_array = np.zeros_like(x)
    K_array[np.where(x**2>700)[0]] = 0
    K_array[np.where(x**2<=700)[0]] = 1/(2*x[np.where(x**2<700)]**2)*((4*x[np.where(x**2<700)]**2+3)*(x[np.where(x**2<700)]**2+1) \
                                        *np.exp(-x[np.where(x[np.where(x**2<700)]**2<700)]**2)-1/x[np.where(x**2<700)]**2*(2*x[np.where(x**2<700)]**2+3)*np.sinh(x[np.where(x**2<700)]**2))
   # if x**2>700:
   #     K = 0
   # else:
   #     K = 1/(2*x**2)*((4*x**2+3)*(x**2+1)*np.exp(-x**2)-1/x**2*(2*x**2+3)*np.sinh(x**2))
    H_array = np.exp(-x**2)*(1-a*2/np.sqrt(np.pi)*K_array)
    return x,a,K_array,H_array




def tau_LS(NHI,wlen,b):
    m_e = 9.1e-28 #g
    c = 3e5 #km/s
    nu = c*1e13/wlen #Hz
    nu_D = c*1e13/LS_tb['Wlen']*(b/c)

    CS = []

    for ll in range(len(LS_tb)):
        H = line_profile(wlen,b,LS_tb,ll)[3]
        const_term = np.sqrt(np.pi) * 4.8e-10**2/m_e/3e10
        cross_secion = const_term * LS_tb['Os_st'][ll]/b* LS_tb['Wlen'][ll]*1e-13 * H
        CS.append(cross_secion)
    
    tauLS = np.sum(CS,axis = 0) * NHI
    return np.sum(CS,axis = 0),tauLS


def tau_abs(abs):
    wlenrest_reg = np.array([700,1300])
    wlen_eachabs_reg = wlenrest_reg * (1+abs['z'])
    wlen_eachabs = np.linspace(int(wlen_eachabs_reg[0]),int(wlen_eachabs_reg[1]),(int(wlen_eachabs_reg[1])-int(wlen_eachabs_reg[0]))*10+1)
    wlenrest = wlen_eachabs/(1+abs['z'])
    tau_eachlyc = np.array([tau_lyc(abs['Nhi'],wlen) for wlen in wlenrest])
    tau_eachLS = np.array([tau_LS(abs['Nhi'],wlen,abs['b'])[1] for wlen in wlenrest])
    tau_eachabs = np.array([tau_lyc(abs['Nhi'],wlen) + tau_LS(abs['Nhi'],wlen,abs['b'])[1] for wlen in wlenrest])
    tau_eachabs = tau_lyc(abs['Nhi'],wlenrest) + tau_LS(abs['Nhi'],wlenrest,abs['b'])[1]

    return wlen_eachabs, tau_eachabs



def interpp(tau_eachabs):
    wlen = tau_eachabs[:,0]
    tau = tau_eachabs[:,1]
    new_wlen = np.round(np.linspace(int(wlen[0]),int(wlen[-1]),num = int((int(wlen[-1])-int(wlen[0]))/0.1 + 1)),decimals=1)
    tau_inter = np.interp(new_wlen,tau_eachabs[:,0],tau_eachabs[:,1])

    return new_wlen,tau_inter



def essemble_fun(absorber,path,filename):
    wlen_essemble = np.round(np.linspace(0,6500,65001),decimals=1)
    tau_essemble = np.zeros_like(wlen_essemble)
   # tau_lyc_essemble = np.zeros_like(wlen_essemble)
   # tau_LS_essemble = np.zeros_like(wlen_essemble)

    for i in range(len(absorber)):
        wlen_eachabs,tau_eachabs = tau_abs(absorber[i])
        wlen_inter,tau_inter = interpp(tau_eachabs)
        wlen_inter,tau_lyc_inter = interpp(tau_eachlyc)
        wlen_inter,tau_LS_inter = interpp(tau_eachLS)
        match_wlen = np.intersect1d(wlen_eachabs,wlen_essemble)
        essemble_position = np.where(np.isin(wlen_essemble,match_wlen))
        inter_position = np.where(np.isin(wlen_eachabs,match_wlen))
        tau_essemble[int(wlen_eachabs[0] * 10):(int(wlen_eachabs[-1] * 10))+1] += tau_eachabs



    wlen_essemble_cut = wlen_essemble[20000:60001]
    tau_essemble_cut = tau_essemble[20000:60001]
    ascii.write(np.c_[wlen_essemble_cut,tau_essemble_cut],path+filename,overwrite=True)
    #print(datetime.datetime.now())
    tau_essemble[essemble_position[0]] += tau_eachabs[inter_position[0]]
    tau_lyc_essemble[essemble_position[0]] += tau_lyc_inter[inter_position[0]]
    tau_LS_essemble[essemble_position[0]] += tau_LS_inter[inter_position[0]]
    return wlen_essemble_cut, tau_essemble_cut 
