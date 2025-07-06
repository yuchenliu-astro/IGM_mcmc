"""
mcmc of IGM 
"""


from multiprocessing import Pool
#import transition
from astropy.io import ascii
from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
import transmission
import time
import datetime

absorber_all = []
for los in range(500):
    absorber = ascii.read('IGM_MC_z3.5/input/absorber_'+str(los)+'.dat')
    absorber_all.append(absorber)



path = 'IGM_MC_z3.5/result_dat/'
filename_head = 'absorber_'
result_nm_list = ['tau_tb_LoS'+filename_head+str(losnum)+'full.dat' for losnum in range(500)]

print(type(absorber_all[0]),type(result_nm_list[0]),type(path))
if __name__ == '__main__':
    pool = Pool(50)
    input = [absorber_all,[path]*500,result_nm_list]
    res = pool.starmap_async(transmission.essemble_fun,input)
    print(datetime.datetime.now())
    pool.close()
    pool.join()

print(datetime.datetime.now())


'''
# for test
for i in range(1000):
    absorber = ascii.read('IGM_MC/absorber_'+str(i)+'.dat')
    if __name__ == '__main__':
        with Pool(64) as p:
            input = zip(absorber)
            tau_oneLoS = p.starmap(transmission.essemble_fun,input)

    tau_abs_oneLoS = [tau_oneLoS[i][0] for i in range(len(tau_oneLoS))]
    tau_lyc_oneLoS = [tau_oneLoS[i][1] for i in range(len(tau_oneLoS))]
    tau_LS_oneLoS = [tau_oneLoS[i][2] for i in range(len(tau_oneLoS))]

    wlen_essemble = np.round(np.linspace(700,5000,43001),decimals=1)
    plt.figure(figsize=[20,5])
#for tau_one in tau_lyc_oneLoS[-10:-1]:
#    plt.plot(wlen_essemble/4.1,np.exp(-tau_one),'k')
    tau_essemble = np.sum(np.array(tau_abs_oneLoS),axis=0)
    tau_ls_essemble = np.sum(np.array(tau_LS_oneLoS),axis=0)
    tau_lyc_essemble = np.sum(np.array(tau_lyc_oneLoS),axis=0)
#plt.plot(wlen_essemble/4.1,np.exp(-tau_ls_essemble))
#plt.plot(wlen_essemble/4.1,np.exp(-tau_lyc_essemble))
    table_oneLoS = Table(np.c_[wlen_essemble,-tau_essemble],names=('wlen_obs','tau'))
    ascii.write(table_oneLoS,'IGM_MC/result_dat/tau_tb_LoS'+str(i)+'.dat',overwrite=True)

    plt.plot(wlen_essemble/4.1,np.exp(-tau_essemble),'k',alpha=0.8)
    plt.xlim(700,1300)
    plt.xlabel(r'$\lambda_{rest}$',fontsize=18)
    plt.ylabel(r'$Transmission$',fontsize=18)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig('IGM_MC/result_img/transmission_LoS'+str(i)+'.jpg',dpi=500,bbox_inches='tight',pad_inches = 0.3)
    '''
